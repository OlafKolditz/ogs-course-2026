"""
Kalibriert einen einzelnen Skalierungsfaktor S fuer die Storage-Funktion
storage_temporal in simulation/stimtec-square-expression.prj, sodass die an
x,y=(0,0) simulierte Druckaufbaukurve moeglichst gut mit der digitalisierten,
gemessenen BH10-Kurve uebereinstimmt.

Die Kurvenform (Zeitfenster/Stufen, Verhaeltnis storage_max/storage_min)
bleibt erhalten; es wird nur ein globaler Faktor S auf die Basiswerte
1e-7 und 1e-6 angewendet:
    storage_min(S) = S * 1e-7
    storage_max(S) = S * 1e-6

Die Permeabilitaet (kappa_temporal) bleibt unveraendert auf dem zuvor
kalibrierten Stand.

Fuer jeden Kandidaten S wird:
  1. der <parameter name="storage_temporal"> Block in der .prj-Datei neu geschrieben
  2. OGS ausgefuehrt
  3. der Druck an (0,0) aus den Ergebnis-VTUs extrahiert
  4. der RMSE gegenueber der (basislinienkorrigierten) Messung berechnet

Ein 1D-Optimierer (scipy.optimize.minimize_scalar) sucht den RMSE-minimalen
Faktor S. Am Ende wird die .prj-Datei mit dem besten gefundenen S final
beschrieben.
"""

import os
import re
import subprocess
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd
import pyvista as pv
from scipy.optimize import minimize_scalar

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
SIM_DIR = os.path.join(SCRIPT_DIR, "..", "simulation")
PRJ_PATH = os.path.join(SIM_DIR, "stimtec-square-expression.prj")
RESULTS_DIR = os.path.join(SIM_DIR, "..", "results")
PVD_PATH = os.path.join(RESULTS_DIR, "stimtec.pvd")
MEAS_CSV = os.path.join(SCRIPT_DIR, "bh10_gemessen_digitalisiert.csv")
OGS_EXE = r"C:\Users\okolditz\AppData\Local\Programs\Python\Python313\Scripts\ogs.exe"

# Basis-Storage-Werte (aktueller Stand vor der Kalibrierung)
BASE_S_MIN = 1e-7
BASE_S_MAX = 1e-6

V_MIN, V_MAX = 0.001, 10.0

# Zeitpunkte/Werte der Stufenfunktion (identisch zu pressure_source_term/kappa_temporal)
COORDS = [2500, 2650, 2660, 2800, 2810, 2950, 2960, 3070, 3080, 3200, 3210,
          3300, 3310, 3450, 3460, 3580, 3590, 3650, 3660, 3730, 3740, 3810,
          3820, 3900, 3910, 4000]
VALUES = [0.001, 0.001, 0.6, 0.6, 1.2, 1.2, 1.7, 1.7, 2.0, 2.0, 3.0, 3.0,
          4.0, 4.0, 5.0, 5.0, 6.0, 6.0, 7.0, 7.0, 8.0, 8.0, 9.0, 9.0, 10.0, 10.0]

WERT_SKALIERUNG = 1.0e6  # angenommene Einheit der Messung: MPa -> Pa

# Auswertungspunkt
POINT = (0.0, 0.0)


def build_storage_expression(s_min, s_max):
    """Baut die verschachtelte ?:-Expression fuer storage_temporal fuer gegebene s_min/s_max."""
    s_values = [
        s_min + (v - V_MIN) / (V_MAX - V_MIN) * (s_max - s_min) for v in VALUES
    ]

    lines = []
    n = len(COORDS)
    for i in range(n - 1):
        t1, t2 = COORDS[i], COORDS[i + 1]
        s1, s2 = s_values[i], s_values[i + 1]
        if i == 0:
            lines.append(f"(t &lt; {t2} ? {s1:.6e} :")
        elif s1 == s2:
            lines.append(f" t &lt; {t2} ? {s1:.6e} :")
        else:
            lines.append(
                f" t &lt; {t2} ? {s1:.6e} + ({s2:.6e} - {s1:.6e}) / ({t2} - {t1}) * (t - {t1}) :"
            )
    lines.append(f" {s_values[-1]:.6e})")
    return "\n                ".join(lines)


def write_storage_block(s_min, s_max):
    with open(PRJ_PATH, "r", encoding="ISO-8859-1") as f:
        content = f.read()

    expression = build_storage_expression(s_min, s_max)
    factor = s_min / BASE_S_MIN

    new_block = f"""<parameter>
            <name>storage_temporal</name>
            <type>Function</type>
            <!-- Kalibrierter Skalierungsfaktor S={factor:.6g} auf Basis-Stufenfunktion
                 (BASE_S_MIN={BASE_S_MIN:.3e}, BASE_S_MAX={BASE_S_MAX:.3e}):
                 storage_min = {s_min:.6e}, storage_max = {s_max:.6e} -->
            <expression>
                {expression}
            </expression>
        </parameter>"""

    pattern = re.compile(
        r"<parameter>\s*<name>storage_temporal</name>.*?</parameter>", re.DOTALL
    )
    content_new, count = pattern.subn(new_block, content)
    if count != 1:
        raise RuntimeError(f"Erwartet genau 1 storage_temporal-Block, gefunden: {count}")

    with open(PRJ_PATH, "w", encoding="ISO-8859-1") as f:
        f.write(content_new)


def run_ogs():
    result = subprocess.run(
        [OGS_EXE, os.path.basename(PRJ_PATH)],
        cwd=SIM_DIR,
        capture_output=True,
        text=True,
        timeout=180,
    )
    if result.returncode != 0:
        raise RuntimeError(f"OGS fehlgeschlagen (exit {result.returncode}):\n{result.stdout[-2000:]}\n{result.stderr[-2000:]}")


def read_pvd(pvd_path):
    tree = ET.parse(pvd_path)
    root = tree.getroot()
    base_dir = os.path.dirname(pvd_path)
    entries = []
    for ds in root.iter("DataSet"):
        t = float(ds.attrib["timestep"])
        fname = ds.attrib["file"]
        entries.append((t, os.path.join(base_dir, fname)))
    entries.sort(key=lambda e: e[0])
    return entries


def extract_pressure_at_point():
    entries = read_pvd(PVD_PATH)
    probe = pv.PolyData(np.array([[POINT[0], POINT[1], 0.0]]))
    times = []
    pressures = []
    for t, vtu_path in entries:
        mesh = pv.read(vtu_path)
        sampled = probe.sample(mesh)
        times.append(t)
        pressures.append(sampled["pressure"][0])
    return np.array(times), np.array(pressures)


def load_measured_pa():
    meas = pd.read_csv(MEAS_CSV)
    baseline = meas.loc[meas["time_s"] < 2650, "value_axis_units"].median()
    meas["pressure_Pa"] = (meas["value_axis_units"] - baseline) * WERT_SKALIERUNG
    return meas["time_s"].values, meas["pressure_Pa"].values


def rmse(sim_t, sim_p, meas_t, meas_p):
    meas_interp = np.interp(sim_t, meas_t, meas_p)
    return float(np.sqrt(np.mean((sim_p - meas_interp) ** 2)))


_eval_cache = {}
_eval_log = []


def objective(log10_S):
    S = 10 ** log10_S
    S_rounded = round(S, 6)
    if S_rounded in _eval_cache:
        return _eval_cache[S_rounded]

    s_min = S * BASE_S_MIN
    s_max = S * BASE_S_MAX
    write_storage_block(s_min, s_max)
    run_ogs()
    sim_t, sim_p = extract_pressure_at_point()
    meas_t, meas_p = load_measured_pa()
    err = rmse(sim_t, sim_p, meas_t, meas_p)

    _eval_cache[S_rounded] = err
    _eval_log.append((S, s_min, s_max, err))
    print(f"S={S:8.4f}  storage_min={s_min:.3e}  storage_max={s_max:.3e}  RMSE={err:.4e} Pa")
    return err


def main():
    print("Starte Kalibrierung des Storage-Skalierungsfaktors S ...")
    print(f"Basis: storage_min={BASE_S_MIN:.3e}, storage_max={BASE_S_MAX:.3e} (S=1)\n")

    res = minimize_scalar(
        objective, bounds=(-2.0, 2.0), method="bounded",
        options={"xatol": 0.03, "maxiter": 25},
    )

    best_log10_S = res.x
    best_S = 10 ** best_log10_S
    best_s_min = best_S * BASE_S_MIN
    best_s_max = best_S * BASE_S_MAX

    print("\n--- Kalibrierung abgeschlossen ---")
    print(f"Bestes S = {best_S:.4f}")
    print(f"storage_min = {best_s_min:.4e}, storage_max = {best_s_max:.4e}")
    print(f"RMSE  = {res.fun:.4e} Pa")

    print("\nAlle getesteten Kandidaten:")
    for S, s_min, s_max, err in sorted(_eval_log, key=lambda e: e[0]):
        print(f"  S={S:8.4f}  storage_min={s_min:.3e}  storage_max={s_max:.3e}  RMSE={err:.4e}")

    # .prj-Datei final mit dem besten Faktor beschreiben und ein letztes Mal rechnen,
    # damit die Ergebnis-VTUs zum finalen Parametersatz passen.
    write_storage_block(best_s_min, best_s_max)
    run_ogs()
    print(f"\n.prj-Datei final aktualisiert: {PRJ_PATH}")
    print(f"Ergebnisse (results/) entsprechen jetzt S={best_S:.4f}.")


if __name__ == "__main__":
    main()
