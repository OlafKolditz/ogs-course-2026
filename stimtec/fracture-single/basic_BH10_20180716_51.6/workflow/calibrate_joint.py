"""
Gemeinsame 2-Parameter-Kalibrierung von Permeabilitaet (kappa_temporal) und
Storage (storage_temporal) in simulation/stimtec-square-expression.prj,
sodass die an x,y=(0,0) simulierte Druckaufbaukurve moeglichst gut mit der
digitalisierten, gemessenen BH10-Kurve uebereinstimmt.

Beide Funktionen behalten ihre Kurvenform (Zeitfenster/Stufen, Verhaeltnis
max/min) bei; es werden zwei globale Skalierungsfaktoren S_k (Permeabilitaet)
und S_s (Storage) auf die urspruenglichen Basiswerte angewendet:
    k_min(S_k)       = S_k * 5e-11,   k_max(S_k)       = S_k * 1e-9
    storage_min(S_s) = S_s * 1e-7,    storage_max(S_s) = S_s * 1e-6

Fuer jeden Kandidaten (S_k, S_s) wird:
  1. sowohl der kappa_temporal- als auch der storage_temporal-Block in der
     .prj-Datei neu geschrieben (ein Dateizugriff)
  2. OGS ausgefuehrt
  3. der Druck an (0,0) aus den Ergebnis-VTUs extrahiert
  4. der RMSE gegenueber der (basislinienkorrigierten) Messung berechnet

Ein 2D-Optimierer (scipy.optimize.minimize, Nelder-Mead mit Bounds) sucht
das RMSE-minimale Paar (S_k, S_s). Am Ende wird die .prj-Datei mit den besten
gefundenen Faktoren final beschrieben.
"""

import os
import re
import subprocess
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd
import pyvista as pv
from scipy.optimize import minimize

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
SIM_DIR = os.path.join(SCRIPT_DIR, "..", "simulation")
PRJ_PATH = os.path.join(SIM_DIR, "stimtec-square-expression.prj")
RESULTS_DIR = os.path.join(SIM_DIR, "..", "results")
PVD_PATH = os.path.join(RESULTS_DIR, "stimtec.pvd")
MEAS_CSV = os.path.join(SCRIPT_DIR, "bh10_gemessen_digitalisiert.csv")
OGS_EXE = r"C:\Users\okolditz\AppData\Local\Programs\Python\Python313\Scripts\ogs.exe"

# Urspruengliche Basiswerte (vor jeglicher Kalibrierung), auf die S_k/S_s
# angewendet werden -- unabhaengig vom aktuellen Stand der .prj-Datei.
BASE_K_MIN = 5e-11
BASE_K_MAX = 1e-9
BASE_S_MIN = 1e-7
BASE_S_MAX = 1e-6

V_MIN, V_MAX = 0.001, 10.0

# Zeitpunkte/Werte der Stufenfunktion (identisch zu pressure_source_term)
COORDS = [2500, 2650, 2660, 2800, 2810, 2950, 2960, 3070, 3080, 3200, 3210,
          3300, 3310, 3450, 3460, 3580, 3590, 3650, 3660, 3730, 3740, 3810,
          3820, 3900, 3910, 4000]
VALUES = [0.001, 0.001, 0.6, 0.6, 1.2, 1.2, 1.7, 1.7, 2.0, 2.0, 3.0, 3.0,
          4.0, 4.0, 5.0, 5.0, 6.0, 6.0, 7.0, 7.0, 8.0, 8.0, 9.0, 9.0, 10.0, 10.0]

WERT_SKALIERUNG = 1.0e6  # angenommene Einheit der Messung: MPa -> Pa

# Auswertungspunkt
POINT = (0.0, 0.0)


def build_expression(v_min_target, v_max_target):
    """Baut die verschachtelte ?:-Expression fuer einen Parameter mit
    denselben Zeitfenstern wie pressure_source_term, linear abgebildet
    zwischen v_min_target und v_max_target."""
    values = [
        v_min_target + (v - V_MIN) / (V_MAX - V_MIN) * (v_max_target - v_min_target)
        for v in VALUES
    ]

    lines = []
    n = len(COORDS)
    for i in range(n - 1):
        t1, t2 = COORDS[i], COORDS[i + 1]
        p1, p2 = values[i], values[i + 1]
        if i == 0:
            lines.append(f"(t &lt; {t2} ? {p1:.6e} :")
        elif p1 == p2:
            lines.append(f" t &lt; {t2} ? {p1:.6e} :")
        else:
            lines.append(
                f" t &lt; {t2} ? {p1:.6e} + ({p2:.6e} - {p1:.6e}) / ({t2} - {t1}) * (t - {t1}) :"
            )
    lines.append(f" {values[-1]:.6e})")
    return "\n                ".join(lines)


def write_parameter_block(name, comment, v_min_target, v_max_target, content):
    """Ersetzt den <parameter name="{name}"> Block im uebergebenen Dateiinhalt
    und gibt den neuen Inhalt zurueck (schreibt noch nicht auf Platte)."""
    expression = build_expression(v_min_target, v_max_target)

    new_block = f"""<parameter>
            <name>{name}</name>
            <type>Function</type>
            <!-- {comment} -->
            <expression>
                {expression}
            </expression>
        </parameter>"""

    pattern = re.compile(
        rf"<parameter>\s*<name>{name}</name>.*?</parameter>", re.DOTALL
    )
    content_new, count = pattern.subn(new_block, content)
    if count != 1:
        raise RuntimeError(f"Erwartet genau 1 {name}-Block, gefunden: {count}")
    return content_new


def write_both_blocks(k_min, k_max, s_min, s_max):
    with open(PRJ_PATH, "r", encoding="ISO-8859-1") as f:
        content = f.read()

    content = write_parameter_block(
        "kappa_temporal",
        f"Gemeinsam kalibriert: k_min={k_min:.6e}, k_max={k_max:.6e} "
        f"(S_k={k_min/BASE_K_MIN:.6g} auf Basis {BASE_K_MIN:.3e}/{BASE_K_MAX:.3e})",
        k_min, k_max, content,
    )
    content = write_parameter_block(
        "storage_temporal",
        f"Gemeinsam kalibriert: storage_min={s_min:.6e}, storage_max={s_max:.6e} "
        f"(S_s={s_min/BASE_S_MIN:.6g} auf Basis {BASE_S_MIN:.3e}/{BASE_S_MAX:.3e})",
        s_min, s_max, content,
    )

    with open(PRJ_PATH, "w", encoding="ISO-8859-1") as f:
        f.write(content)


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


def objective(x):
    log10_Sk, log10_Ss = x
    S_k = 10 ** log10_Sk
    S_s = 10 ** log10_Ss
    key = (round(S_k, 6), round(S_s, 6))
    if key in _eval_cache:
        return _eval_cache[key]

    k_min = S_k * BASE_K_MIN
    k_max = S_k * BASE_K_MAX
    s_min = S_s * BASE_S_MIN
    s_max = S_s * BASE_S_MAX

    write_both_blocks(k_min, k_max, s_min, s_max)
    run_ogs()
    sim_t, sim_p = extract_pressure_at_point()
    meas_t, meas_p = load_measured_pa()
    err = rmse(sim_t, sim_p, meas_t, meas_p)

    _eval_cache[key] = err
    _eval_log.append((S_k, S_s, k_min, k_max, s_min, s_max, err))
    print(f"S_k={S_k:8.4f}  S_s={S_s:8.4f}  "
          f"k=[{k_min:.3e},{k_max:.3e}]  storage=[{s_min:.3e},{s_max:.3e}]  RMSE={err:.4e} Pa")
    return err


def main():
    print("Starte gemeinsame Kalibrierung von S_k (Permeabilitaet) und S_s (Storage) ...")
    print(f"Basis Permeabilitaet: k_min={BASE_K_MIN:.3e}, k_max={BASE_K_MAX:.3e} (S_k=1)")
    print(f"Basis Storage:        s_min={BASE_S_MIN:.3e}, s_max={BASE_S_MAX:.3e} (S_s=1)\n")

    # Startpunkt: zuvor einzeln kalibrierte Werte (S_k=1.3, S_s=3.43) als
    # guter Ausgangspunkt fuer die gemeinsame Suche.
    x0 = [np.log10(1.3), np.log10(3.43)]
    bounds = [(-1.0, 2.0), (-2.0, 2.0)]

    res = minimize(
        objective, x0, method="Nelder-Mead", bounds=bounds,
        options={"xatol": 0.02, "fatol": 1e3, "maxiter": 60, "maxfev": 60},
    )

    best_log10_Sk, best_log10_Ss = res.x
    best_Sk = 10 ** best_log10_Sk
    best_Ss = 10 ** best_log10_Ss
    best_k_min = best_Sk * BASE_K_MIN
    best_k_max = best_Sk * BASE_K_MAX
    best_s_min = best_Ss * BASE_S_MIN
    best_s_max = best_Ss * BASE_S_MAX

    print("\n--- Kalibrierung abgeschlossen ---")
    print(f"Bestes S_k = {best_Sk:.4f}  (k_min={best_k_min:.4e}, k_max={best_k_max:.4e})")
    print(f"Bestes S_s = {best_Ss:.4f}  (storage_min={best_s_min:.4e}, storage_max={best_s_max:.4e})")
    print(f"RMSE  = {res.fun:.4e} Pa")
    print(f"Funktionsauswertungen: {len(_eval_log)}")

    print("\nAlle getesteten Kandidaten (sortiert nach RMSE):")
    for S_k, S_s, k_min, k_max, s_min, s_max, err in sorted(_eval_log, key=lambda e: e[-1])[:15]:
        print(f"  S_k={S_k:7.3f}  S_s={S_s:7.3f}  RMSE={err:.4e}")

    # .prj-Datei final mit den besten Faktoren beschreiben und ein letztes
    # Mal rechnen, damit die Ergebnis-VTUs zum finalen Parametersatz passen.
    write_both_blocks(best_k_min, best_k_max, best_s_min, best_s_max)
    run_ogs()
    print(f"\n.prj-Datei final aktualisiert: {PRJ_PATH}")
    print(f"Ergebnisse (results/) entsprechen jetzt S_k={best_Sk:.4f}, S_s={best_Ss:.4f}.")


if __name__ == "__main__":
    main()
