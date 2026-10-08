"""
Gemeinsame Kalibrierung der konstanten Medium-Eigenschaften
<permeability> und <storage> in simulation/stimtec.prj, so dass der
simulierte Druck (process_variable "pressure") am Brunnen (0,0) bestmoeglich
mit der gemessenen BH10-Kurve (2. Spalte in bh10_gemessen_digitalisiert.csv)
uebereinstimmt.

Messung: Spalte 2 in MPa, Baseline (Median t<2650 s) wird abgezogen, -> Pa.
Zielfunktion: RMSE [Pa] zwischen Simulation und auf die Simulationszeiten
interpolierter Messung (1750..2670 s, Pumpversuch 40.6).
Optimierer: Log-Raster (Vorsuche) + Nelder-Mead in (log10 k, log10 S).
"""
import bh40  # 40.6: Messung/Baseline
import os, re, subprocess, sys, shutil, xml.etree.ElementTree as ET
import numpy as np, pandas as pd, pyvista as pv
from scipy.optimize import minimize

HERE = os.path.dirname(os.path.abspath(__file__))
SIM_DIR = os.path.join(HERE, "..", "simulation")
PRJ_NAME = os.environ.get("PRJ_NAME", "stimtec.prj")   # z.B. stimtec_bc_right.prj
PRJ = os.path.join(SIM_DIR, PRJ_NAME)
TAG = os.path.splitext(PRJ_NAME)[0]
RES = os.path.join(HERE, "..", "results")
PVD = os.path.join(RES, TAG + ".pvd")
MEAS = os.path.join(HERE, "bh10_40_gemessen_digitalisiert.csv")
OGS = shutil.which("ogs") or "ogs"
SCALE = 1e6
WELL = np.array([0.0, 0.0, 0.0])
P_REF = float(os.environ.get("P_REF", "0"))  # Referenzdruck (Anfangsdruck), wird vom Simulationsdruck abgezogen

with open(PRJ, encoding="ISO-8859-1") as f:
    TEMPLATE = f.read()


def set_prop(content, name, value):
    pat = re.compile(rf"(<name>{name}</name>\s*<type>Constant</type>\s*<value>)([^<]*)(</value>)")
    new, n = pat.subn(rf"\g<1>{value:.6e}\g<3>", content)
    if n != 1:
        raise RuntimeError(f"{name}: {n} Treffer")
    return new


def write_prj(k, S, path=PRJ):
    c = set_prop(set_prop(TEMPLATE, "permeability", k), "storage", S)
    with open(path, "w", encoding="ISO-8859-1") as f:
        f.write(c)


def run_ogs():
    r = subprocess.run([OGS, "-l", "warn", PRJ_NAME], cwd=SIM_DIR, capture_output=True, text=True, timeout=600)
    if r.returncode != 0:
        raise RuntimeError(r.stdout[-2000:])


_well_id = None
def sim_pressure():
    global _well_id
    root = ET.parse(PVD).getroot()
    ents = sorted((float(d.attrib["timestep"]), os.path.join(RES, d.attrib["file"])) for d in root.iter("DataSet"))
    t, p = [], []
    for ti, fn in ents:
        m = pv.read(fn)
        if _well_id is None:
            _well_id = int(np.argmin(np.linalg.norm(m.points - WELL, axis=1)))
        t.append(ti); p.append(m.point_data["pressure"][_well_id] - P_REF)
    return np.array(t), np.array(p)


def load_meas():
    return bh40.load_meas()   # 40.6: Baseline Median t < 1780 s


MT, MP = load_meas()
LOG = []


def evaluate(k, S):
    write_prj(k, S)
    try:
        run_ogs()
    except RuntimeError:
        LOG.append((k, S, np.nan))
        print(f"k={k:.4e}  S={S:.4e}  OGS nicht konvergiert", flush=True)
        return 1e12, None, None
    t, p = sim_pressure()
    err = float(np.sqrt(np.mean((p - np.interp(t, MT, MP)) ** 2)))
    LOG.append((k, S, err))
    print(f"k={k:.4e}  S={S:.4e}  RMSE={err:.4e} Pa", flush=True)
    return err, t, p


def obj(x):
    return evaluate(10 ** x[0], 10 ** x[1])[0]


def main():
    # 1) Vorsuche auf log-Raster
    grid = [(a, b) for a in np.arange(-16, -9.9, 1.0) for b in np.arange(-12, -3.9, 1.0)]
    for a, b in grid:
        obj((a, b))
    best = min((e for e in LOG if np.isfinite(e[2])), key=lambda e: e[2])
    x0 = [np.log10(best[0]), np.log10(best[1])]
    # 2) Nelder-Mead Verfeinerung
    res = minimize(obj, x0, method="Nelder-Mead",
                   options={"xatol": 1e-3, "fatol": 1.0, "maxfev": 200,
                            "initial_simplex": [x0, [x0[0] + 0.3, x0[1]], [x0[0], x0[1] + 0.3]]})
    k, S = 10 ** res.x[0], 10 ** res.x[1]
    err, t, p = evaluate(k, S)  # finaler Lauf -> results/ und prj passen zu Optimum
    pd.DataFrame(LOG, columns=["permeability", "storage", "rmse_Pa"]).to_csv(
        os.path.join(HERE, f"calibrate_k_S_{TAG}_log.csv"), index=False)
    pd.DataFrame({"time_s": t, "p_sim_Pa": p, "p_meas_Pa": np.interp(t, MT, MP)}).to_csv(
        os.path.join(HERE, f"calibrate_k_S_{TAG}_best.csv"), index=False)
    print(f"\nOPTIMUM: permeability={k:.4e} m2, storage={S:.4e} 1/Pa, RMSE={err:.4e} Pa")


if __name__ == "__main__":
    main()
