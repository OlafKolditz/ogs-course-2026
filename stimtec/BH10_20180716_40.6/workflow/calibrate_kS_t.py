"""
Kalibrierung zeit-/injektionsabhaengiger Permeabilitaet und Storage.

Ausgangsmodell: simulation/stimtec_bc_right.prj (p = 1 MPa am rechten Rand).
Neues Modell:   simulation/stimtec_bc_right_kS_t.prj

permeability und storage werden als <type>Parameter</type> an die Function-
Parameter "permeability_t" und "storage_t" gekoppelt. Deren <expression> ist
exakt die verschachtelte ?:-Expression von "pressure_source_term" (gleiche
Zeitgrenzen), nur die Werte je Stufe werden linear abgebildet:

    X(t) = X_min + (X_max - X_min) * (q(t) - q_min) / (q_max - q_min)

q(t) = Injektionsrate der Stufe, q_min/q_max = kleinster/groesster Wert.
Damit sind die zeitlichen Abhaengigkeiten identisch zur Injektion.
Kalibriert werden 4 Parameter: log10 k_min, log10 k_max, log10 S_min, log10 S_max.
"""
import bh40  # 40.6: Stufen, Messung
import os, re, subprocess, shutil, xml.etree.ElementTree as ET
import numpy as np, pandas as pd, pyvista as pv
from scipy.optimize import minimize

HERE = os.path.dirname(os.path.abspath(__file__))
SIM_DIR = os.path.join(HERE, "..", "simulation")
SRC_PRJ = os.path.join(SIM_DIR, "stimtec_bc_right.prj")
# Variante per Umgebungsvariablen:
#   PRJ_NAME    Name der erzeugten prj (Default stimtec_bc_right_kS_t.prj)
#   BC_GEOMETRY Geometrie der Dirichlet-RB (Default "right" = ganze Kante; "right_mid" = Punkt (50,0,0))
#   X0          Startwerte "log10kmin,log10kmax,log10Smin,log10Smax" -> Rastersuche entfaellt
PRJ_NAME = os.environ.get("PRJ_NAME", "stimtec_bc_right_kS_t.prj")
BC_GEOMETRY = os.environ.get("BC_GEOMETRY", "right")
X0 = os.environ.get("X0")
PRJ = os.path.join(SIM_DIR, PRJ_NAME)
TAG = os.path.splitext(PRJ_NAME)[0]
RES = os.path.join(HERE, "..", "results")
PVD = os.path.join(RES, TAG + ".pvd")
MEAS = os.path.join(HERE, "bh10_40_gemessen_digitalisiert.csv")
OGS = shutil.which("ogs") or "ogs"
P_REF = 1.0e6                      # Anfangs-/Randdruck
WELL = np.array([0.0, 0.0, 0.0])

with open(SRC_PRJ, encoding="ISO-8859-1") as f:
    SRC = f.read()

# --- Injektions-Expression extrahieren ------------------------------------
m = re.search(r"(<name>pressure_source_term</name>\s*<type>Function</type>\s*<expression>)(.*?)(</expression>)",
              SRC, re.DOTALL)
Q_EXPR = m.group(2)
# Werte stehen jeweils nach "?" bzw. als letzter Wert (": 0)")
VAL_RE = re.compile(r"(\?\s*)([0-9.eE+-]+)")
Q_VALS = [float(v) for _, v in VAL_RE.findall(Q_EXPR)]
Q_MIN, Q_MAX = min(Q_VALS), max(Q_VALS)


def mapped_expression(x_min, x_max):
    """Gleiche Expression wie pressure_source_term, Werte linear auf [x_min, x_max]."""
    return bh40.expression(bh40.mapped_values(x_min, x_max))  # 40.6: Stufen aus injektion_stufen_40.csv

    def rep(mt):
        q = float(mt.group(2))
        x = x_min + (x_max - x_min) * (q - Q_MIN) / (Q_MAX - Q_MIN)
        return f"{mt.group(1)}{x:.6e}"
    e = VAL_RE.sub(rep, Q_EXPR)
    # letzter Else-Wert (nach Injektionsende, t >= 4149.97 s): x_min statt 0
    e = re.sub(r":\s*0\s*\)\s*$", f": {x_min:.6e})\n\t\t\t", e.rstrip() )
    return e


def build_template():
    c = SRC.replace("<prefix>../results/stimtec_bc_right</prefix>", f"<prefix>../results/{TAG}</prefix>")
    if BC_GEOMETRY != "right":
        c, n = re.subn(r"<geometry>right</geometry>", f"<geometry>{BC_GEOMETRY}</geometry>", c)
        assert n == 1, "BC-Geometrie"
    for prop, par in (("permeability", "permeability_t"), ("storage", "storage_t")):
        pat = re.compile(rf"<name>{prop}</name>\s*<type>Constant</type>\s*<value>[^<]*</value>")
        c, n = pat.subn(f"<name>{prop}</name>\n                    <type>Parameter</type>\n"
                        f"                    <parameter_name>{par}</parameter_name>", c)
        assert n == 1, prop
    blocks = ""
    for par in ("permeability_t", "storage_t"):
        blocks += (f"        <parameter>\n            <!-- @{par}@ -->\n            <name>{par}</name>\n"
                   f"            <type>Function</type>\n            <expression>@{par}_EXPR@</expression>\n"
                   f"        </parameter>\n")
    c = c.replace("    <parameters>\n", "    <parameters>\n" + blocks, 1)
    return c


TEMPLATE = build_template()


def write_prj(kmin, kmax, smin, smax):
    c = TEMPLATE
    c = c.replace("@permeability_t_EXPR@", mapped_expression(kmin, kmax))
    c = c.replace("@storage_t_EXPR@", mapped_expression(smin, smax))
    c = c.replace("<!-- @permeability_t@ -->",
                  f"<!-- k(t): gleiche Zeitfenster wie pressure_source_term, k_min={kmin:.4e}, k_max={kmax:.4e} m2 -->")
    c = c.replace("<!-- @storage_t@ -->",
                  f"<!-- S(t): gleiche Zeitfenster wie pressure_source_term, S_min={smin:.4e}, S_max={smax:.4e} 1/Pa -->")
    with open(PRJ, "w", encoding="ISO-8859-1") as f:
        f.write(c)


def run_ogs():
    r = subprocess.run([OGS, "-l", "warn", PRJ_NAME], cwd=SIM_DIR, capture_output=True, text=True, timeout=600)
    if r.returncode != 0:
        raise RuntimeError(r.stdout[-2000:])


_well = None
def sim_pressure():
    global _well
    root = ET.parse(PVD).getroot()
    ents = sorted((float(d.attrib["timestep"]), os.path.join(RES, d.attrib["file"])) for d in root.iter("DataSet"))
    t, p = [], []
    for ti, fn in ents:
        mm = pv.read(fn)
        if _well is None:
            _well = int(np.argmin(np.linalg.norm(mm.points - WELL, axis=1)))
        t.append(ti); p.append(mm.point_data["pressure"][_well] - P_REF)
    return np.array(t), np.array(p)


def load_meas():
    return bh40.load_meas()


MT, MP = load_meas()
LOG = []


def evaluate(x):
    kmin, kmax, smin, smax = 10 ** np.asarray(x)
    write_prj(kmin, kmax, smin, smax)
    try:
        run_ogs()
    except RuntimeError:
        LOG.append((kmin, kmax, smin, smax, np.nan)); return 1e12, None, None
    t, p = sim_pressure()
    err = float(np.sqrt(np.mean((p - np.interp(t, MT, MP)) ** 2)))
    LOG.append((kmin, kmax, smin, smax, err))
    print(f"k=[{kmin:.3e},{kmax:.3e}] S=[{smin:.3e},{smax:.3e}] RMSE={err:.4e}", flush=True)
    return err, t, p


def obj(x):
    return evaluate(x)[0]


def main():
    if X0:
        best = (None, [float(v) for v in X0.split(",")])
    else:
        best = grid_search()
    x0 = np.array(best[1], float)
    simplex = [x0] + [x0 + 0.4 * np.eye(4)[i] for i in range(4)]
    return finish(x0, simplex)


def grid_search():
    # Vorsuche: k_min, k_max, S (S_min = S_max) auf Raster
    best = (np.inf, None)
    for a in np.arange(-16, -12.9, 0.5):
        for b in np.arange(a, -11.9, 0.5):
            for s in (-13, -12, -11, -10, -9):
                e = obj([a, b, s, s])
                if e < best[0]:
                    best = (e, [a, b, s, s])
    return best


def finish(x0, simplex):
    res = minimize(obj, x0, method="Nelder-Mead",
                   options={"xatol": 2e-3, "fatol": 10.0, "maxfev": 400, "initial_simplex": simplex})
    err, t, p = evaluate(res.x)
    kmin, kmax, smin, smax = 10 ** res.x
    pd.DataFrame(LOG, columns=["k_min", "k_max", "S_min", "S_max", "rmse_Pa"]).to_csv(
        os.path.join(HERE, f"calibrate_{TAG}_log.csv"), index=False)
    pd.DataFrame({"time_s": t, "p_sim_Pa": p, "p_meas_Pa": np.interp(t, MT, MP)}).to_csv(
        os.path.join(HERE, f"calibrate_{TAG}_best.csv"), index=False)
    print(f"\nOPTIMUM: k_min={kmin:.4e} k_max={kmax:.4e} m2, S_min={smin:.4e} S_max={smax:.4e} 1/Pa, RMSE={err:.4e} Pa")


if __name__ == "__main__":
    main()
