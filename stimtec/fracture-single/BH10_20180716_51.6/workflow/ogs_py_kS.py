"""
STIMTEC – OGS-Simulation ueber die Python-Schnittstelle (ogs.OGSSimulator)
mit in Python berechneter Permeabilitaet und Storage.

Projektdatei: simulation/stimtec_bc_rightmid_py.prj
  - Mesh simulation/stimtec-square-py.vtu mit Zelldaten permeability_py, storage_py
  - permeability/storage im Medium: <type>Parameter</type> -> MeshElement-Parameter
    auf diesen Zelldaten
  - Dirichlet p = 1 MPa im Punkt right_mid (50,0,0), Anfangsdruck 1 MPa

Ablauf je Zeitschritt:
  1. Python berechnet k und S fuer den naechsten Zeitpunkt t_{n+1}
     (optional abhaengig vom aktuellen Druckfeld p^n)
  2. Werte werden direkt in den OGS-Speicher (numpy-View) geschrieben
  3. OGS rechnet den Zeitschritt (execute_time_step)
  4. Druck am Brunnen wird direkt aus dem OGS-Speicher gelesen

Die Funktionen injection_rate / permeability / storage ersetzen die
?:-Expressions der prj-Dateien; eigene Ansaetze (z.B. k(p)) koennen durch
Ersetzen dieser Funktionen eingebaut werden (siehe MODELS unten).
"""
import os
import numpy as np
import ogs.OGSMesh  # noqa: F401  (registriert den OGSMesh-Typ fuer sim.mesh())
from ogs import OGSSimulator

HERE = os.path.dirname(os.path.abspath(__file__))
SIM_DIR = os.path.normpath(os.path.join(HERE, "..", "simulation"))
RES_DIR = os.path.normpath(os.path.join(HERE, "..", "results"))
PRJ = os.path.join(SIM_DIR, "stimtec_bc_rightmid_py.prj")
MESH_NAME = "stimtec-square-py"
P_REF = 1.0e6                     # Anfangs-/Randdruck [Pa]
WELL = np.array([0.0, 0.0, 0.0])

# ---------------------------------------------------------------------------
# Injektionsrate: identisch zu <parameter> pressure_source_term in der prj
# ---------------------------------------------------------------------------
T_BOUNDS = np.array([2663.923, 2689.323, 2706.473, 2800.123, 2942.223, 3081.073,
                     3195.523, 3308.323, 3466.023, 3575.273, 3665.573, 3740.823,
                     3819.173, 3908.773, 3991.373, 4149.97])
Q_STEPS = np.array([1.03227e-06, 4.20422763e-06, 1.04732717e-06, 2.70980080e-06,
                    2.20319288e-05, 2.73892046e-05, 3.44347213e-05, 5.11701849e-05,
                    6.94618453e-05, 8.71195366e-05, 1.00001857e-04, 1.11927032e-04,
                    1.31219025e-04, 1.46891195e-04, 1.61563301e-04, 9.76304533e-07])
Q_MIN, Q_MAX = Q_STEPS.min(), Q_STEPS.max()


def injection_rate(t):
    """q(t) wie pressure_source_term: Stufe i gilt fuer t < T_BOUNDS[i], danach 0."""
    i = np.searchsorted(T_BOUNDS, t, side="right")
    return Q_STEPS[i] if i < len(Q_STEPS) else 0.0


def rate_fraction(t):
    """Normierte Rate (q - q_min)/(q_max - q_min) in [0,1]; nach Injektionsende 0."""
    i = np.searchsorted(T_BOUNDS, t, side="right")
    if i >= len(Q_STEPS):
        return 0.0
    return (Q_STEPS[i] - Q_MIN) / (Q_MAX - Q_MIN)


# ---------------------------------------------------------------------------
# Modelle fuer Permeabilitaet und Storage
#   Signatur: f(t, dp_cells, par) -> Wert oder Array je Zelle
#   dp_cells = Ueberdruck p - P_REF je Zelle [Pa] (aus Zeitschritt n)
# ---------------------------------------------------------------------------
def k_rate(t, dp, par):
    """Wie prj-Expression permeability_t: linear in normierter Injektionsrate."""
    return par["k_min"] + (par["k_max"] - par["k_min"]) * rate_fraction(t)


def S_rate(t, dp, par):
    """Wie prj-Expression storage_t."""
    return par["S_min"] + (par["S_max"] - par["S_min"]) * rate_fraction(t)


def k_pressure_exp(t, dp, par):
    """Beispiel echter Druckabhaengigkeit (Kluftoeffnung): k = k0 * exp(beta * dp),
    begrenzt auf k_max. dp je Zelle -> raeumlich variable Permeabilitaet."""
    return np.minimum(par["k0"] * np.exp(par["beta_k"] * np.maximum(dp, 0.0)), par["k_cap"])


def S_pressure_exp(t, dp, par):
    return np.minimum(par["S0"] * np.exp(par["beta_S"] * np.maximum(dp, 0.0)), par["S_cap"])


MODELS = {
    "rate": (k_rate, S_rate),                    # identisch zu stimtec_bc_rightmid_kS_t.prj
    "pressure_exp": (k_pressure_exp, S_pressure_exp),
}

# Kalibriertes Optimum aus stimtec_bc_rightmid_kS_t.prj (RMSE 0.322 MPa)
PAR_RATE = dict(k_min=3.3679e-15, k_max=2.4110e-14, S_min=2.4790e-14, S_max=1.4138e-11)


# ---------------------------------------------------------------------------
# Simulation
# ---------------------------------------------------------------------------
def run(model="rate", par=None, write_vtu=True, verbose=False):
    """Fuehrt die OGS-Simulation Zeitschritt fuer Zeitschritt aus.
    Rueckgabe: (t, dp_well, k_well, S_well) als numpy-Arrays."""
    par = PAR_RATE if par is None else par
    k_fun, S_fun = MODELS[model] if isinstance(model, str) else model

    os.makedirs(RES_DIR, exist_ok=True)
    args = ["ogs", PRJ, "-o", RES_DIR]
    if not verbose:
        args += ["-l", "error"]
    sim = OGSSimulator.OGSSimulation(args)
    try:
        mesh = sim.mesh(MESH_NAME)
        k_arr = mesh.data_array("permeability_py", "double")   # direkte Views in OGS-Speicher
        S_arr = mesh.data_array("storage_py", "double")
        pts = np.asarray(mesh.points()).reshape(-1, 3)
        well = int(np.argmin(np.linalg.norm(pts - WELL, axis=1)))
        cell_mean = _CellAverager(mesh)

        def update(t_next):
            p = np.asarray(mesh.data_array("pressure", "double"))
            dp_cells = cell_mean(p) - P_REF
            k_arr[:] = k_fun(t_next, dp_cells, par)
            S_arr[:] = S_fun(t_next, dp_cells, par)

        dt = 10.0  # wie FixedTimeStepping in der prj
        t = sim.current_time()
        T, P, K, S = [t], [mesh.data_array("pressure", "double")[well] - P_REF], [np.nan], [np.nan]
        while sim.current_time() < sim.end_time() - 1e-9:
            t_next = sim.current_time() + dt
            update(t_next)
            status = sim.execute_time_step()
            if status != 0:
                raise RuntimeError(f"OGS-Zeitschritt bei t={t_next} fehlgeschlagen")
            T.append(sim.current_time())
            P.append(mesh.data_array("pressure", "double")[well] - P_REF)
            K.append(float(np.mean(k_arr))); S.append(float(np.mean(S_arr)))
    finally:
        sim.close()
    return np.array(T), np.array(P), np.array(K), np.array(S)


_NODES_PER_VTK_TYPE = {1: 1, 3: 2, 5: 3, 9: 4, 10: 4, 12: 8, 13: 6, 14: 5}


class _CellAverager:
    """Mittelwert eines Knotenfeldes je Zelle. OGSMesh.cells() liefert
    (flache Konnektivitaet, VTK-Zelltypen); das Mesh enthaelt Dreiecke und
    Randlinien (gemischte Zelltypen)."""

    def __init__(self, mesh):
        conn, types = (np.asarray(a) for a in mesh.cells())
        n = np.array([_NODES_PER_VTK_TYPE[int(ti)] for ti in types])
        self.conn = conn
        self.offsets = np.concatenate([[0], np.cumsum(n)[:-1]])
        self.n = n

    def __call__(self, nodal):
        return np.add.reduceat(nodal[self.conn], self.offsets) / self.n


if __name__ == "__main__":
    import pandas as pd
    t, dp, k, S = run("rate", PAR_RATE)
    meas = pd.read_csv(os.path.join(HERE, "bh10_gemessen_digitalisiert.csv"))
    base = meas.loc[meas["time_s"] < 2650, "value_axis_units"].median()
    mp = np.interp(t, meas["time_s"], (meas["value_axis_units"] - base) * 1e6)
    print(f"RMSE = {np.sqrt(np.mean((dp - mp) ** 2)):.4e} Pa")
    pd.DataFrame({"time_s": t, "p_sim_Pa": dp, "p_meas_Pa": mp, "k_mean": k, "S_mean": S}).to_csv(
        os.path.join(HERE, "ogs_py_kS_rate.csv"), index=False)
