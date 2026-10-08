"""
Gemeinsame Definitionen fuer den Pumpversuch BH10_20180718_40.6
(Grafik BH10_20180716_40.6/workflow/BH10_20180718_40.6.pdf).

 - Messung: bh10_40_gemessen_digitalisiert.csv (Druck P int d, MPa), Baseline = Median t < 1780 s
 - Injektion: injektion_stufen_40.csv (aus injektion_stufen_40.py)
 - Simulationszeitraum 1750..2670 s, dt = 5 s
"""
import os
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
T_START, T_END, DT = 1750.0, 2670.0, 5.0
T_BASE = 1780.0
MEAS_FILE = os.path.join(HERE, "bh10_40_gemessen_digitalisiert.csv")

_st = pd.read_csv(os.path.join(HERE, "injektion_stufen_40.csv"))
T_BOUNDS = _st["t_end_s"].values[:-1].astype(float)   # Stufe i gilt fuer t < T_BOUNDS[i]
Q_STEPS = _st["q_m3s"].values.astype(float)            # letzte Stufe: nach Shut-in bis Simulationsende
Q_MIN, Q_MAX = Q_STEPS.min(), Q_STEPS.max()


def load_meas():
    """Zeit [s], Druckaufbau [Pa] (Baseline abgezogen)."""
    d = pd.read_csv(MEAS_FILE)
    base = d.loc[d["time_s"] < T_BASE, "value_axis_units"].median()
    return d["time_s"].values, (d["value_axis_units"].values - base) * 1e6


def expression(values):
    """?:-Expression im Stil von pressure_source_term fuer Stufenwerte 'values'."""
    lines = []
    for i, (tb, v) in enumerate(zip(T_BOUNDS, values[:-1])):
        pre = "(" if i == 0 else " "
        lines.append(f"\t\t\t\t{pre}t &lt; {tb:.1f} ? {v:.8e} :")
    lines.append(f"\t\t\t\t {values[-1]:.8e})")
    return "\n" + "\n".join(lines) + "\n\t\t\t"


def mapped_values(x_min, x_max):
    """Lineare Abbildung der Stufenraten auf [x_min, x_max] (wie Variante 3/4 bei 51.6)."""
    return x_min + (x_max - x_min) * (Q_STEPS - Q_MIN) / (Q_MAX - Q_MIN)


def rate_fraction(t):
    i = np.searchsorted(T_BOUNDS, t, side="right")
    return (Q_STEPS[i] - Q_MIN) / (Q_MAX - Q_MIN)


def injection_rate(t):
    return Q_STEPS[np.searchsorted(T_BOUNDS, t, side="right")]
