"""
Plot der zeitabhaengigen <expression>s aus einer OGS-Projektdatei
(Default: simulation/stimtec_bc_rightmid_kS_t.prj):
    pressure_source_term, permeability_t, storage_t
im gleichen Design wie plot_messung.py.

Die Expressions sind verschachtelte ?:-Stufenfunktionen der Form
    (t < t1 ? v1 : t < t2 ? v2 : ... : v_else)
und werden hier direkt aus der prj gelesen und ausgewertet.

Aufruf von der Konsole:
    python plot_expressions.py                       # -> expressions_kS_t.png
    python plot_expressions.py --prj ../simulation/stimtec_bc_right_kS_t.prj --show
"""
import argparse
import os
import re
import xml.etree.ElementTree as ET

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PRJ = os.path.join(HERE, "..", "simulation", "stimtec_bc_rightmid_kS_t.prj")
T_BASELINE = 2650.0  # s, Stimulationsbeginn (wie plot_messung.py)

NUM = r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?"
PARAMETER = [  # (Name in prj, Achsenbeschriftung, Farbe)
    ("pressure_source_term", "Injektionsrate q [m³/s]", "#eb6834"),
    ("permeability_t", "Permeabilität k [m²]", "#2a78d6"),
    ("storage_t", "Storage S [1/Pa]", "#1baf7a"),
]


def lese_expression(prj, name):
    """Liest die <expression> des Function-Parameters `name` aus der prj."""
    root = ET.parse(prj).getroot()
    for p in root.iter("parameter"):
        if p.findtext("name") == name:
            return p.findtext("expression")
    raise KeyError(f"Parameter {name} nicht in {prj}")


def parse_stufen(expr):
    """Zerlegt '(t < t1 ? v1 : t < t2 ? v2 : ... : v_else)' in Grenzen, Werte, Else-Wert."""
    paare = re.findall(rf"t\s*<\s*({NUM})\s*\?\s*({NUM})", expr)
    rest = re.search(rf":\s*({NUM})\s*\)\s*$", expr.strip())
    grenzen = np.array([float(a) for a, _ in paare])
    werte = np.array([float(b) for _, b in paare])
    return grenzen, werte, float(rest.group(1))


def werte_aus(expr, t):
    """Wertet die Stufen-Expression fuer ein Zeit-Array t aus (wie exprtk in OGS)."""
    grenzen, werte, sonst = parse_stufen(expr)
    i = np.searchsorted(grenzen, t, side="right")
    return np.where(i < len(werte), werte[np.minimum(i, len(werte) - 1)], sonst)


def plot_expressions(prj=PRJ, t_bereich=(2500.0, 4280.0), datei=None):
    t = np.linspace(*t_bereich, 20001)
    fig, axes = plt.subplots(len(PARAMETER), 1, figsize=(10, 8.5), sharex=True, facecolor="#fcfcfb")
    for ax, (name, ylabel, farbe) in zip(axes, PARAMETER):
        y = werte_aus(lese_expression(prj, name), t)
        ax.set_facecolor("#fcfcfb")
        ax.plot(t, y, color=farbe, lw=2, label=name)
        ax.axvline(T_BASELINE, color="#52514e", lw=0.8, ls=":")
        ax.set_ylabel(ylabel)
        ax.grid(alpha=0.25)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        ax.legend(frameon=False, loc="upper left")
        ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    axes[0].annotate("Stimulationsbeginn", xy=(T_BASELINE, 0.12),
                     xycoords=("data", "axes fraction"), xytext=(-4, 0), textcoords="offset points",
                     rotation=90, ha="right", va="bottom", fontsize=9, color="#52514e")
    axes[0].set_title(f"Zeitabhängige Expressions – {os.path.basename(prj)}", loc="left", pad=16)
    axes[-1].set_xlabel("Zeit [s]")
    fig.tight_layout()
    if datei:
        fig.savefig(datei, dpi=150)
    return axes


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--prj", default=PRJ)
    ap.add_argument("--out", default=os.path.join(HERE, "expressions_kS_t.png"))
    ap.add_argument("--show", action="store_true", help="Plotfenster oeffnen")
    a = ap.parse_args()
    if not a.show:
        matplotlib.use("Agg")
    plot_expressions(a.prj, datei=a.out)
    print(f"Plot gespeichert: {a.out}")
    if a.show:
        plt.show()
