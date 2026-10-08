"""
Vergleich Messung (40.6) vs. kalibrierte Simulationen aller Varianten.
Liest die *_best.csv der Kalibrierskripte (Spalten time_s, p_sim_Pa, p_meas_Pa).

Aufruf: python plot_kalibrierung_40.py  -> kalibrierung_40_varianten.png
        python plot_kalibrierung_40.py calibrate_k_S_stimtec_best.csv  -> Einzelplot
"""
import os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd
import bh40

HERE = bh40.HERE
VARIANTEN = [  # (Datei, Beschriftung, Farbe)
    ("calibrate_k_S_stimtec_best.csv", "V1 k, S konstant, geschlossen", "#9a9892"),
    ("calibrate_k_S_stimtec_bc_right_best.csv", "V2 k, S konstant, RB rechts", "#c9a227"),
    ("calibrate_stimtec_bc_right_kS_t_best.csv", "V3 k(t), S(t), RB rechts", "#2a78d6"),
    ("calibrate_stimtec_bc_rightmid_kS_t_best.csv", "V4 k(t), S(t), RB Punkt (50,0)", "#eb6834"),
    ("calibrate_py_pressure_best.csv", "V5 Python k(p), S(p)", "#1baf7a"),
]


def main(dateien=None, out="kalibrierung_40_varianten.png"):
    MT, MP = bh40.load_meas()
    fig, ax = plt.subplots(figsize=(10, 5.5), facecolor="#fcfcfb")
    ax.set_facecolor("#fcfcfb")
    ax.plot(MT, MP / 1e6, color="#52514e", lw=1.0, label="Messung 40.6 (Baseline abgezogen)")
    for fn, lab, col in VARIANTEN:
        if dateien and fn not in dateien:
            continue
        f = os.path.join(HERE, fn)
        if not os.path.exists(f):
            continue
        d = pd.read_csv(f)
        rmse = np.sqrt(np.mean((d.p_sim_Pa - d.p_meas_Pa) ** 2)) / 1e6
        ax.plot(d.time_s, d.p_sim_Pa / 1e6, color=col, lw=2, label=f"{lab}: RMSE {rmse:.3f} MPa")
    ax.set_xlabel("Zeit [s]"); ax.set_ylabel("Druckaufbau am Brunnen [MPa]")
    ax.set_title("BH10_20180718_40.6 – Kalibrierung permeability / storage", loc="left")
    ax.grid(alpha=0.25)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.legend(frameon=False, loc="upper left", fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, out), dpi=150)
    print("Plot:", out)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        main(sys.argv[1:], out=os.path.splitext(sys.argv[1])[0].replace("calibrate_", "kalibrierung_40_").replace("_best", "") + ".png")
    else:
        main()
