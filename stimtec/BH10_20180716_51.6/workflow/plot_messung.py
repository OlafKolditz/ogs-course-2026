"""
Plot der digitalisierten BH10-Messkurve 51.6 (bh10_gemessen_digitalisiert.csv)
plus Injektionsrate (bh10_rate_digitalisiert.csv, rechte Achse).
Layout identisch zu BH10_20180716_40.6/workflow/plot_messung.py.

Spalten: time_s [s], value_axis_units [MPa, digitalisierte y-Achse]

Aufruf von der Konsole:
    python plot_messung.py                 # Rohwerte + Baseline-korrigiert -> bh10_gemessen.png
    python plot_messung.py --nur-roh       # nur Rohwerte
    python plot_messung.py --show          # zusaetzlich Fenster oeffnen

Aus Python:
    from plot_messung import plot_messung
    ax = plot_messung()                    # gibt die matplotlib-Achse zurueck
"""
import argparse
import os

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.join(HERE, "bh10_gemessen_digitalisiert.csv")
CSV_RATE = os.path.join(HERE, "bh10_rate_digitalisiert.csv")
T_BASELINE = 2650.0  # s, Ende der Ruhephase vor Stimulationsbeginn (51.6)


def lade_messung(csv=CSV, t_baseline=T_BASELINE):
    """Liest die Messung; ergaenzt die Spalte dp_MPa (Baseline = Median fuer t < t_baseline)."""
    d = pd.read_csv(csv)
    base = d.loc[d["time_s"] < t_baseline, "value_axis_units"].median()
    d["dp_MPa"] = d["value_axis_units"] - base
    d.attrs["baseline_MPa"] = base
    return d


def plot_messung(csv=CSV, baseline_korrigiert=True, roh=True, ax=None, datei=None):
    """Plottet die digitalisierte Messkurve.

    baseline_korrigiert: Druckaufbau dp = p - Baseline (wie in den Kalibrierungen)
    roh:                 digitalisierte Rohwerte (value_axis_units)
    ax:                  vorhandene matplotlib-Achse (optional)
    datei:               PNG-Dateiname zum Speichern (optional)
    """
    d = lade_messung(csv)
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 5), facecolor="#fcfcfb")
    ax.set_facecolor("#fcfcfb")
    if roh:
        ax.plot(d["time_s"], d["value_axis_units"], color="#52514e", lw=1.2,
                label="digitalisiert, Rohwerte")
    if baseline_korrigiert:
        ax.plot(d["time_s"], d["dp_MPa"], color="#eb6834", lw=2,
                label=f"Druckaufbau (Baseline {d.attrs['baseline_MPa']:.3f} MPa abgezogen)")
    if roh and baseline_korrigiert:
        ax.axhline(d.attrs["baseline_MPa"], color="#52514e", lw=0.8, ls="--")
    ax.axvline(T_BASELINE, color="#52514e", lw=0.8, ls=":")
    # Beschriftung senkrecht entlang der Linie, links davon in der freien Flaeche vor Pumpbeginn
    ax.annotate(f"Stimulationsbeginn (t = {T_BASELINE:.0f} s)", xy=(T_BASELINE, 0.5),
                xycoords=("data", "axes fraction"), xytext=(-4, 0), textcoords="offset points",
                rotation=90, ha="right", va="center", fontsize=9, color="#52514e")
    ax.set_xlabel("Zeit [s]")
    ax.set_ylabel("Druck [MPa]")
    if os.path.exists(CSV_RATE):
        r = pd.read_csv(CSV_RATE)
        ax2 = ax.twinx()
        ax2.plot(r["time_s"], r["q_lpm"], color="#f2a37f", lw=0.6, alpha=0.7)
        ax2.plot(r["time_s"], r["q_smoothed_lpm"], color="#1baf7a", lw=1.5, label="Injektionsrate (geglättet)")
        ax2.set_ylabel("Injektionsrate [l/min]")
        ax2.spines["top"].set_visible(False)
        ax2.legend(frameon=False, loc="upper left", bbox_to_anchor=(0.0, 0.84))
    ax.set_title("BH10, Versuch BH10_20180716_51.6 – digitalisierte Messkurve", loc="left")
    ax.grid(alpha=0.25)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(0.0, 1.0))
    if datei:
        ax.figure.tight_layout()
        ax.figure.savefig(datei, dpi=150)
    return ax


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--csv", default=CSV)
    ap.add_argument("--out", default=os.path.join(HERE, "bh10_gemessen.png"))
    ap.add_argument("--nur-roh", action="store_true", help="nur Rohwerte plotten")
    ap.add_argument("--show", action="store_true", help="Plotfenster oeffnen")
    a = ap.parse_args()
    if not a.show:
        matplotlib.use("Agg")
    plot_messung(a.csv, baseline_korrigiert=not a.nur_roh, datei=a.out)
    print(f"Plot gespeichert: {a.out}")
    if a.show:
        plt.show()
