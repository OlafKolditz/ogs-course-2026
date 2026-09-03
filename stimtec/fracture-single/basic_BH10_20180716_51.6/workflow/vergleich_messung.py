"""
Vergleicht die simulierten Druckaufbaukurven (druckaufbau.csv, erzeugt von
auswertung_druckaufbau.py) mit der digitalisierten gemessenen Kurve
(bh10_gemessen_digitalisiert.csv, erzeugt von digitize_bh10.py, gelbe Kurve
aus BH10_20180716_51.6.png).

Annahme: Die y-Achse des gemessenen Diagramms ist in MPa (Grossenordnung
passt zu den simulierten Druecken in Pa/1e6 nahe der Injektionsstelle).
Falls die Messung tatsaechlich in bar vorliegt, WERT_SKALIERUNG anpassen.
"""

import os
import pandas as pd
import matplotlib.pyplot as plt

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
SIM_CSV = os.path.join(SCRIPT_DIR, "druckaufbau.csv")
MEAS_CSV = os.path.join(SCRIPT_DIR, "bh10_gemessen_digitalisiert.csv")
PLOT_OUT = os.path.join(SCRIPT_DIR, "vergleich_messung.png")
CSV_OUT = os.path.join(SCRIPT_DIR, "vergleich_messung.csv")

# Angenommene Einheit der gemessenen Kurve: MPa -> Pa
WERT_SKALIERUNG = 1.0e6
MEAS_LABEL = "gemessen BH10 20180716 (digitalisiert, angenommen MPa)"


def main():
    sim = pd.read_csv(SIM_CSV)
    meas = pd.read_csv(MEAS_CSV)

    # Baseline vor Stimulationsbeginn (Sensor-Nullpunktversatz) abziehen,
    # damit "Druckaufbau" (Delta ggue. Ausgangszustand) mit der Simulation
    # vergleichbar ist, die bei 0 startet.
    baseline_mask = meas["time_s"] < 2650
    baseline = meas.loc[baseline_mask, "value_axis_units"].median()
    print(f"Basislinie der Messung (t<2650s): {baseline:.3f} (Achseneinheiten) -> wird abgezogen")
    meas["value_relative"] = meas["value_axis_units"] - baseline
    meas["pressure_Pa"] = meas["value_relative"] * WERT_SKALIERUNG

    fig, ax = plt.subplots(figsize=(10, 6.5))

    pressure_cols = [c for c in sim.columns if c == "pressure_(0,0)_Pa"]
    for col in pressure_cols:
        label = col.replace("pressure_", "").replace("_Pa", "")
        ax.plot(sim["time_s"], sim[col], label=f"Simulation x,y = {label}")

    ax.plot(
        meas["time_s"],
        meas["pressure_Pa"],
        color="goldenrod",
        linewidth=2,
        label=MEAS_LABEL,
    )

    ax.set_xlabel("Zeit [s]")
    ax.set_ylabel("Druck [Pa]")
    ax.set_title("Vergleich: simulierter vs. gemessener Druckaufbau")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(PLOT_OUT, dpi=150)
    print(f"Plot geschrieben: {PLOT_OUT}")

    # Gemessene Kurve auf die Simulationszeitpunkte interpolieren, um eine
    # gemeinsame Vergleichstabelle zu erzeugen.
    import numpy as np

    merged = sim.copy()
    merged["pressure_gemessen_BH10_Pa"] = np.interp(
        sim["time_s"], meas["time_s"], meas["pressure_Pa"]
    )
    merged.to_csv(CSV_OUT, index=False)
    print(f"CSV geschrieben: {CSV_OUT}")


if __name__ == "__main__":
    main()
