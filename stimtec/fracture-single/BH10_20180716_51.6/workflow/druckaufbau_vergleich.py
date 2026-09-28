"""
Auswertung und Validierung des Druckaufbaus einer OGS-Simulation
(PVD/VTU-Zeitreihe) gegen eine digitalisierte Feldmessung.

Kombiniert die vormals getrennten Skripte auswertung_druckaufbau.py und
vergleich_messung.py:

  1. Liest results/stimtec.pvd, extrahiert an den Koordinaten
     [0,0], [1,0], [2,0], [5,0], [10,0] den Punktdatensatz "pressure"
     fuer jeden Zeitschritt (FEM-Interpolation) und schreibt:
       - druckaufbau.csv   (Zeit, Druck je Punkt)
       - druckaufbau.png   (Plot Druck vs. Zeit, alle Punkte)

  2. Vergleicht die simulierte Kurve an x,y=(0,0) mit der digitalisierten
     gemessenen Kurve (bh10_gemessen_digitalisiert.csv, erzeugt von
     digitize_bh10.py, gelbe Kurve aus BH10_20180716_51.6.png) und schreibt:
       - vergleich_messung.csv
       - vergleich_messung.png

Annahme: Die y-Achse des gemessenen Diagramms ist in MPa (Groessenordnung
passt zu den simulierten Druecken in Pa/1e6 nahe der Injektionsstelle).
Falls die Messung tatsaechlich in bar vorliegt, WERT_SKALIERUNG anpassen.
Die Messung hat zudem eine Sensor-Baseline vor Stimulationsbeginn
(t<2650s), die abgezogen wird, damit "Druckaufbau" (Delta ggue.
Ausgangszustand) mit der Simulation vergleichbar ist.
"""

import os
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd
import pyvista as pv
import matplotlib.pyplot as plt

# --- Konfiguration -----------------------------------------------------

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(SCRIPT_DIR, "../results")
PVD_FILE = os.path.join(OUT_DIR, "stimtec.pvd")

# Auszuwertende Koordinaten [x, y] (z = 0)
POINTS = [
    (0.0, 0.0),
    (1.0, 0.0),
    (2.0, 0.0),
    (5.0, 0.0),
    (10.0, 0.0),
]
VERGLEICHSPUNKT = "(0,0)"  # welcher der obigen Punkte mit der Messung verglichen wird

PRESSURE_FIELD = "pressure"

DRUCKAUFBAU_CSV = os.path.join(SCRIPT_DIR, "druckaufbau.csv")
DRUCKAUFBAU_PLOT = os.path.join(SCRIPT_DIR, "druckaufbau.png")

MEAS_CSV = os.path.join(SCRIPT_DIR, "bh10_gemessen_digitalisiert.csv")
VERGLEICH_CSV = os.path.join(SCRIPT_DIR, "vergleich_messung.csv")
VERGLEICH_PLOT = os.path.join(SCRIPT_DIR, "vergleich_messung.png")

WERT_SKALIERUNG = 1.0e6  # angenommene Einheit der Messung: MPa -> Pa
MEAS_LABEL = "gemessen BH10 20180716 (digitalisiert, angenommen MPa)"


# --- PVD einlesen --------------------------------------------------------

def read_pvd(pvd_path):
    """Liest eine PVD-Datei und gibt Liste von (zeit, vtu_pfad) zurueck."""
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


# --- Schritt 1: Druckaufbau an ausgewaehlten Punkten auswerten ----------

def make_probe_points(target_points):
    """Erzeugt eine PolyData-Punktwolke aus den Zielkoordinaten (fuer FEM-Interpolation)."""
    coords = np.array([[x, y, 0.0] for (x, y) in target_points])
    return pv.PolyData(coords)


def auswertung_druckaufbau():
    if not os.path.isfile(PVD_FILE):
        raise FileNotFoundError(f"PVD-Datei nicht gefunden: {PVD_FILE}")

    entries = read_pvd(PVD_FILE)
    if not entries:
        raise RuntimeError("Keine Zeitschritte in der PVD-Datei gefunden.")

    print(f"{len(entries)} Zeitschritte gefunden.")

    first_mesh = pv.read(entries[0][1])
    if PRESSURE_FIELD not in first_mesh.point_data:
        raise KeyError(
            f"Feld '{PRESSURE_FIELD}' nicht in Punktdaten gefunden. "
            f"Verfuegbare Felder: {list(first_mesh.point_data.keys())}"
        )

    labels = [f"({x:g},{y:g})" for (x, y) in POINTS]
    probe_points = make_probe_points(POINTS)

    # Testlauf am ersten Zeitschritt, um sicherzustellen, dass alle Punkte
    # innerhalb der Netzgebietsgrenzen liegen (sample() interpoliert per FEM-Formfunktionen).
    # Beachte: probe_points.sample(mesh) sampelt Daten von "mesh" auf die
    # Punkte von "probe_points" (nicht umgekehrt).
    test_sample = probe_points.sample(first_mesh)
    valid = test_sample["vtkValidPointMask"]
    for label, ok in zip(labels, valid):
        if not ok:
            print(f"Warnung: Zielpunkt {label} liegt ausserhalb des Netzes!")

    # Zeitreihe auswerten (FEM-Interpolation an den exakten Zielkoordinaten)
    times = []
    pressures = {label: [] for label in labels}

    for t, vtu_path in entries:
        mesh = pv.read(vtu_path)
        sampled = probe_points.sample(mesh)
        p = sampled[PRESSURE_FIELD]
        times.append(t)
        for label, val in zip(labels, p):
            pressures[label].append(val)

    df = pd.DataFrame({"time_s": times})
    for label in labels:
        df[f"pressure_{label}_Pa"] = pressures[label]

    df.to_csv(DRUCKAUFBAU_CSV, index=False)
    print(f"\nCSV geschrieben: {DRUCKAUFBAU_CSV}")

    fig, ax = plt.subplots(figsize=(9, 6))
    for label in labels:
        ax.plot(df["time_s"], df[f"pressure_{label}_Pa"], marker="o", markersize=3, label=f"x,y = {label}")

    ax.set_xlabel("Zeit [s]")
    ax.set_ylabel("Druck [Pa]")
    ax.set_title("Druckaufbau ueber die Zeit an ausgewaehlten Punkten")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(DRUCKAUFBAU_PLOT, dpi=150)
    plt.close(fig)
    print(f"Plot geschrieben: {DRUCKAUFBAU_PLOT}")

    return df


# --- Schritt 2: Vergleich mit der Messung -------------------------------

def vergleich_messung(sim_df):
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

    sim_col = f"pressure_{VERGLEICHSPUNKT}_Pa"
    ax.plot(sim_df["time_s"], sim_df[sim_col], label=f"Simulation x,y = {VERGLEICHSPUNKT}")

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
    fig.savefig(VERGLEICH_PLOT, dpi=150)
    plt.close(fig)
    print(f"Plot geschrieben: {VERGLEICH_PLOT}")

    # Gemessene Kurve auf die Simulationszeitpunkte interpolieren, um eine
    # gemeinsame Vergleichstabelle zu erzeugen.
    merged = sim_df.copy()
    merged["pressure_gemessen_BH10_Pa"] = np.interp(
        sim_df["time_s"], meas["time_s"], meas["pressure_Pa"]
    )
    merged.to_csv(VERGLEICH_CSV, index=False)
    print(f"CSV geschrieben: {VERGLEICH_CSV}")


def main():
    sim_df = auswertung_druckaufbau()
    vergleich_messung(sim_df)


if __name__ == "__main__":
    main()
