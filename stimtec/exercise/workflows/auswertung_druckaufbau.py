"""
Auswertung des Druckaufbaus ueber die Zeit an ausgewaehlten Punkten
einer OGS-Simulation (PVD/VTU-Zeitreihe).

Liest _out/stimtec.pvd, extrahiert an den Koordinaten
[0,0], [1,0], [2,0], [5,0], [10,0] den Punktdatensatz "pressure"
fuer jeden Zeitschritt und erzeugt:
  - druckaufbau.csv   (Zeit, Druck je Punkt)
  - druckaufbau.png   (Plot Druck vs. Zeit)
"""

import os
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd
import pyvista as pv
import matplotlib.pyplot as plt

# --- Konfiguration -----------------------------------------------------

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../results")
PVD_FILE = os.path.join(OUT_DIR, "stimtec.pvd")

# Auszuwertende Koordinaten [x, y] (z = 0)
POINTS = [
    (0.0, 0.0),
    (1.0, 0.0),
    (2.0, 0.0),
    (5.0, 0.0),
    (10.0, 0.0),
]

PRESSURE_FIELD = "pressure"

CSV_OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "druckaufbau.csv")
PLOT_OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "druckaufbau.png")


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


# --- Auswertung ------------------------------------------------------

def make_probe_points(target_points):
    """Erzeugt eine PolyData-Punktwolke aus den Zielkoordinaten (fuer FEM-Interpolation)."""
    coords = np.array([[x, y, 0.0] for (x, y) in target_points])
    return pv.PolyData(coords)


def main():
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

    # Ergebnisse als DataFrame / CSV
    df = pd.DataFrame({"time_s": times})
    for label in labels:
        df[f"pressure_{label}_Pa"] = pressures[label]

    df.to_csv(CSV_OUT, index=False)
    print(f"\nCSV geschrieben: {CSV_OUT}")

    # Plot
    fig, ax = plt.subplots(figsize=(9, 6))
    for label in labels:
        ax.plot(df["time_s"], df[f"pressure_{label}_Pa"], marker="o", markersize=3, label=f"x,y = {label}")

    ax.set_xlabel("Zeit [s]")
    ax.set_ylabel("Druck [Pa]")
    ax.set_title("Druckaufbau ueber die Zeit an ausgewaehlten Punkten")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(PLOT_OUT, dpi=150)
    print(f"Plot geschrieben: {PLOT_OUT}")


if __name__ == "__main__":
    main()
