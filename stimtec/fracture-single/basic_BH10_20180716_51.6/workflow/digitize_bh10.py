"""
Digitalisiert die gelbe Kurve (gemessene Werte) aus dem Diagramm-Screenshot
BH10_20180716_51.6.png und schreibt sie als Zeitreihe (Zeit, Wert) in eine CSV.

Vorgehen:
 - Farberkennung der gelben Kurve (Excel-Standardfarbe RGB 255,192,0)
 - Kalibrierung der Pixel->Datenwert-Achsen anhand der bekannten,
   automatisch erkannten Gitterlinien (x: 2700..4500 s in 200-er Schritten,
   y: 0..10 in 2er-Schritten)
 - je Bildspalte wird der Median der gelben Pixelzeilen als Kurvenwert
   verwendet
"""

import os
import numpy as np
from PIL import Image

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
IMG_FILE = os.path.join(SCRIPT_DIR, "BH10_20180716_51.6.png")
CSV_OUT = os.path.join(SCRIPT_DIR, "bh10_gemessen_digitalisiert.csv")

YELLOW = np.array([255, 192, 0])
COLOR_TOL = 60  # Summe der abs. RGB-Differenzen

# Bekannte Gitterlinienwerte fuer die Kalibrierung (robuster als die
# aeusseren Plotraender, die je nach Achsenpadding leicht abweichen).
X_GRID_VALUES = [2700, 2900, 3100, 3300, 3500, 3700, 3900, 4100, 4300, 4500]
Y_GRID_VALUES = [10, 8, 6, 4, 2, 0]

GRID_GRAY = np.array([217, 217, 217])
GRID_TOL = 10


def find_gridlines(mask_axis_sum, threshold, count_expected):
    """Findet Gruppen zusammenhaengender Indizes mit hoher Trefferzahl."""
    idx = np.where(mask_axis_sum > threshold)[0]
    groups = []
    cur = [idx[0]]
    for i in idx[1:]:
        if i - cur[-1] <= 3:
            cur.append(i)
        else:
            groups.append(cur)
            cur = [i]
    groups.append(cur)
    centers = [float(np.mean(g)) for g in groups]
    return centers


def calibrate(arr):
    h, w, _ = arr.shape
    grid_mask = np.abs(arr.astype(int) - GRID_GRAY).sum(axis=2) < GRID_TOL

    col_counts = grid_mask.sum(axis=0)
    x_centers_all = find_gridlines(col_counts, h * 0.3, None)
    # die letzten N (N = len(X_GRID_VALUES)) Zentren vor dem rechten Rand
    # entsprechen den Werten in X_GRID_VALUES (Randlinien werden ignoriert,
    # da sie nicht exakt auf einem Gitterwert liegen)
    x_centers = x_centers_all[-(len(X_GRID_VALUES) + 1):-1] if len(x_centers_all) > len(X_GRID_VALUES) else x_centers_all[: len(X_GRID_VALUES)]
    if len(x_centers) != len(X_GRID_VALUES):
        # Fallback: nimm die mittleren, am gleichmaessigsten verteilten Linien
        x_centers = x_centers_all[1:1 + len(X_GRID_VALUES)]

    row_counts = grid_mask.sum(axis=1)
    y_centers_all = find_gridlines(row_counts, w * 0.3, None)
    y_centers = y_centers_all[1:1 + len(Y_GRID_VALUES)]

    ax, bx = np.polyfit(x_centers, X_GRID_VALUES, 1)
    ay, by = np.polyfit(y_centers, Y_GRID_VALUES, 1)

    return (ax, bx), (ay, by)


def extract_curve(arr, color, tol=COLOR_TOL):
    diff = np.abs(arr.astype(int) - color).sum(axis=2)
    mask = diff < tol
    cols = np.where(mask.any(axis=0))[0]
    px_x = []
    px_y = []
    for c in cols:
        rows = np.where(mask[:, c])[0]
        px_x.append(c)
        px_y.append(np.median(rows))
    return np.array(px_x, dtype=float), np.array(px_y, dtype=float)


def main():
    img = Image.open(IMG_FILE).convert("RGB")
    arr = np.array(img)

    (ax, bx), (ay, by) = calibrate(arr)
    print(f"x-Kalibrierung: value = {ax:.6f} * px + {bx:.3f}")
    print(f"y-Kalibrierung: value = {ay:.6f} * px + {by:.3f}")

    px_x, px_y = extract_curve(arr, YELLOW)
    time_s = ax * px_x + bx
    value = ay * px_y + by

    order = np.argsort(time_s)
    time_s = time_s[order]
    value = value[order]

    print(f"{len(time_s)} digitalisierte Punkte")
    print(f"Zeitbereich: {time_s.min():.1f} .. {time_s.max():.1f} s")
    print(f"Wertebereich: {value.min():.3f} .. {value.max():.3f}")

    import pandas as pd

    df = pd.DataFrame({"time_s": time_s, "value_axis_units": value})
    df.to_csv(CSV_OUT, index=False)
    print(f"CSV geschrieben: {CSV_OUT}")


if __name__ == "__main__":
    main()
