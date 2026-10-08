"""
Digitalisiert die blaue Kurve (Injektionsrate) aus dem Diagramm-Screenshot
BH10_20180716_51.6.png (Ergaenzung zu digitize_bh10.py, das die gelbe Druckkurve liefert).

Vorgehen (identisch zu digitize_bh10.py):
 - Farberkennung der blauen Kurve (Excel-Standardfarbe RGB 91,155,213)
 - Achsenkalibrierung ueber die automatisch erkannten Gitterlinien (calibrate aus digitize_bh10.py);
   Druck und Rate teilen sich im Excel-Diagramm dieselbe y-Achse (Rate in l/min)
 - je Bildspalte Median der blauen Pixelzeilen = Rohrate (Mitte des Rauschbandes)
 - Spalten, in denen Blau von der gelben Kurve verdeckt ist, werden linear interpoliert
 - Glaettung: gleitender Median ueber 15 s (erhaelt die Stufenflanken)
Ausgabe (gleiches Format wie 40.6):
  bh10_rate_digitalisiert.csv   time_s, q_lpm, q_smoothed_lpm
"""
import os
import numpy as np
import pandas as pd
from PIL import Image

from digitize_bh10 import IMG_FILE, calibrate, extract_curve

HERE = os.path.dirname(os.path.abspath(__file__))
CSV_OUT = os.path.join(HERE, "bh10_rate_digitalisiert.csv")
BLUE = np.array([91, 155, 213])
COLOR_TOL = 60
SMOOTH_S = 15.0   # s, Fensterbreite gleitender Median


def main():
    arr = np.array(Image.open(IMG_FILE).convert("RGB"))
    (ax, bx), (ay, by) = calibrate(arr)
    print(f"x: t = {ax:.6f}*px + {bx:.3f}   y: v = {ay:.6f}*px + {by:.3f}")
    px_x, px_y = extract_curve(arr, BLUE, COLOR_TOL)
    # alle Bildspalten zwischen erster und letzter blauer Spalte (Luecken interpolieren)
    cols = np.arange(px_x.min(), px_x.max() + 1)
    py = np.interp(cols, px_x, px_y)
    t = ax * cols + bx
    q = np.clip(ay * py + by, 0.0, None)
    dt = np.median(np.diff(t))
    win = max(int(round(SMOOTH_S / dt)) | 1, 3)
    qs = pd.Series(q).rolling(win, center=True, min_periods=1).median().values
    pd.DataFrame({"time_s": t, "q_lpm": q, "q_smoothed_lpm": qs}).to_csv(CSV_OUT, index=False)
    print(f"{len(t)} Punkte, t {t.min():.1f}..{t.max():.1f} s, q {q.min():.3f}..{q.max():.3f} l/min, "
          f"{len(cols) - len(px_x)} Spalten interpoliert")
    print(f"CSV geschrieben: {CSV_OUT}")


if __name__ == "__main__":
    main()
