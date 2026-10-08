"""
Digitalisierung des Pumpversuchs BH10_20180718_40.6 (Grafik BH10_20180718_40.6.pdf).

Die PDF ist eine Matplotlib-Vektorgrafik. Statt Pixel-Farberkennung werden die
Linienpfade direkt aus der PDF gelesen (pdfplumber) -> exakte Stuetzpunkte.
 - Achsenkalibrierung ueber die Tick-Marken (x: 1800..2600 s, y: 0..8)
 - Kurven nach Linienfarbe (Matplotlib tab10):
     blau   P int d            Druck [MPa]
     orange Q downhole         Rate [l/min]
     gruen  Q downhole_smoothed Rate [l/min]
 - rote gestrichelte Linien = Stufengrenzen [s]
Ausgabe:
  bh10_40_gemessen_digitalisiert.csv   time_s, value_axis_units  (Druck, MPa; gleiches Format wie 51.6)
  bh10_40_rate_digitalisiert.csv       time_s, q_lpm, q_smoothed_lpm
  bh10_40_stufen.csv                   Stufengrenzen (rote Linien)
"""
import os
import numpy as np
import pandas as pd
import pdfplumber

HERE = os.path.dirname(os.path.abspath(__file__))
PDF = os.path.join(HERE, "..", "BH10_20180718_40.6.pdf")   # liegt im Versuchsordner

COLORS = {"p": (0.1215686275, 0.4666666667, 0.7058823529),
          "q": (1.0, 0.4980392157, 0.0549019608),
          "qs": (0.1725490196, 0.6274509804, 0.1725490196)}
RED = (0.8392156863, 0.1529411765, 0.1568627451)
X_TICKS = [1800, 2000, 2200, 2400, 2600]
Y_TICKS = [0, 1, 2, 3, 4, 5, 6, 7, 8]


def close(a, b):
    return isinstance(a, tuple) and len(a) == 3 and np.allclose(a, b, atol=1e-3)


def main():
    page = pdfplumber.open(PDF).pages[0]
    # Tick-Marken: kurze schwarze Linien ausserhalb der Achsenbox
    xt = sorted(l["x0"] for l in page.lines if abs(l["x1"] - l["x0"]) < 1e-6 and 3 < l["bottom"] - l["top"] < 4)
    yt = sorted((l["top"] for l in page.lines if abs(l["bottom"] - l["top"]) < 1e-6 and 3 < l["x1"] - l["x0"] < 4),
                reverse=True)
    ax, bx = np.polyfit(xt, X_TICKS, 1)
    ay, by = np.polyfit(yt, Y_TICKS, 1)
    print(f"x: t = {ax:.6f}*px + {bx:.3f}   y: v = {ay:.6f}*px + {by:.3f}")

    curves = {}
    for c in page.curves:
        for k, col in COLORS.items():
            if close(c["stroking_color"], col) and len(c["pts"]) > 100:   # Legendenlinien ignorieren
                pts = np.array(c["pts"])
                curves[k] = (ax * pts[:, 0] + bx, ay * pts[:, 1] + by)
    for k, (t, v) in curves.items():
        print(f"{k}: {len(t)} Punkte, t {t.min():.1f}..{t.max():.1f} s, Wert {v.min():.3f}..{v.max():.3f}")

    red = sorted(ax * l["x0"] + bx for l in page.lines if close(l["stroking_color"], RED))
    print("Stufengrenzen (rot):", [round(r, 1) for r in red])

    t, p = curves["p"]
    pd.DataFrame({"time_s": t, "value_axis_units": p}).to_csv(
        os.path.join(HERE, "bh10_40_gemessen_digitalisiert.csv"), index=False)
    tq, q = curves["q"]
    tqs, qs = curves["qs"]
    pd.DataFrame({"time_s": tq, "q_lpm": q, "q_smoothed_lpm": np.interp(tq, tqs, qs)}).to_csv(
        os.path.join(HERE, "bh10_40_rate_digitalisiert.csv"), index=False)
    pd.DataFrame({"t_red_s": red}).to_csv(os.path.join(HERE, "bh10_40_stufen.csv"), index=False)


if __name__ == "__main__":
    main()
