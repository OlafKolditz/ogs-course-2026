# Pumpversuch BH10_20180718_40.6 – Kalibrierung permeability / storage (30.09.2026)

Gleiche Methodik wie BH10_20180716_51.6 (Varianten 1–5). Modell: LiquidFlow 2D, 100x100 m, Punktquelle "well" (0,0), OGS 6.5.x.

## 1. Digitalisierung (workflow/)
- Quelle: `BH10_20180718_40.6.pdf` (Matplotlib-Vektorgrafik) -> Kurven direkt aus den PDF-Pfaden gelesen (`digitize_bh10_40.py`, pdfplumber), Achsen über Tick-Marken kalibriert (keine Pixel-Digitalisierung nötig)
- `bh10_40_gemessen_digitalisiert.csv`: Druck "P int d" [MPa], 2852 Punkte, 1750–2670 s (gleiches Format wie 51.6)
- `bh10_40_rate_digitalisiert.csv`: Q downhole und Q downhole_smoothed [l/min]
- `bh10_40_stufen.csv`: rote Stufenmarkierungen der Grafik
- `injektion_stufen_40.py` -> `injektion_stufen_40.csv`: Stufengrenzen aus Sprungflanken der geglätteten Rate, Stufenwert = Mittel der Rohrate, q[m³/s] = q[l/min]/60000

| bis t [s] | 1785 | 1900.5 | 2071.5 | 2151.8 | 2322.8 | 2422.3 | 2510.3 | 2670 |
|---|---|---|---|---|---|---|---|---|
| q [l/min] | 0.003 | 0.548 | 1.118 | 2.066 | 4.061 | 6.086 | 7.994 | 0.011 |

- Baseline: Median t < 1780 s = 0.087 MPa; Plot `plot_messung.py` -> `bh10_40_gemessen.png`

## 2. Modelle (simulation/, erzeugt mit `workflow/make_prj_40.py` aus den 51.6-prj)
- Zeitraum 1750–2670 s, dt = 5 s (51.6: 10 s), pressure_source_term aus Tabelle oben
- Gemeinsame Definitionen (Stufen, Messung, Baseline): `workflow/bh40.py`
- Zielfunktion: RMSE des Drucks am Brunnen über 1750–2670 s (inkl. Druckabbau nach Shut-in)

## 3. Ergebnisse

| Variante | prj | Parameter | RMSE |
|---|---|---|---|
| V1 k, S konstant, geschlossen | stimtec.prj | k = 9.15e-15 m², S = 1.08e-9 1/Pa | 1.163 MPa |
| V2 k, S konstant, RB rechts 1 MPa | stimtec_bc_right.prj | k = 9.15e-15, S = 1.08e-9 | 1.163 MPa |
| V3 k(t), S(t), RB rechts | stimtec_bc_right_kS_t.prj | k 6.16e-16…1.73e-14, S 2.01e-10…7.95e-10 | **0.0975 MPa** |
| V4 k(t), S(t), RB Punkt (50,0) | stimtec_bc_rightmid_kS_t.prj | k 6.16e-16…1.73e-14, S 2.01e-10…7.98e-10 | **0.0975 MPa** |
| V5a Python "rate" (= V4) | stimtec_bc_rightmid_py.prj | wie V4 | 0.0975 MPa (identisch) |
| V5b Python "pressure_exp" | stimtec_bc_rightmid_py.prj | k0 = 3.68e-15, β_k = 1.06 1/MPa, S0 = 4.73e-14, β_S = 0.79 1/MPa | 0.211 MPa |

Aufrufe (im Ordner workflow):
- V1: `python calibrate_k_S.py`; V2: `PRJ_NAME=stimtec_bc_right.prj P_REF=1e6 python calibrate_k_S.py`
- V3: `python calibrate_kS_t.py` (Raster + Nelder-Mead, ~460 Läufe)
- V4: `PRJ_NAME=stimtec_bc_rightmid_kS_t.prj BC_GEOMETRY=right_mid X0="-15.2105,-13.7615,-9.6970,-9.0998" python calibrate_kS_t.py`
- V5: `py -3.13 ogs_py_kS.py` (rate, PAR_RATE = Optimum V4), `py -3.13 calibrate_py_pressure.py`
- Plots: `python plot_kalibrierung_40.py` -> `kalibrierung_40_varianten.png`; `python plot_expressions.py` -> `expressions_kS_t.png`

## 4. Befunde
- Wie bei 51.6 dominiert mit konstanten k/S der Strukturfehler (Rate ×14.6, Druck nur ×2.3). Mit k(t), S(t) gekoppelt an die Injektionsrate wird die Messung sehr gut getroffen (RMSE 0.10 MPa, Faktor 12 besser als V1; 51.6: 0.32 MPa).
- Randbedingung wirkungslos (V1 = V2, V3 = V4): die Diffusivität ist klein (k/(μS) ≈ 0.02 m²/s), die Druckfront erreicht bis zum Shut-in nur ~10 m (Δp(10 m) ≈ 0.04 MPa, Δp(20 m) < 1 Pa).
- Permeabilität steigt mit der Rate um Faktor ~28 (51.6: ~7); Storage ist etwa 15–50× größer als bei 51.6 und nur schwach ratenabhängig (Faktor 4) -> langsamer, transienter Druckanstieg innerhalb der Stufen und langsamer Druckabbau nach Shut-in.
- Restabweichungen V4 (mittl. Fehler je Stufe): Stufe 4 l/min (2152–2323 s) −0.12 MPa unterschätzt, Stufe 8 l/min (2422–2510 s) +0.10 MPa überschätzt (Messung stagniert bei ~3.75 MPa -> Hinweis auf zusätzliche Öffnung/Leckage bei hoher Rate), sonst |Fehler| < 0.06 MPa.
- pressure_exp (V5b): gute Gesamttendenz, aber Spitzen bei jedem Ratensprung durch explizite Kopplung k(p^n) – wie bei 51.6.
