this is the single fracture model approach (rectangle fracture)

claude instructions

26.09.2026: Projekt: stimtec
Opus5.5

(OK1)
Optimiere die <properties> <name>permeability</name> und <name>storage</name> so dass die berechneten Druckwerte <process_variable>pressure</process_variable> bestmöglich mit den gemessenen Werten (2. Spalte in bh10_gemessen_digitalisiert.csv (Verzeichnis workflow) übereinstimmen.


# Kalibrierung permeability / storage (stimtec.prj), 26.09.2026

- Modell: simulation/stimtec.prj (LiquidFlow 2D, 100x100 m, Punktquelle "well" (0,0), keine Randbedingungen = geschlossene Ränder), OGS 6.5.9
- Messung: workflow/bh10_gemessen_digitalisiert.csv, Spalte 2 in MPa, Baseline (Median t<2650 s) abgezogen
- Zielfunktion: RMSE des Drucks am Brunnen (0,0), t = 2500–4000 s
- Methode: Log-Raster k 1e-16…1e-10, S 1e-12…1e-4, dann Nelder-Mead in log10 (Skript workflow/calibrate_k_S.py)
- **Optimum: permeability = 1.409e-14 m², storage = 1.419e-11 1/Pa, RMSE = 1.69 MPa** (in stimtec.prj eingetragen)
- Befund: Strukturfehler dominiert. Messung erreicht je Raten-Stufe schnell ein Plateau (quasi-stationär, nichtlinear in der Rate); das Modell mit konstanten k/S und geschlossenen Rändern akkumuliert Druck linear -> unterschätzt früh, überschätzt spät (Ende ~9.3 statt 6.7 MPa).
- Nächste Schritte: Fernfeld-Dirichlet-RB (p=0 am Rand), druck-/zeitabhängige Permeabilität (Fracture-Öffnung) oder Kalibrierung nur ausgewählter Zeitfenster.
- Numerik: Picard abstol=1e-6 (absolut) scheitert bei sehr kleinem S (hohe Drücke) -> in diesen Fällen OGS-Abbruch; ggf. reltol verwenden.


(OK2)
Ja, eine gute Idee. Setzte bitte eine Druckrandbedingung an der rechten Seite, <polyline id="1" name="right">, mit einem konstanten Druckwert von 1 MPa. Erstelle bitte eine neue Projektdatei (prj), damit wir vorherige Ergebnisse behalten.

## Variante 2: Dirichlet-RB rechts (simulation/stimtec_bc_right.prj)

- p = 1 MPa auf polyline "right" (x = 50 m); Anfangsdruck ebenfalls 1 MPa (p_init), Vergleich mit Messung über p − 1 MPa; Output-Prefix ../results/stimtec_bc_right
- Aufruf: `PRJ_NAME=stimtec_bc_right.prj P_REF=1e6 python calibrate_k_S.py`
- **Optimum: permeability = 1.465e-14 m², storage = 1.004e-11 1/Pa, RMSE = 1.685 MPa** (praktisch keine Verbesserung ggü. 1.69 MPa)
- Grund: bei optimaler Diffusivität (~1.5 m²/s) erreicht die Druckfront den 50 m entfernten Rand kaum; quasi-stationärer Zweig (k≈3e-14, S≤1e-13) liefert nur RMSE ≈ 2.05 MPa, weil der gemessene Druck nicht proportional zur Rate ist (Rate ×7, Druck ×2.2).
- Folgerung: lineares Modell mit konstanten k/S reicht nicht; nächster Schritt druckabhängige Permeabilität k(p) (Fracture-Öffnung) oder zeitabhängige k(t)/S(t).

(OK3)
Ja, Permeabilität und Storage müssen abhängig mit Injektionsdruck, <name>pressure_source_term</name>, gemacht werden. Verwendende bitte für die Druckabhängigkeiten von Permeabilität und Storage die gleiche <expression> wie für den Injektionsdruck, damit die zeitlichen Abhängigkeiten identisch sind.

## Variante 3: k(t), S(t) wie Injektion (simulation/stimtec_bc_right_kS_t.prj) – beste Variante

- Basis: stimtec_bc_right.prj. permeability/storage als `<type>Parameter</type>` -> Function-Parameter `permeability_t`, `storage_t`
- Deren expression = identische ?:-Struktur/Zeitgrenzen wie `pressure_source_term`; Stufenwerte linear abgebildet: X(t) = X_min + (X_max − X_min)·(q(t) − q_min)/(q_max − q_min); nach Injektionsende X_min
- 4 Parameter kalibriert (Raster + Nelder-Mead in log10), Skript workflow/calibrate_kS_t.py, ~450 OGS-Läufe
- **Optimum: k_min = 3.30e-15, k_max = 2.61e-14 m²; S_min = 3.02e-14, S_max = 5.22e-12 1/Pa; RMSE = 0.334 MPa** (Faktor 5 besser)
- Restabweichungen: Anfangsphase < 2800 s (Sim. ~0.35 MPa vor Pumpbeginn wegen kleiner Grundrate bei niedrigem k; erster Peak 2663–2689 s unterschätzt); leichte Überschätzung 2950–3080 s
- Interpretation: Permeabilität steigt mit der Rate um Faktor ~8 (Kluftöffnung); Storage stark variabel, schwächer bestimmt

(OK4)
Das passt schon hervorragend. Kannst Du als weitere Variante die Druckrandbedingung auf der rechten Seite anpassen. Die Randbedingung soll nicht auf der gesamten rechten Kante sondern nur im Mittelpunkt der Kante (50,0,0) gesetzt werden. Hierfür kannst die gml Datei verändern. Bitte erstelle hierfür auch ein neues Projektfile prj, damit wir vorherige Ergebnisse behalten.

## Variante 4: Punkt-RB (50,0,0) + k(t), S(t) (simulation/stimtec_bc_rightmid_kS_t.prj) – aktuell beste Variante

- gml: neuer Punkt id=5 (50,0,0) name="right_mid" in stimtec-square.gml (bestehende Geometrien unverändert); Dirichlet p = 1 MPa nur dort
- Aufruf: `PRJ_NAME=stimtec_bc_rightmid_kS_t.prj BC_GEOMETRY=right_mid X0="..." python calibrate_kS_t.py` (Nelder-Mead ab Optimum Variante 3, ~150 Läufe)
- Mit Parametern von Variante 3: RMSE 0.717 MPa -> RB-Geometrie ist bei diesen Parametern relevant
- **Optimum: k_min = 3.37e-15, k_max = 2.41e-14 m²; S_min = 2.48e-14, S_max = 1.41e-11 1/Pa; RMSE = 0.322 MPa** (leicht besser als 0.334)
- k nahezu unverändert, S_max ~2.7× größer (Punkt-RB entlastet weniger -> mehr Speicherung nötig)
- Restabweichungen wie Variante 3 (Anfangsphase, erster Peak)

(OK5)
Kannst Du die expressions für druckabhängige Permeabilitäten und Storage auch in Python-Funktionen überführen und dann die Python-Schnittstelle von OGS für die numerische Simulation hierfür verwenden. Damit können die druckabhängige Permeabilitäten und Storage noch flexibler berechnet werden.

# OGS-Python-Schnittstelle für k und S (26.09.2026)

Variante 5 der STIMTEC-Kalibrierung (Fortsetzung von claude/kalibrierung_k_S_2026-09-26.md).

## Aufbau
- prj: simulation/stimtec_bc_rightmid_py.prj (Basis Variante 4: Dirichlet 1 MPa im Punkt right_mid (50,0,0), p_init 1 MPa)
- Mesh: simulation/stimtec-square-py.vtu = stimtec-square.vtu + Zelldaten `permeability_py`, `storage_py`
- Medium: permeability/storage `<type>Parameter</type>` -> `<type>MeshElement</type>`-Parameter auf diesen Zelldaten
- Treiber: workflow/ogs_py_kS.py mit `ogs.OGSSimulator.OGSSimulation`; je Zeitschritt: k, S für t_{n+1} in Python berechnen (optional aus Zelldruck p^n), per numpy-View direkt in OGS-Speicher schreiben, `execute_time_step()`, Brunnendruck direkt aus Speicher lesen
- Wichtig: `import ogs.OGSMesh` vor `sim.mesh()` (sonst "Unregistered type: OGSMesh"); Mesh enthält gemischte Zelltypen (1796 Dreiecke + 80 Randlinien), cells() = (flache Konnektivität, VTK-Typen)
- Ein Lauf ~1.5 s, mehrfach im selben Prozess lauffähig

## Modelle (Funktionen f(t, dp_cells, par))
- `rate`: identisch zur prj-Expression (lineare Abbildung der Injektionsrate). Verifikation: RMSE 3.2240e5 Pa wie Variante 4, max. Abweichung 10 Pa (Rundung)
- `pressure_exp`: k = k0·exp(β_k·Δp), S = S0·exp(β_S·Δp) je Zelle (räumlich variabel). Kalibriert (workflow/calibrate_py_pressure.py, ~150 Läufe): k0 = 3.50e-15 m², β_k = 0.49 1/MPa, S0 = 3.63e-14 1/Pa, β_S = 1.04 1/MPa, **RMSE 0.469 MPa** (schlechter als rate-Modell 0.322)
- Schwächen pressure_exp: explizite Kopplung (k aus p^n) -> Spitze bei Ratensprung t≈2810 s; k bleibt nach Ratenrückgang hoch (keine Hysterese/Schließung). Mögliche Verbesserungen: Picard-Iteration k(p^{n+1}) oder kleinere Zeitschritte, Hysterese, Kombination k(q, p).

(OK6)
Wie starte ich OGS mit den Python-Funktionen von der Konsole?


Du startest nicht ogs.exe, sondern das Python-Skript. Es startet OGS über die Python-Schnittstelle und berechnet vor jedem Zeitschritt k und S.

Voraussetzung: das OGS-Python-Paket in derselben Python-Installation, aus der dein ogs.exe kommt (Python 3.13):

py -3.13 -m pip install ogs numpy pandas scipy

Ich habe mit OGS 6.5.9 getestet. Deine Version zeigt py -3.13 -m pip show ogs.

1. Modell „rate“ (wie die prj-Expression) rechnen:

cd C:\User\94-CLD\stimtec\workflow
py -3.13 ogs_py_kS.py

Das Skript gibt den RMSE aus (0,322 MPa) und schreibt ogs_py_kS_rate.csv. Die VTU-Ergebnisse landen in results\stimtec_bc_rightmid_py_*.vtu und results\stimtec_bc_rightmid_py.pvd.

2. Druckabhängiges Modell mit eigenen Parametern rechnen:

py -3.13 -c "import ogs_py_kS as o; t,dp,k,S = o.run('pressure_exp', dict(k0=3.4968e-15, beta_k=0.4897e-6, k_cap=1e-12, S0=3.6278e-14, beta_S=1.0434e-6, S_cap=1e-9)); print(dp.max())"

Mit o.run(..., verbose=True) siehst du die vollständige OGS-Ausgabe.

3. Druckabhängiges Modell neu kalibrieren:

py -3.13 calibrate_py_pressure.py

Das dauert etwa 5 bis 10 Minuten.

Wichtig: Wenn du ogs stimtec_bc_rightmid_py.prj direkt aufrufst, läuft die Rechnung zwar durch, aber ohne die Python-Funktionen. k und S bleiben dann die festen Startwerte aus stimtec-square-py.vtu (1e-14 m² bzw. 1e-11 1/Pa).

Soll ich dem Skript noch Kommandozeilen-Optionen geben, zum Beispiel py -3.13 ogs_py_kS.py --model pressure_exp --k0 3.5e-15 ...? Dann brauchst du den -c-Aufruf nicht mehr.

(OK: Präsentationsmaterial)

Plotte bitte die digitalisierte Messkurve bh10_gemessen_digitalisiert.csv mit einer Python Funktion.
Platziere den Text "Simulationsbeginn" besser in der Grafik.

Plotte die die expressions für permeability_t, storage_t und pressure_source_term in in stimtec_bc_rightmid_kS_t.prj im gleichen Design wie die digitalisierte Messkurve.


28.09.2026
- Übertragen von CLD nach EDU
- Auswertung testen

Daten:
- workflow/plot_messung.py > bh10_gemessen.png
- workflow/plot_expression > expressions_kS_t.png

Rechnen: (1)
- simulation/stimtec.prj OK
- workflow/auswertung_druckaufbau.py > druckaufbau.png
- workflow/vergleich_messung.py > vergleich_messung.png

Rechnen: (2)
- simulation/stimtec_bc_right.prj
- workflow/auswertung_druckaufbau_bc_right.py > druckaufbau.png
- workflow/vergleich_messung.py > vergleich_messung.png ! shift 1 MPa
> Änderung von AB und BC von 1e6 auf 0 Pa

Rechnen: (4) (best fit)
- simulation/stimtec_bc_rightmid_kS_t.prj
- workflow/auswertung_druckaufbau_bc_rightmid_kS_t.py > druckaufbau.png
- workflow/vergleich_messung.py > vergleich_messung.png ! shift 1 MPa
> Änderung von AB und BC von 1e6 auf 0 Pa

Rechnen: (5) (python version for best fit)
- workflow: py -3.13 ogs_py_kS.py
- workflow/auswertung_druckaufbau_python_version.py > druckaufbau.png
- workflow/vergleich_messung.py > vergleich_messung.png ! shift 1 MPa
> Änderung von AB und BC von 1e6 auf 0 Pa
