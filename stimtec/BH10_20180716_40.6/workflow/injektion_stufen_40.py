"""
Leitet aus der digitalisierten Rate (bh10_40_rate_digitalisiert.csv) die
Injektionsstufen fuer pressure_source_term ab.
 - Stufengrenzen: Mitte der Sprungflanken der geglaetteten Rate (|dq/dt| > 0.05 lpm/s)
 - Stufenwert: Mittel der Rohrate (Q downhole) im Fenster, ohne +-5 s um die Flanken
 - Umrechnung l/min -> m^3/s: q / 60000 (wie 51.6)
Ausgabe: injektion_stufen_40.csv (t_end_s, q_lpm, q_m3s); letzter Eintrag gilt bis t_end der Simulation
"""
import os
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
d = pd.read_csv(os.path.join(HERE, "bh10_40_rate_digitalisiert.csv"))
t, q, qs = d.time_s.values, d.q_lpm.values, d.q_smoothed_lpm.values
tg = np.arange(t.min(), t.max(), 0.5)
g = np.interp(tg, t, qs)
idx = np.where(np.abs(np.gradient(g, tg)) > 0.05)[0]
groups, cur = [], [idx[0]]
for i in idx[1:]:
    if i - cur[-1] <= 4:
        cur.append(i)
    else:
        groups.append(cur); cur = [i]
groups.append(cur)
bounds = [0.5 * (tg[c[0]] + tg[c[-1]]) for c in groups]
edges = [t.min()] + bounds + [t.max()]
rows = []
for a, b in zip(edges[:-1], edges[1:]):
    m = (t > a + 5) & (t < b - 5)
    ql = max(float(np.mean(q[m])), 0.0)
    rows.append((round(b, 1), ql, ql / 60000.0))
df = pd.DataFrame(rows, columns=["t_end_s", "q_lpm", "q_m3s"])
df.to_csv(os.path.join(HERE, "injektion_stufen_40.csv"), index=False)
print(df.to_string())
