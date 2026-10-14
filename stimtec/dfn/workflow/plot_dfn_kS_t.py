import numpy as np, pandas as pd, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
o = np.load('eval.npy', allow_pickle=True).item(); a = o['o_kst']; c = o['o_test_a1']
m = pd.read_csv('bh10_40_gemessen_digitalisiert.csv'); base = m.loc[m.time_s < 1780, 'value_axis_units'].median()
tm, pm = m.time_s.values, m.value_axis_units.values - base
st = pd.read_csv('injektion_stufen_40.csv'); tb = np.r_[1750, st.t_end_s.values]
sim = np.interp(tm, a[:, 0], a[:, 1:].mean(1)) / 1e6
print('Stufe [l/min] | t | mittl. Fehler Sim-Mess [MPa]')
for i in range(len(st)):
    s = (tm >= tb[i]) & (tm < tb[i+1]); print(f'{st.q_lpm[i]:6.3f} {tb[i]:7.1f}-{tb[i+1]:7.1f}  {np.mean(sim[s]-pm[s]):+.3f}')
fig, (ax, axr) = plt.subplots(2, 1, figsize=(9, 6.4), sharex=True, gridspec_kw=dict(height_ratios=[3, 1]))
ax.plot(tm, pm, color='0.45', lw=1, label='Messung 40.6 (Baseline abgezogen)')
ax.plot(c[:, 0], c[:, 1:].mean(1)/1e6, color='0.7', lw=1.4, ls='--', label='DFN k, S konstant (Mittel V4 2D): RMSE 1.24 MPa')
ax.plot(a[:, 0], a[:, 1]/1e6, color='#1f6feb', lw=1.8, label='DFN V4 k(t), S(t): Kluft 1 (Knoten 270)')
ax.plot(a[:, 0], a[:, 2]/1e6, color='#d4762c', lw=1.8, label='DFN V4 k(t), S(t): Kluft 2 (Knoten 271)')
ax.set_ylabel('Druckaufbau [MPa]'); ax.grid(alpha=.3); ax.set_xlim(1750, 2670)
ax2 = ax.twinx(); q = st.q_lpm.values
ax2.step(tb, np.r_[q, q[-1]], where='post', color='#2a9d5c', lw=1, alpha=.7, label='Injektionsrate'); ax2.set_ylabel('Rate [l/min]', color='#2a9d5c'); ax2.set_ylim(0, 25)
h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels(); ax.legend(h1+h2, l1+l2, loc='upper left', fontsize=8)
ax.set_title('DFN-Modell V4, Pumpversuch BH10 40.6: RMSE 0.099 MPa (Mittel beider Knoten)')
ax.text(0.99, 0.03, 'k: 2.98e-16 … 8.25e-15 m²\nS: 3.14e-10 … 1.77e-09 1/Pa\nApertur 1 m, je 50 % Rate', transform=ax.transAxes, ha='right', va='bottom', fontsize=8.5)
axr.plot(tm, sim - pm, color='#7a4fd1', lw=1); axr.axhline(0, color='0.3', lw=.8); axr.set_ylabel('Sim − Mess [MPa]'); axr.set_xlabel('Zeit [s]'); axr.grid(alpha=.3)
fig.tight_layout(); fig.savefig('dfn_40_kS_t_kalibrierung.png', dpi=130)
