import numpy as np, pandas as pd, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
o = np.load('eval.npy', allow_pickle=True).item(); a = o['o_split']; b = o['o_kst']
m = pd.read_csv('bh10_40_gemessen_digitalisiert.csv'); base = m.loc[m.time_s < 1780, 'value_axis_units'].median()
tm, pm = m.time_s.values, m.value_axis_units.values - base
st = pd.read_csv('injektion_stufen_40.csv'); tb = np.r_[1750, st.t_end_s.values]
p1 = np.interp(tm, a[:, 0], a[:, 1]) / 1e6; p2 = np.interp(tm, a[:, 0], a[:, 2]) / 1e6
print('Stufe [l/min] | Fehler Kluft1 | Kluft2 | Mittel [MPa]')
for i in range(len(st)):
    s = (tm >= tb[i]) & (tm < tb[i+1])
    print(f'{st.q_lpm[i]:6.3f}  {np.mean(p1[s]-pm[s]):+.3f}  {np.mean(p2[s]-pm[s]):+.3f}  {np.mean((p1[s]+p2[s])/2-pm[s]):+.3f}')
fig, (ax, axr) = plt.subplots(2, 1, figsize=(9, 6.4), sharex=True, gridspec_kw=dict(height_ratios=[3, 1]))
ax.plot(tm, pm, color='0.45', lw=1, label='Messung 40.6 (Baseline abgezogen)')
ax.plot(b[:, 0], b[:, 1]/1e6, color='#1f6feb', lw=1, ls=':', label='vorher 50/50: Kluft 1')
ax.plot(b[:, 0], b[:, 2]/1e6, color='#d4762c', lw=1, ls=':', label='vorher 50/50: Kluft 2')
ax.plot(a[:, 0], a[:, 1]/1e6, color='#1f6feb', lw=1.8, label='Aufteilung kalibriert: Kluft 1 (46.5 %)')
ax.plot(a[:, 0], a[:, 2]/1e6, color='#d4762c', lw=1.8, label='Aufteilung kalibriert: Kluft 2 (53.5 %)')
ax.set_ylabel('Druckaufbau [MPa]'); ax.grid(alpha=.3); ax.set_xlim(1750, 2670)
ax2 = ax.twinx(); q = st.q_lpm.values
ax2.step(tb, np.r_[q, q[-1]], where='post', color='#2a9d5c', lw=1, alpha=.7, label='Injektionsrate'); ax2.set_ylabel('Rate [l/min]', color='#2a9d5c'); ax2.set_ylim(0, 25)
h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels(); ax.legend(h1+h2, l1+l2, loc='upper left', fontsize=8)
ax.set_title('DFN V4 + Ratenaufteilung, Pumpversuch BH10 40.6: RMSE 0.118 MPa (beide Knoten)')
ax.text(0.99, 0.03, 'k: 3.40e-16 … 9.06e-15 m²\nS: 2.50e-10 … 9.15e-10 1/Pa\nKluft 1 / 2: 46.5 % / 53.5 %', transform=ax.transAxes, ha='right', va='bottom', fontsize=8.5)
axr.plot(tm, p1 - pm, color='#1f6feb', lw=1, label='Kluft 1 − Mess'); axr.plot(tm, p2 - pm, color='#d4762c', lw=1, label='Kluft 2 − Mess')
axr.axhline(0, color='0.3', lw=.8); axr.set_ylabel('Sim − Mess [MPa]'); axr.set_xlabel('Zeit [s]'); axr.grid(alpha=.3); axr.legend(fontsize=8, loc='upper left', ncol=2)
fig.tight_layout(); fig.savefig('dfn_40_kS_split_kalibrierung.png', dpi=130)
