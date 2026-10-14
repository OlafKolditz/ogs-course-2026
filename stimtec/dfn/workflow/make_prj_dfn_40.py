"""stimtec_dfn_40.prj: DFN-Modell mit Injektionsstufen und mittleren k/S des Pumpversuchs 40.6.

Basis: stimtec_dfn_st_rim.prj (aus make_prj_dfn_st.py --rim): nodale Quellterme an den
Schnittpunkten BH10 x Kluft 1/2, Dirichlet p0 (hydrostatisch) auf den Kluftraendern.

Aufruf:  python make_prj_dfn_40.py [basis.prj] [ziel.prj] [--aperture=1e-6] [--split1=0.5]
 - Injektion: ../BH10_20180716_40.6/workflow/injektion_stufen_40.csv (q_m3s), Massenrate = 1000*q [kg/s]
 - k, S: zeitgewichtete Mittel von k(t), S(t) der Kalibrierung V4 (1750-2670 s)
 - Zeit 1750-2670 s, dt = 5 s (wie 40.6)
"""
import os, re, sys
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
args = [a for a in sys.argv[1:] if not a.startswith('--')]
opts = dict(a[2:].split('=') for a in sys.argv[1:] if a.startswith('--') and '=' in a)
src = args[0] if len(args) > 0 else os.path.join(HERE, 'stimtec_dfn_st_rim.prj')
dst = args[1] if len(args) > 1 else os.path.join(HERE, 'stimtec_dfn_40.prj')
APERTURE = float(opts.get('aperture', 1e-6))
SPLIT1 = float(opts.get('split1', 0.5))
RHO = 1000.0

STUFEN = os.path.join(HERE, '..', 'BH10_20180716_40.6', 'workflow', 'injektion_stufen_40.csv')
if not os.path.exists(STUFEN):
    STUFEN = os.path.join(HERE, 'injektion_stufen_40.csv')
st = pd.read_csv(STUFEN)
T_BOUNDS = st['t_end_s'].values[:-1]
Q = st['q_m3s'].values

# V4-Kalibrierung 40.6 (stimtec_bc_rightmid_kS_t.prj), Stufenwerte k(t), S(t)
K_T = np.array([6.15800298e-16, 1.75458723e-15, 2.94568919e-15, 4.92619269e-15,
                9.09311722e-15, 1.33221033e-14, 1.73077985e-14, 6.33437583e-16])
S_T = np.array([2.00999695e-10, 2.41751712e-10, 2.84375847e-10, 3.55249082e-10,
                5.04364406e-10, 6.55700629e-10, 7.98330578e-10, 2.01630853e-10])
T_START, T_END, DT = 1750.0, 2670.0, 5.0
w = np.diff(np.r_[T_START, T_BOUNDS, T_END])
K_MEAN = float((K_T * w).sum() / w.sum())
S_MEAN = float((S_T * w).sum() / w.sum())

def expr(vals):
    s = '('
    for tb, v in zip(T_BOUNDS, vals[:-1]):
        s += f't &lt; {tb:.1f} ? {v:.8e} : '
    return s + f'{vals[-1]:.8e})'

s = open(src, encoding='ISO-8859-1').read()

# Quellterme: Stufen 40.6, Massenrate [kg/s]
for name, f in [('q_frac1', SPLIT1), ('q_frac2', 1 - SPLIT1)]:
    s = re.sub(rf'(<!--)[^>]*(-->\s*<name>{name}</name>\s*<type>Function</type>\s*<expression>).*?(</expression>)',
               lambda m: f'<!-- {name}: {f:g} * rho * q(t), q(t) = Injektionsstufen 40.6 [m3/s] -> [kg/s] -->'
                         + m.group(2)[3:] + expr(f * RHO * Q) + m.group(3), s, count=1, flags=re.S)

# Permeabilitaet aller Kluefte
for i in range(1, 5):
    s = re.sub(rf'(<name>kappa_frac{i}</name>\s*<type>Constant</type>\s*<value>)[^<]*(</value>)',
               rf'\g<1>{K_MEAN:.4e}\g<2>', s)
# Speicherkoeffizient (Media 1-4: storage Constant)
s = re.sub(r'(<name>storage</name>\s*<type>Constant</type>\s*<value>)\s*0\.0\s*(</value>)',
           rf'\g<1>{S_MEAN:.4e}\g<2>', s)
# Apertur
s = re.sub(r'(<name>fracture_thickness_const</name>.*?</parameter>)',
           lambda m: re.sub(r'(<index>[1-4]</index>\s*<value>)[^<]*(</value>)', rf'\g<1>{APERTURE:g}\g<2>', m.group(1)),
           s, count=1, flags=re.S)
# Zeitschrittsteuerung wie 40.6
s = re.sub(r'<time_stepping>.*?</time_stepping>', f'''<time_stepping>
                    <type>FixedTimeStepping</type>
                    <t_initial>{T_START:g}</t_initial>
                    <t_end>{T_END:g}</t_end>
                    <timesteps>
                        <pair>
                            <repeat>{int(round((T_END - T_START) / DT))}</repeat>
                            <delta_t>{DT:g}</delta_t>
                        </pair>
                    </timesteps>
                </time_stepping>''', s, count=1, flags=re.S)
s = s.replace('<prefix>Stimtec_DFN_{:meshname}</prefix>', '<prefix>Stimtec_DFN_40_{:meshname}</prefix>')
s = s.replace('<OpenGeoSysProject>', f'''<OpenGeoSysProject>
    <!-- DFN-Modell Pumpversuch BH10_20180716_40.6 (erzeugt mit make_prj_dfn_40.py)
         k = {K_MEAN:.4e} m2, S = {S_MEAN:.4e} 1/Pa (zeitgew. Mittel V4 40.6, 1750-2670 s),
         Apertur = {APERTURE:g} m, Ratenanteil Kluft 1 = {SPLIT1:g} -->''', 1)
open(dst, 'w', encoding='ISO-8859-1').write(s)
print(f'geschrieben: {dst}  k={K_MEAN:.4e} S={S_MEAN:.4e} a={APERTURE:g} split1={SPLIT1:g}')
