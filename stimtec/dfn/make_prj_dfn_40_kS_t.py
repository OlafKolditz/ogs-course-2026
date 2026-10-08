"""V4 im DFN-Modell: k(t), S(t) gekoppelt an die Injektionsrate (wie 2D-Variante V4 bei 40.6).

Basis: stimtec_dfn_40_a1.prj (Apertur 1 m = 2D-Einheitsdicke, Dirichlet p0 auf Kluftrand,
nodale Quellterme an Knoten 270/271, Stufen 40.6).
X(t) = X_min + (X_max - X_min) * (q(t) - q_min) / (q_max - q_min), gleiche Zeitfenster wie q(t);
gilt fuer alle 4 Kluefte.

Aufruf: python make_prj_dfn_40_kS_t.py [basis.prj] [ziel.prj] [--kmin= --kmax= --smin= --smax=] [--split1=] [--calib]
 --split1: Anteil der Gesamtrate in Kluft 1 (Knoten 270), Rest Kluft 2 (Knoten 271); ohne Angabe 0.5
 --calib: Output nur auf BH10_intersections (schneller, fuer die Kalibrierung)
"""
import os, re, sys
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))

STUFEN = os.path.join(HERE, '..', 'BH10_20180716_40.6', 'workflow', 'injektion_stufen_40.csv')
if not os.path.exists(STUFEN):
    STUFEN = os.path.join(HERE, 'injektion_stufen_40.csv')
_st = pd.read_csv(STUFEN)
T_BOUNDS = _st['t_end_s'].values[:-1]
Q = _st['q_m3s'].values
QMIN, QMAX = Q.min(), Q.max()

# Startwerte = Optimum V4 2D (40.6)
DEFAULT = dict(kmin=6.158e-16, kmax=1.7308e-14, smin=2.0100e-10, smax=7.9833e-10)


def mapped(xmin, xmax):
    return xmin + (xmax - xmin) * (Q - QMIN) / (QMAX - QMIN)


def expr(vals):
    s = '('
    for tb, v in zip(T_BOUNDS, vals[:-1]):
        s += f't &lt; {tb:.1f} ? {v:.8e} : '
    return s + f'{vals[-1]:.8e})'


def build(base_text, kmin, kmax, smin, smax, calib=False, prefix='Stimtec_DFN_40_kS_t', split1=None):
    s = base_text
    if split1 is not None:
        # Ratenaufteilung: Kluft 1 erhaelt split1, Kluft 2 (1 - split1) der Gesamtrate
        for name, f in [('q_frac1', split1), ('q_frac2', 1.0 - split1)]:
            s = re.sub(rf'<!--[^\n]*?-->(\s*<name>{name}</name>\s*<type>Function</type>\s*<expression>).*?(</expression>)',
                       lambda m: f'<!-- {name}: {f:.6f} * rho * q(t), q(t) = Injektionsstufen 40.6 [m3/s] -> [kg/s] -->'
                                 + m.group(1) + expr(f * 1000.0 * Q) + m.group(2), s, count=1, flags=re.S)
    # Permeabilitaet: kappa_frac1..4 -> Verweis auf permeability_t
    s = re.sub(r'<parameter_name>kappa_frac[1-4]</parameter_name>', '<parameter_name>permeability_t</parameter_name>', s)
    # Speicherung: Constant -> Parameter storage_t
    s = re.sub(r'<name>storage</name>\s*<type>Constant</type>\s*<value>[^<]*</value>',
               '<name>storage</name>\n                    <type>Parameter</type>\n                    <parameter_name>storage_t</parameter_name>', s)
    pars = f'''        <parameter>
            <!-- k(t): Zeitfenster wie Injektion 40.6, k_min={kmin:.4e}, k_max={kmax:.4e} m2 -->
            <name>permeability_t</name>
            <type>Function</type>
            <expression>{expr(mapped(kmin, kmax))}</expression>
        </parameter>
        <parameter>
            <!-- S(t): Zeitfenster wie Injektion 40.6, S_min={smin:.4e}, S_max={smax:.4e} 1/Pa -->
            <name>storage_t</name>
            <type>Function</type>
            <expression>{expr(mapped(smin, smax))}</expression>
        </parameter>
'''
    s = re.sub(r'\s*<parameter>\s*<!-- k\(t\).*?</parameter>\s*<parameter>\s*<!-- S\(t\).*?</parameter>', '', s, flags=re.S)
    s = s.replace('<parameters>', '<parameters>\n' + pars.rstrip('\n'), 1)
    s = re.sub(r'<prefix>[^<]*</prefix>', f'<prefix>{prefix}_{{:meshname}}</prefix>', s, count=1)
    if calib:
        s = re.sub(r'<meshes>\s*<mesh>mesh</mesh>\s*<mesh>BH10_intersections</mesh>\s*</meshes>',
                   '<meshes>\n                <mesh>BH10_intersections</mesh>\n            </meshes>', s)
        s = s.replace('<prefix>', '<prefix>', 1)
    s = re.sub(r'(<!-- DFN-Modell Pumpversuch.*?)-->',
               lambda m: m.group(1) + f'\n         V4: k(t), S(t) an Injektionsrate gekoppelt (make_prj_dfn_40_kS_t.py) -->', s, count=1, flags=re.S)
    return s


if __name__ == '__main__':
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    opts = dict(a[2:].split('=') for a in sys.argv[1:] if a.startswith('--') and '=' in a)
    src = args[0] if args else os.path.join(HERE, 'stimtec_dfn_40_a1.prj')
    dst = args[1] if len(args) > 1 else os.path.join(HERE, 'stimtec_dfn_40_kS_t.prj')
    p = {k: float(opts.get(k, v)) for k, v in DEFAULT.items()}
    if 'split1' in opts:
        p['split1'] = float(opts['split1'])
    txt = build(open(src, encoding='ISO-8859-1').read(), **p, calib='--calib' in sys.argv)
    open(dst, 'w', encoding='ISO-8859-1').write(txt)
    print('geschrieben:', dst, p)
