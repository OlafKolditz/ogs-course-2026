"""Kalibrierung V4 im DFN-Modell (Pumpversuch 40.6): k_min, k_max, S_min, S_max (log10, Nelder-Mead).
Zielfunktion: RMSE zwischen gemessenem Druckaufbau (Baseline = Median t < 1780 s) und dem
mittleren simulierten Druckaufbau an den Schnittpunkten BH10 x Kluft 1/2 (Knoten 270, 271), 1750-2670 s.
Aufruf: python calibrate_dfn_kS_t.py   (Umgebung: X0="log10 kmin,kmax,smin,smax", MAXFEV=250, OGS=ogs)
"""
import os, re, glob, shutil, subprocess, tempfile, time
import numpy as np, pandas as pd
from scipy.optimize import minimize
import vtk
from vtk.util.numpy_support import vtk_to_numpy as v2n
from make_prj_dfn_40_kS_t import build, DEFAULT

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = open(os.path.join(HERE, 'stimtec_dfn_40_a1.prj'), encoding='ISO-8859-1').read()
MESHES = ['mesh.vtu', 'BH10_frac1.vtu', 'BH10_frac2.vtu', 'BH10_intersections.vtu', 'bc_rim.vtu']
OGS = os.environ.get('OGS', 'ogs')
MEAS = os.path.join(HERE, '..', 'BH10_20180716_40.6', 'workflow', 'bh10_40_gemessen_digitalisiert.csv')
if not os.path.exists(MEAS):
    MEAS = os.path.join(HERE, 'bh10_40_gemessen_digitalisiert.csv')
m = pd.read_csv(MEAS)
base = m.loc[m.time_s < 1780, 'value_axis_units'].median()
TM, PM = m.time_s.values, (m.value_axis_units.values - base) * 1e6
LOG = os.path.join(HERE, 'calibrate_dfn_kS_t_log.csv')
BUILD_KW = {}   # z. B. dict(split1=0.4), gesetzt von calibrate_dfn_kS_split.py


def tt(f):
    return float(re.search(r'_t_([0-9.e+-]+)\.vtu', f).group(1))


def simulate(kmin, kmax, smin, smax):
    d = tempfile.mkdtemp(prefix='dfn_')
    try:
        for f in MESHES:
            os.symlink(os.path.join(HERE, f), os.path.join(d, f)) if os.name != 'nt' else shutil.copy(os.path.join(HERE, f), d)
        open(os.path.join(d, 'run.prj'), 'w', encoding='ISO-8859-1').write(build(BASE, kmin, kmax, smin, smax, calib=True, **BUILD_KW))
        r = subprocess.run([OGS, 'run.prj', '-o', 'out'], cwd=d, capture_output=True, text=True)
        fs = sorted(glob.glob(os.path.join(d, 'out', '*BH10_intersections_ts_*.vtu')), key=tt)
        if r.returncode != 0 or tt(fs[-1]) < 2669:
            return None
        rows = []
        for f in fs:
            rd = vtk.vtkXMLUnstructuredGridReader(); rd.SetFileName(f); rd.Update()
            rows.append([tt(f)] + list(v2n(rd.GetOutput().GetPointData().GetArray('pressure'))))
        a = np.array(rows); a[:, 1:] -= a[0, 1:]
        return a
    finally:
        shutil.rmtree(d, ignore_errors=True)


def rmse(a):
    return np.sqrt(np.mean((np.interp(TM, a[:, 0], a[:, 1:].mean(1)) - PM) ** 2))


n = [0]
def obj(x):
    kmin, kmax, smin, smax = 10.0 ** np.asarray(x)
    a = simulate(kmin, kmax, smin, smax)
    f = 1e9 if a is None else rmse(a)
    n[0] += 1
    with open(LOG, 'a') as fh:
        fh.write(f'{n[0]},{kmin:.6e},{kmax:.6e},{smin:.6e},{smax:.6e},{f:.6e}\n')
    print(f'{n[0]:4d} k=[{kmin:.3e},{kmax:.3e}] S=[{smin:.3e},{smax:.3e}] RMSE={f/1e6:.4f} MPa', flush=True)
    return f


if __name__ == '__main__':
    if not os.path.exists(LOG):
        open(LOG, 'w').write('run,k_min,k_max,S_min,S_max,rmse_Pa\n')
    x0 = os.environ.get('X0')
    x0 = np.array([float(v) for v in x0.split(',')]) if x0 else np.log10([DEFAULT[k] for k in ('kmin', 'kmax', 'smin', 'smax')])
    simplex = np.vstack([x0] + [x0 + 0.5 * np.eye(4)[i] for i in range(4)])
    t0 = time.time()
    res = minimize(obj, x0, method='Nelder-Mead',
                   options=dict(initial_simplex=simplex, maxfev=int(os.environ.get('MAXFEV', 250)), xatol=0.005, fatol=100.0))
    kmin, kmax, smin, smax = 10 ** res.x
    print(f'Optimum: k_min={kmin:.4e} k_max={kmax:.4e} S_min={smin:.4e} S_max={smax:.4e} RMSE={res.fun/1e6:.4f} MPa '
          f'({res.nfev} Laeufe, {time.time()-t0:.0f} s)')
    open(os.path.join(HERE, 'calibrate_dfn_kS_t_best.txt'), 'w').write(
        f'k_min={kmin:.6e}\nk_max={kmax:.6e}\nS_min={smin:.6e}\nS_max={smax:.6e}\nrmse_MPa={res.fun/1e6:.6f}\nnfev={res.nfev}\n')
