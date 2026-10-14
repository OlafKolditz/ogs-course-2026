"""Kalibrierung DFN-V4 + Ratenaufteilung (Pumpversuch 40.6): k_min, k_max, S_min, S_max (log10) und
split1 = Anteil der Injektionsrate in Kluft 1 (Knoten 270), Rest in Kluft 2 (Knoten 271); split1 ueber logit.

Zielfunktion (OBJ, Standard 'nodes'): RMSE ueber beide Knoten einzeln gegen die Messung,
  sqrt(mean((p270 - p_mess)^2 + (p271 - p_mess)^2)/2)  -> beide Schnittpunkte liegen im selben
  Packerintervall und muessen denselben Druck zeigen. Alternative OBJ=mean: RMSE des Knotenmittels.
Aufruf: python calibrate_dfn_kS_split.py   (Umgebung: X0="log10 kmin,kmax,smin,smax,split1", MAXFEV, OBJ)
"""
import os, sys, time
import numpy as np
from scipy.optimize import minimize
import calibrate_dfn_kS_t as C

HERE = os.path.dirname(os.path.abspath(__file__))
OBJ = os.environ.get('OBJ', 'nodes')
LOG = os.path.join(HERE, f'calibrate_dfn_kS_split_{OBJ}_log.csv')


def errors(a):
    p1 = np.interp(C.TM, a[:, 0], a[:, 1]); p2 = np.interp(C.TM, a[:, 0], a[:, 2])
    r_nodes = np.sqrt(0.5 * np.mean((p1 - C.PM) ** 2 + (p2 - C.PM) ** 2))
    r_mean = np.sqrt(np.mean((0.5 * (p1 + p2) - C.PM) ** 2))
    r_diff = np.sqrt(np.mean((p1 - p2) ** 2))
    return r_nodes, r_mean, r_diff


def run(kmin, kmax, smin, smax, split1):
    C.BUILD_KW = dict(split1=split1)
    return C.simulate(kmin, kmax, smin, smax)


n = [0]
def obj(x):
    kmin, kmax, smin, smax = 10.0 ** np.asarray(x[:4]); split1 = 1 / (1 + np.exp(-x[4]))
    a = run(kmin, kmax, smin, smax, split1)
    rn, rm, rd = (1e9, 1e9, 1e9) if a is None else errors(a)
    f = rn if OBJ == 'nodes' else rm
    n[0] += 1
    with open(LOG, 'a') as fh:
        fh.write(f'{n[0]},{kmin:.6e},{kmax:.6e},{smin:.6e},{smax:.6e},{split1:.6f},{rn:.6e},{rm:.6e},{rd:.6e}\n')
    print(f'{n[0]:4d} k=[{kmin:.3e},{kmax:.3e}] S=[{smin:.3e},{smax:.3e}] split1={split1:.3f} '
          f'RMSE_nodes={rn/1e6:.4f} RMSE_mean={rm/1e6:.4f} |p1-p2|={rd/1e6:.4f} MPa', flush=True)
    return f


if __name__ == '__main__':
    if not os.path.exists(LOG):
        open(LOG, 'w').write('run,k_min,k_max,S_min,S_max,split1,rmse_nodes_Pa,rmse_mean_Pa,rms_p1_p2_Pa\n')
    if len(sys.argv) > 1 and sys.argv[1] == 'scan':        # Sensitivitaet split1 bei festen k, S
        x0 = [float(v) for v in os.environ['X0'].split(',')][:4]
        for s1 in [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]:
            obj(np.r_[x0, np.log(s1 / (1 - s1))])
        sys.exit()
    v = [float(v) for v in os.environ['X0'].split(',')]
    x0 = np.r_[v[:4], np.log(v[4] / (1 - v[4]))]
    steps = [0.3, 0.3, 0.3, 0.3, 0.6]
    simplex = np.vstack([x0] + [x0 + steps[i] * np.eye(5)[i] for i in range(5)])
    t0 = time.time()
    res = minimize(obj, x0, method='Nelder-Mead',
                   options=dict(initial_simplex=simplex, maxfev=int(os.environ.get('MAXFEV', 300)), xatol=0.005, fatol=100.0))
    kmin, kmax, smin, smax = 10 ** res.x[:4]; s1 = 1 / (1 + np.exp(-res.x[4]))
    print(f'Optimum ({OBJ}): k_min={kmin:.4e} k_max={kmax:.4e} S_min={smin:.4e} S_max={smax:.4e} split1={s1:.4f} '
          f'RMSE={res.fun/1e6:.4f} MPa ({res.nfev} Laeufe, {time.time()-t0:.0f} s)')
    open(os.path.join(HERE, f'calibrate_dfn_kS_split_{OBJ}_best.txt'), 'w').write(
        f'k_min={kmin:.6e}\nk_max={kmax:.6e}\nS_min={smin:.6e}\nS_max={smax:.6e}\nsplit1={s1:.6f}\nrmse_MPa={res.fun/1e6:.6f}\nnfev={res.nfev}\n')
