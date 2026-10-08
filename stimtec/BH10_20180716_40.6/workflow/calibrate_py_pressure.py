"""
Kalibrierung des druckabhaengigen Modells "pressure_exp" (ogs_py_kS.py) ueber
die OGS-Python-Schnittstelle:
    k(p) = k0 * exp(beta_k * dp),  S(p) = S0 * exp(beta_S * dp),  dp = p - 1 MPa je Zelle
Parameter: log10 k0, beta_k [1/MPa], log10 S0, beta_S [1/MPa]
"""
import os
import numpy as np, pandas as pd
from scipy.optimize import minimize
import ogs_py_kS as o

HERE = o.HERE
import bh40
MT, MP = bh40.load_meas()   # 40.6
LOG = []


def par_of(x):
    return dict(k0=10 ** x[0], beta_k=x[1] * 1e-6, k_cap=1e-12,
                S0=10 ** x[2], beta_S=x[3] * 1e-6, S_cap=1e-9)


def obj(x):
    try:
        t, dp, k, S = o.run("pressure_exp", par_of(x))
        err = float(np.sqrt(np.mean((dp - np.interp(t, MT, MP)) ** 2)))
    except RuntimeError:
        err = 1e12
    LOG.append((*x, err))
    print(f"log10k0={x[0]:.3f} beta_k={x[1]:.3f}/MPa log10S0={x[2]:.3f} beta_S={x[3]:.3f}/MPa RMSE={err:.4e}", flush=True)
    return err


def main():
    # Startwert aus Variante 4: k steigt ~Faktor 7 bis ~7 MPa -> beta_k ~ 0.3/MPa
    x0 = np.array([np.log10(3.4e-15), 0.3, np.log10(2.5e-14), 0.8])
    simplex = [x0] + [x0 + d for d in np.diag([0.3, 0.1, 0.5, 0.2])]
    res = minimize(obj, x0, method="Nelder-Mead",
                   options={"xatol": 1e-3, "fatol": 10.0, "maxfev": 350, "initial_simplex": simplex})
    t, dp, k, S = o.run("pressure_exp", par_of(res.x))
    err = float(np.sqrt(np.mean((dp - np.interp(t, MT, MP)) ** 2)))
    pd.DataFrame(LOG, columns=["log10_k0", "beta_k_perMPa", "log10_S0", "beta_S_perMPa", "rmse_Pa"]).to_csv(
        os.path.join(HERE, "calibrate_py_pressure_log.csv"), index=False)
    pd.DataFrame({"time_s": t, "p_sim_Pa": dp, "p_meas_Pa": np.interp(t, MT, MP), "k_mean": k, "S_mean": S}).to_csv(
        os.path.join(HERE, "calibrate_py_pressure_best.csv"), index=False)
    p = par_of(res.x)
    print(f"\nOPTIMUM: k0={p['k0']:.4e} m2, beta_k={res.x[1]:.4f} 1/MPa, "
          f"S0={p['S0']:.4e} 1/Pa, beta_S={res.x[3]:.4f} 1/MPa, RMSE={err:.4e} Pa")


if __name__ == "__main__":
    main()
