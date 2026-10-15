"""
Comparison of the OGS simulation (results/*.pvd) with the measured pressure
(workflow/bh10_40_gemessen_digitalisiert.csv) - pumping test BH10_20180718_40.6.

 - Simulation: pressure at the mesh node closest to the observation point
   (default: well (0,0,0)), relative to the initial pressure (first time step),
   so that p0 = 0 (V1) and p0 = 1 MPa (V2..V5) are treated the same way.
 - Measurement: pressure "P int d" [MPa], baseline = median t < 1780 s subtracted (as in bh40.py).
 - Comparison at the simulation output times (measurement linearly interpolated),
   i.e. identical to the objective function of the calibration scripts.
 - Metrics: RMSE, MAE, bias, max |error|, NSE, Pearson r - total and per injection stage.

Usage (in the workflow folder):
    python compare_sim_meas_en.py                       # all *.pvd in ../results
    python compare_sim_meas_en.py ../results/stimtec_bc_rightmid_py.pvd
    python compare_sim_meas_en.py --point 0 0 0 --tmin 1785 --tmax 2670
    python compare_sim_meas_en.py --pref 1e6            # fixed reference pressure [Pa]
    python compare_sim_meas_en.py --no-show             # save only, no windows (e.g. batch)

Output (in the workflow folder):
    comparison_<name>.csv       time series time_s, p_sim_Pa, p_meas_Pa, residual_Pa
    comparison_metrics.csv      metrics per simulation (total + per stage)
    comparison_sim_meas.png     pressure curves + residuals
    comparison_stages.png       RMSE and bias per injection stage
The figures are also shown on screen (plt.show()).
"""
import argparse, glob, os, sys
import xml.etree.ElementTree as ET
import numpy as np, pandas as pd
import matplotlib
if "--no-show" in sys.argv:          # batch mode without display
    matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pyvista as pv

HERE = os.path.dirname(os.path.abspath(__file__))
RES_DIR = os.path.normpath(os.path.join(HERE, "..", "results"))
MEAS_FILE = os.path.join(HERE, "bh10_40_gemessen_digitalisiert.csv")
STEPS_FILE = os.path.join(HERE, "injektion_stufen_40.csv")
T_BASE = 1780.0
T_TOL = 1.0   # [s] tolerance at the ends of the measured series (measurement starts at 1750.009 s)
COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#c9a227", "#8e5bd0", "#d64545", "#9a9892"]


# ------------------------------------------------------------ measurement
def load_meas():
    """Time [s], pressure build-up [Pa] (baseline = median t < T_BASE subtracted)."""
    d = pd.read_csv(MEAS_FILE).sort_values("time_s")
    d = d.drop_duplicates("time_s")
    base = d.loc[d["time_s"] < T_BASE, "value_axis_units"].median()
    return d["time_s"].values, (d["value_axis_units"].values - base) * 1e6


def load_steps():
    """Injection stages (t_end_s, q_l_min) or None."""
    if not os.path.exists(STEPS_FILE):
        return None
    st = pd.read_csv(STEPS_FILE)
    qcol = next((c for c in st.columns if c.lower() in ("q_lpm", "q_l_min", "q_lmin")), None)
    return st["t_end_s"].values.astype(float), (st[qcol].values if qcol else None)


# ------------------------------------------------------------- simulation
def read_pvd(pvd, point, pref):
    """Time [s], pressure build-up [Pa] at the observation point from an OGS pvd series."""
    base_dir = os.path.dirname(os.path.abspath(pvd))
    ents = sorted((float(ds.attrib["timestep"]), os.path.normpath(os.path.join(base_dir, ds.attrib["file"])))
                  for ds in ET.parse(pvd).getroot().iter("DataSet"))
    t, p, node, dist, missing = [], [], None, None, 0
    for ti, fn in ents:
        if not os.path.exists(fn):
            missing += 1
            continue
        m = pv.read(fn)
        if node is None:
            d = np.linalg.norm(m.points - point, axis=1)
            node, dist = int(np.argmin(d)), float(d.min())
        t.append(ti)
        p.append(float(m.point_data["pressure"][node]))
    if not t:
        raise FileNotFoundError(f"no vtu files found for {pvd}")
    if missing:
        print(f"  Warning: {missing} of {len(ents)} vtu files are missing")
    t, p = np.array(t), np.array(p)
    p_ref = p[0] if pref == "initial" else float(pref)
    return t, p - p_ref, node, dist, p_ref


# ---------------------------------------------------------------- metrics
def metrics(sim, meas):
    e = sim - meas
    sst = np.sum((meas - meas.mean()) ** 2)
    return {
        "n": len(e),
        "RMSE_MPa": np.sqrt(np.mean(e ** 2)) / 1e6,
        "MAE_MPa": np.mean(np.abs(e)) / 1e6,
        "Bias_MPa": np.mean(e) / 1e6,
        "MaxAbs_MPa": np.max(np.abs(e)) / 1e6,
        "NSE": 1 - np.sum(e ** 2) / sst if sst > 0 else np.nan,
        "r": np.corrcoef(sim, meas)[0, 1] if len(e) > 1 and np.std(sim) > 0 else np.nan,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("pvd", nargs="*", help="pvd files (default: all in ../results)")
    ap.add_argument("--point", nargs=3, type=float, default=[0.0, 0.0, 0.0], metavar=("X", "Y", "Z"),
                    help="observation point [m] (default: well 0 0 0)")
    ap.add_argument("--pref", default="initial",
                    help="reference pressure [Pa] or 'initial' (= pressure at first time step, default)")
    ap.add_argument("--tmin", type=float, default=None, help="start of evaluation window [s]")
    ap.add_argument("--tmax", type=float, default=None, help="end of evaluation window [s]")
    ap.add_argument("--out", default="comparison_sim_meas.png", help="name of the figure")
    ap.add_argument("--no-show", action="store_true", help="save figures only, do not display them")
    a = ap.parse_args()

    pvds = a.pvd or sorted(glob.glob(os.path.join(RES_DIR, "*.pvd")))
    if not pvds:
        sys.exit(f"no pvd files in {RES_DIR}")
    point = np.array(a.point)
    MT, MP = load_meas()
    steps = load_steps()

    fig, (ax, axr) = plt.subplots(2, 1, figsize=(10, 7.5), sharex=True,
                                  gridspec_kw={"height_ratios": [3, 1.3]}, facecolor="#fcfcfb")
    ax.plot(MT, MP / 1e6, color="#52514e", lw=1.0, label="Measurement 40.6 (baseline subtracted)")

    rows = []
    for i, pvd in enumerate(pvds):
        name = os.path.splitext(os.path.basename(pvd))[0]
        print(f"{name}:")
        t, p, node, dist, p_ref = read_pvd(pvd, point, a.pref)
        sel = (t >= MT.min() - T_TOL) & (t <= MT.max() + T_TOL)
        if a.tmin is not None: sel &= t >= a.tmin
        if a.tmax is not None: sel &= t <= a.tmax
        ts, ps = t[sel], p[sel]
        if len(ts) < 2:
            print("  skipped: fewer than 2 time points in the evaluation window")
            continue
        pm = np.interp(ts, MT, MP)
        pd.DataFrame({"time_s": ts, "p_sim_Pa": ps, "p_meas_Pa": pm, "residual_Pa": ps - pm}) \
          .to_csv(os.path.join(HERE, f"comparison_{name}.csv"), index=False)

        m = metrics(ps, pm)
        rows.append({"simulation": name, "stage": "total", "t_from": ts.min(), "t_to": ts.max(), **m})
        print(f"  node {node} (distance {dist:.3g} m), p_ref = {p_ref:.4g} Pa, {len(ts)} time points")
        print("  RMSE {RMSE_MPa:.4f} MPa | MAE {MAE_MPa:.4f} | Bias {Bias_MPa:+.4f} | "
              "max|e| {MaxAbs_MPa:.4f} MPa | NSE {NSE:.4f} | r {r:.4f}".format(**m))

        if steps is not None:  # metrics per injection stage
            t_end, q = steps
            t0 = ts.min()
            for k, te in enumerate(t_end):
                s = (ts >= t0) & (ts < te) if k < len(t_end) - 1 else (ts >= t0)
                if s.sum() >= 2:
                    lab = f"{q[k]:.3f} l/min" if q is not None else f"stage {k + 1}"
                    rows.append({"simulation": name, "stage": lab, "t_from": t0, "t_to": te,
                                 **metrics(ps[s], pm[s])})
                t0 = te

        col = COLORS[i % len(COLORS)]
        ax.plot(t, p / 1e6, color=col, lw=2, label=f"{name}: RMSE {m['RMSE_MPa']:.3f} MPa, NSE {m['NSE']:.3f}")
        axr.plot(ts, (ps - pm) / 1e6, color=col, lw=1.5)

    res = pd.DataFrame(rows)
    res.to_csv(os.path.join(HERE, "comparison_metrics.csv"), index=False, float_format="%.6g")
    if steps is not None and not res.empty:
        print("\nMetrics per stage:")
        with pd.option_context("display.width", 140, "display.float_format", "{:.4f}".format):
            print(res[["simulation", "stage", "t_from", "t_to", "RMSE_MPa", "Bias_MPa", "MaxAbs_MPa"]]
                  .to_string(index=False))

    if steps is not None:  # mark stage boundaries
        for te in steps[0][:-1]:
            for x in (ax, axr):
                x.axvline(te, color="#c8c6c0", lw=0.8, ls="--", zorder=0)
    axr.axhline(0, color="#52514e", lw=0.8)
    ax.set_ylabel("Pressure build-up [MPa]")
    axr.set_ylabel("Sim − Meas [MPa]")
    axr.set_xlabel("Time [s]")
    ax.set_title(f"BH10_20180718_40.6 – simulation vs. measurement at point "
                 f"({point[0]:g}, {point[1]:g}, {point[2]:g})", loc="left")
    for x in (ax, axr):
        x.set_facecolor("#fcfcfb"); x.grid(alpha=0.25)
        for s in ("top", "right"):
            x.spines[s].set_visible(False)
    ax.legend(frameon=False, loc="upper left", fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, a.out), dpi=150)
    print(f"\nFigure: {a.out}\nMetrics: comparison_metrics.csv")

    if steps is not None and not res.empty:
        plot_stages(res, os.path.join(HERE, "comparison_stages.png"))
        print("Figure: comparison_stages.png")
    if not a.no_show:
        plt.show()


def plot_stages(res, out):
    """Bar chart of RMSE and bias per injection stage for all simulations."""
    st = res[res["stage"] != "total"]
    sims = list(dict.fromkeys(st["simulation"]))
    labels = list(dict.fromkeys(st["stage"]))
    x = np.arange(len(labels))
    w = 0.8 / max(len(sims), 1)
    fig, axs = plt.subplots(2, 1, figsize=(10, 6.5), sharex=True, facecolor="#fcfcfb")
    for j, sim in enumerate(sims):
        d = st[st["simulation"] == sim].set_index("stage").reindex(labels)
        col = COLORS[j % len(COLORS)]
        axs[0].bar(x + (j - (len(sims) - 1) / 2) * w, d["RMSE_MPa"], w, color=col, label=sim)
        axs[1].bar(x + (j - (len(sims) - 1) / 2) * w, d["Bias_MPa"], w, color=col)
    axs[0].set_ylabel("RMSE [MPa]")
    axs[1].set_ylabel("Bias Sim − Meas [MPa]")
    axs[1].axhline(0, color="#52514e", lw=0.8)
    axs[1].set_xticks(x, labels)
    axs[1].set_xlabel("Injection stage (in chronological order)")
    axs[0].set_title("Deviation per injection stage", loc="left")
    axs[0].legend(frameon=False, fontsize=9)
    for ax in axs:
        ax.set_facecolor("#fcfcfb"); ax.grid(axis="y", alpha=0.25)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(out, dpi=150)


if __name__ == "__main__":
    main()
