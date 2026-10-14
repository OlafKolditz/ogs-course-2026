#!/usr/bin/env python3
"""
compare_dfn_injection.py  -  STIMTEC DFN: pressure evolution at one injection point

Compares the pressure history of all (or selected) OGS runs in ./results at a
single BH10 / DFN intersection node. Only results that actually exist are used;
missing projects are skipped with a note.

Results layout accepted (per project, inside the results folder):
    results/<project>/            extracted OGS output (*.pvd + *.vtu)
    results/<project>.zip         zipped OGS output (read on the fly, not extracted permanently)

Source mesh per project (first one found is used):
    *_BH10_intersections.pvd      small boundary mesh, point picked via bulk_node_ids
    *_mesh.pvd / any other *.pvd  full mesh, point picked by node id directly

Injection points (BH10_intersections.csv):
    node 270  ->  fracture 1  (MaterialID 1, 5.4 m from BH10 deep end)
    node 271  ->  fracture 2  (MaterialID 2, 15.5 m from BH10 deep end)

Examples
    python compare_dfn_injection.py                         # all projects, node 271
    python compare_dfn_injection.py --node 270
    python compare_dfn_injection.py -p stimtec_dfn_40 stimtec_dfn_40_kS_t stimtec_dfn_40_kS_split
    python compare_dfn_injection.py -p 40_kS               # substring match also works
    python compare_dfn_injection.py --measured other_test.csv --baseline-end 1780
    python compare_dfn_injection.py --no-measured --absolute  # simulation only, absolute pressure
    python compare_dfn_injection.py --ylim -1 6             # zoom (plot units)
    python compare_dfn_injection.py --save cmp_node271.png --csv cmp_node271.csv --no-show

Measured data
    By default the digitised pressure of test BH10_20180716_40.6 is added:
        ../BH10_20180716_40.6/workflow/bh10_40_gemessen_digitalisiert.csv  (or ./ as fallback)
    columns time_s, value_axis_units [MPa]. As in calibrate_dfn_kS_t.py, the baseline
    (median for t < 1780 s) is subtracted from the measurement and the simulation is shown
    as pressure build-up p(t) - p(t0). The RMSE simulation vs. measurement (on the
    overlapping time window) is printed and shown in the legend.

Requires: numpy, matplotlib, pyvista   (pip install numpy matplotlib pyvista)
"""
import argparse
import glob
import os
import sys
import tempfile
import xml.etree.ElementTree as ET
import zipfile

import numpy as np
import matplotlib.pyplot as plt
import pyvista as pv

HERE = os.path.dirname(os.path.abspath(__file__))
INJECTION_NODES = {270: "Kluft 1 (node 270)", 271: "Kluft 2 (node 271)"}
UNITS = {"Pa": 1.0, "kPa": 1e3, "bar": 1e5, "MPa": 1e6}
MEAS_DEFAULT = [
    os.path.join(HERE, "..", "BH10_20180716_40.6", "workflow", "bh10_40_gemessen_digitalisiert.csv"),
    os.path.join(HERE, "bh10_40_gemessen_digitalisiert.csv"),
]


# --------------------------------------------------------------------------- #
def read_pvd(pvd_path):
    """Return list of (time, vtu_path) from an OGS .pvd collection."""
    root = ET.parse(pvd_path).getroot()
    base = os.path.dirname(pvd_path)
    out = []
    for ds in root.iter("DataSet"):
        f = os.path.join(base, ds.get("file"))
        if os.path.isfile(f):                      # skip time steps that are missing
            out.append((float(ds.get("timestep")), f))
    return sorted(out)


def pick_pvd(folder):
    """Prefer the intersection boundary mesh, else the full mesh."""
    pvds = sorted(glob.glob(os.path.join(folder, "*.pvd")))
    if not pvds:
        return None
    for key in ("_BH10_intersections.pvd", "_mesh.pvd"):
        hits = [p for p in pvds if p.endswith(key)]
        if hits:
            return hits[0]
    return pvds[0]


def point_index(mesh, node, is_bulk):
    if is_bulk:
        ids = np.asarray(mesh.point_data["bulk_node_ids"])
        hit = np.where(ids == node)[0]
        if hit.size == 0:
            raise KeyError(f"node {node} not in bulk_node_ids {ids.tolist()}")
        return int(hit[0])
    if node >= mesh.n_points:
        raise KeyError(f"node {node} >= n_points {mesh.n_points}")
    return node


def pressure_history(folder, node, var="pressure"):
    pvd = pick_pvd(folder)
    if pvd is None:
        raise FileNotFoundError("no .pvd found")
    steps = read_pvd(pvd)
    if not steps:
        raise FileNotFoundError(f"{os.path.basename(pvd)}: no vtu files present")
    first = pv.read(steps[0][1])
    is_bulk = "bulk_node_ids" in first.point_data
    idx = point_index(first, node, is_bulk)
    t, p = [], []
    for ti, f in steps:
        m = first if f == steps[0][1] else pv.read(f)
        t.append(ti)
        p.append(float(np.asarray(m.point_data[var])[idx]))
    return np.array(t), np.array(p), os.path.basename(pvd)


# --------------------------------------------------------------------------- #
def find_projects(results_dir):
    """Map project name -> ('dir'|'zip', path) for everything that exists."""
    found = {}
    for entry in sorted(os.listdir(results_dir)):
        path = os.path.join(results_dir, entry)
        if os.path.isdir(path):
            # also accept nested folders, e.g. results/<p>/<p>/*.pvd from "Extract All"
            pvds = sorted(glob.glob(os.path.join(path, "**", "*.pvd"), recursive=True),
                          key=lambda f: f.count(os.sep))
            if pvds:
                found[entry] = ("dir", os.path.dirname(pvds[0]))
        elif entry.lower().endswith(".zip"):
            found.setdefault(entry[:-4], ("zip", path))   # folder wins over zip
    return found


def select(found, wanted):
    if not wanted:
        return list(found)
    sel = []
    for w in wanted:
        if w in found:
            sel.append(w)
        else:
            hits = [k for k in found if w in k]
            if not hits:
                print(f"  [skip] '{w}': no results found")
            sel += [h for h in hits if h not in sel]
    return sel


def read_measured(path):
    """CSV with time [s] in col 0 and pressure [MPa] in col 1 (header allowed, ',' or ';')."""
    import csv
    with open(path, newline="") as fh:
        sample = fh.read(2048)
        fh.seek(0)
        delim = ";" if sample.count(";") > sample.count(",") else ","
        rows = list(csv.reader(fh, delimiter=delim))
    t, p, header = [], [], None
    for r in rows:
        try:
            t.append(float(r[0])); p.append(float(r[1]))
        except (ValueError, IndexError):
            header = header or r
    t, p = np.array(t), np.array(p)
    o = np.argsort(t, kind="stable")
    return t[o], p[o], header


def rmse(t_sim, p_sim, t_meas, p_meas):
    """RMSE on the overlapping time window (simulation interpolated to measured times)."""
    s = (t_meas >= t_sim[0]) & (t_meas <= t_sim[-1])
    if s.sum() < 2:
        return None
    return float(np.sqrt(np.mean((np.interp(t_meas[s], t_sim, p_sim) - p_meas[s]) ** 2)))


# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-r", "--results", default=os.path.join(HERE, "results"),
                    help="results folder (default: ./results next to this script)")
    ap.add_argument("-n", "--node", type=int, default=271,
                    help="bulk mesh node id of the injection point (270 or 271; default 271)")
    ap.add_argument("-p", "--projects", nargs="*",
                    help="project names or substrings; default: all existing results")
    ap.add_argument("--var", default="pressure", help="point variable (default: pressure)")
    ap.add_argument("--unit", default="MPa", choices=UNITS, help="pressure unit for plot")
    ap.add_argument("--time-unit", default="s", choices=["s", "min", "h"],
                    help="time unit for plot (OGS times are taken as seconds)")
    ap.add_argument("--relative", action="store_true",
                    help="plot pressure change p(t) - p(t0) per run "
                         "(default when measured data are shown)")
    ap.add_argument("--absolute", action="store_true",
                    help="plot absolute simulated pressure even with measured data")
    ap.add_argument("--measured", help="CSV with measured data (time_s, pressure in MPa); "
                                       "default: digitised test BH10_20180716_40.6")
    ap.add_argument("--no-measured", action="store_true", help="do not plot measured data")
    ap.add_argument("--baseline-end", type=float, default=1780.0,
                    help="measured baseline = median for t < this time [s] (default 1780); "
                         "use a negative value to keep the raw measurement")
    ap.add_argument("--ylim", nargs=2, type=float, metavar=("YMIN", "YMAX"),
                    help="y-axis limits in plot units (useful if one run diverges)")
    ap.add_argument("--show-diverged", action="store_true",
                    help="also plot runs with |p| > 1 GPa (hidden by default)")
    ap.add_argument("--save", help="save figure (png/pdf/svg)")
    ap.add_argument("--csv", help="write all extracted curves to one CSV (long format)")
    ap.add_argument("--no-show", action="store_true", help="do not open a plot window")
    a = ap.parse_args()

    if not os.path.isdir(a.results):
        sys.exit(f"results folder not found: {a.results}")
    found = find_projects(a.results)
    if not found:
        sys.exit(f"no OGS results (folders with .pvd or .zip) in {a.results}")
    names = select(found, a.projects)
    if not names:
        sys.exit("none of the requested projects has results")

    # ---- measured data ---------------------------------------------------- #
    meas = None
    if not a.no_measured:
        mpath = a.measured or next((f for f in MEAS_DEFAULT if os.path.isfile(f)), None)
        if mpath and os.path.isfile(mpath):
            tm, pm, _ = read_measured(mpath)
            pm = pm * 1e6                                  # MPa -> Pa
            if a.baseline_end >= 0 and np.any(tm < a.baseline_end):
                pm = pm - np.median(pm[tm < a.baseline_end])
            meas = (tm, pm, os.path.basename(mpath))
            print(f"  [meas] {mpath}: {len(tm)} points, t={tm[0]:.0f}..{tm[-1]:.0f} s")
        else:
            print(f"  [meas] no measured data found ({a.measured or 'default path'}) - skipped")
    relative = (a.relative or meas is not None) and not a.absolute

    tscale = {"s": 1.0, "min": 60.0, "h": 3600.0}[a.time_unit]
    pscale = UNITS[a.unit]
    curves = {}
    with tempfile.TemporaryDirectory() as tmp:
        for name in names:
            kind, path = found[name]
            folder = path
            if kind == "zip":
                folder = os.path.join(tmp, name)
                with zipfile.ZipFile(path) as z:
                    z.extractall(folder)
                sub = [d for d in glob.glob(os.path.join(folder, "**", "*.pvd"), recursive=True)]
                if sub:
                    folder = os.path.dirname(sub[0])
            try:
                t, p, src = pressure_history(folder, a.node, a.var)
            except Exception as e:                      # noqa: BLE001
                print(f"  [skip] {name}: {e}")
                continue
            if np.nanmax(np.abs(p)) > 1e9 and not a.show_diverged:   # > 1 GPa: not physical
                print(f"  [warn] {name}: |p| > 1 GPa - run looks diverged, not plotted "
                      f"(use --show-diverged to include it)")
                continue
            curves[name] = (t, p)
            print(f"  [ok]   {name:28s} {len(t):4d} steps  t={t[0]:.0f}..{t[-1]:.0f} s  "
                  f"p_max={p.max()/pscale:.3f} {a.unit}  ({kind}: {src})")

    if not curves:
        sys.exit("no curves could be extracted")

    # ---- CSV export ------------------------------------------------------- #
    if a.csv:
        with open(a.csv, "w") as fh:
            fh.write(f"project,node,time_s,{a.var}_Pa\n")
            for name, (t, p) in curves.items():
                for ti, pi in zip(t, p):
                    fh.write(f"{name},{a.node},{ti:.6g},{pi:.10g}\n")
        print(f"  wrote {a.csv}")

    # ---- plot ------------------------------------------------------------- #
    fig, ax = plt.subplots(figsize=(10, 5.5))
    if meas is not None:
        tm, pm, mname = meas
        ax.plot(tm / tscale, pm / pscale, color="0.35", lw=1.0, zorder=1,
                label=f"Messung ({mname}" + (", Baseline abgezogen)" if a.baseline_end >= 0 else ")"))
    for name, (t, p) in curves.items():
        y = (p - p[0]) if relative else p
        label = name
        if meas is not None:
            e = rmse(t, y, meas[0], meas[1])
            if e is not None:
                label += f"  (RMSE {e/pscale:.3g} {a.unit})"
                print(f"  [rmse] {name:28s} {e/1e6:.4f} MPa")
            else:
                print(f"  [rmse] {name:28s} no time overlap with measurement")
        ax.plot(t / tscale, y / pscale, lw=1.6, zorder=2, label=label)
    where = INJECTION_NODES.get(a.node, f"node {a.node}")
    ax.set_title(f"STIMTEC DFN – Druckverlauf am Injektionspunkt BH10 / {where}")
    ax.set_xlabel(f"Zeit [{a.time_unit}]")
    ax.set_ylabel(("Druckaufbau Δp" if relative else a.var) + f" [{a.unit}]")
    if a.ylim:
        ax.set_ylim(*a.ylim)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc="best")
    fig.tight_layout()
    if a.save:
        fig.savefig(a.save, dpi=150)
        print(f"  wrote {a.save}")
    if not a.no_show:
        plt.show()


if __name__ == "__main__":
    main()
