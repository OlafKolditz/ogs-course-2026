import subprocess
from pathlib import Path

import ogstools as ot
import matplotlib.pyplot as plt
import numpy as np

#-##########################################################
# 1 Plotting framework
plt.rcParams.update({
    'axes.titlesize': 14,
    'axes.labelsize': 14,
    'xtick.labelsize': 14,
    'ytick.labelsize': 14,
    'legend.fontsize': 14,
})
fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(20, 10))

BASE_DIR = Path(__file__).resolve().parent

#-##########################################################
# 2 Load numerical results: basic vs. leakage
# NOTE: this compares the two numerical solutions (basic, leakage) against
# each other. The analytical Theis solution is intentionally NOT plotted.
xn = 100
time_vals = [8.64, 86.4, 1728.0, 24192.0, 172800.0, 604800.0, 864000.0]

# Each case writes its results next to its project file (see <prefix> in
# the .prj files). If results are missing, OGS is run in that folder first.
CASES = {
    "basic": (BASE_DIR / "basic", "ogs-basic.prj", "theis2Dcircle_basic.pvd"),
    "leakage": (BASE_DIR / "leakage", "ogs-leakage.prj", "theis2Dcircle_leakage.pvd"),
}


def load_case(name, rerun=False):
    folder, prj, pvd = CASES[name]
    if rerun or not (folder / pvd).exists():
        subprocess.run(["ogs", prj, "-o", str(folder)], cwd=folder, check=True)
    return ot.MeshSeries(str(folder / pvd))


ms_basic = load_case("basic")
ms_leakage = load_case("leakage")
pressure = ot.variables.pressure.replace(data_unit="m", output_unit="m", output_name="hydraulic head")

#-#############################
# 2.1 Profiles along a line at selected times: basic (solid) vs leakage (dashed)
xaxis = np.column_stack((np.linspace(0.0, xn, xn), np.zeros((xn, 2))))

probe_basic = ms_basic.probe(xaxis)
labels = [f"$t={np.round(x, 2)}s$" for x in probe_basic[1:].timevalues]
ot.plot.line(probe_basic[1:], "x", pressure, labels=labels, ax=ax[0],
             linestyle="-", fontsize=14)

# overlay leakage case with the same color cycle, dashed, unlabeled
# (matching color = matching time; solid vs dashed = basic vs leakage)
probe_leakage = ms_leakage.probe(xaxis)
ot.plot.line(probe_leakage[1:], "x", pressure, ax=ax[0],
             linestyle="--", fontsize=14)

ax[0].set_title("Profiles at different times (basic solid, leakage dashed)")
ax[0].set_ylabel("hydraulic head / drawdown [m]")

#-#############################
# 2.2 Temporal evolution at selected points: basic (solid) vs leakage (dashed)
x_vals = [0.3048, 1, 10, 20]
points_observation = np.array([[x, 0.0, 0.0] for x in x_vals])

point_series_basic = ms_basic.probe(points_observation)
labels = [f"$x={x}m$" for x in x_vals]
ot.plot.line(point_series_basic, "time", pressure, labels=labels, ax=ax[1],
             linestyle="-", fontsize=14)

point_series_leakage = ms_leakage.probe(points_observation)
ot.plot.line(point_series_leakage, "time", pressure, ax=ax[1],
             linestyle="--", fontsize=14)

ax[1].set_title("Temporal evolution (basic solid, leakage dashed)")
ax[1].set_xscale("log")
ax[1].set_ylabel("hydraulic head / drawdown [m]")

plt.tight_layout()
plt.savefig(BASE_DIR / "basic_vs_leakage_comparison.png", dpi=150)
plt.show()
