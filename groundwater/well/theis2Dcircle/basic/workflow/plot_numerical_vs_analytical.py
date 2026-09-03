import ogstools as ot
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from scipy.special import exp1

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

# Aquifer properties (must match theis-2a-run.prj)
S = 0.001
T = 9.2903e-4
Q = 0.016  # Pumping rate from well (m3/s)

#-##########################################################
# 2 Analytical Theis solution
def calc_u(r, S, T, t):
    """Calculate and return the dimensionless time parameter, u."""
    return r**2 * S / 4 / T / t

def theis_drawdown(t, S, T, Q, r):
    """Calculate and return the drawdown s(r,t) for parameters S, T.

    Theis equation, s(r,t) = Q * W(u) / (4.pi.T),
    where W(u) is the Well function for u = Sr^2 / (4Tt).
    """
    u = calc_u(r, S, T, t)
    return Q / 4 / np.pi / T * exp1(u)

#-##########################################################
# 3 Load numerical results
xn = 100
rn = 100
r = np.arange(1, rn + 1, 1)
time_vals = [8.64, 86.4, 1728.0, 24192.0, 172800.0, 604800.0, 864000.0]

ms = ot.MeshSeries("../results/theis2Dcircle.pvd")
pressure = ot.variables.pressure.replace(data_unit="m", output_unit="m", output_name="hydraulic head")

#-#############################
# 3.1 Profiles along a line at selected times: numerical vs analytical
xaxis = np.column_stack((np.linspace(0.0, xn, xn), np.zeros((xn, 2))))
ms_probe = ms.probe(xaxis)

labels = [f"$t={np.round(x, 2)}s$" for x in ms_probe[1:].timevalues]
ot.plot.line(ms_probe[1:], "x", pressure, labels=labels, ax=ax[0], fontsize=14)

# overlay analytical Theis drawdown for the same times, reusing each numerical line's own color
numerical_colors = [line.get_color() for line in ax[0].get_lines()]
for color, t_val in zip(numerical_colors, time_vals):
    s_theis = theis_drawdown(t_val, S, T, Q, r)
    ax[0].plot(r, s_theis, "--", color=color, linewidth=1.5)

ax[0].set_title("Profiles at different times (numerical solid, Theis dashed)")
ax[0].set_ylabel("hydraulic head / drawdown [m]")

#-#############################
# 3.2 Temporal evolution at selected points: numerical vs analytical
x_vals = [0.3048, 1, 10, 20]
points_observation = np.array([[x, 0.0, 0.0] for x in x_vals])
point_series = ms.probe(points_observation)

labels = [f"$x={x}m$" for x in x_vals]
ot.plot.line(point_series, "time", pressure, labels=labels, ax=ax[1], fontsize=14)

# overlay analytical Theis drawdown for the same observation points
t_continuous = np.logspace(np.log10(time_vals[0]), np.log10(time_vals[-1]), 200)
numerical_colors_t = [line.get_color() for line in ax[1].get_lines()]
for color, x_obs in zip(numerical_colors_t, x_vals):
    s_theis_t = theis_drawdown(t_continuous, S, T, Q, x_obs)
    ax[1].plot(t_continuous, s_theis_t, "--", color=color, linewidth=1.5)

ax[1].set_title("Temporal evolution (numerical solid, Theis dashed)")
ax[1].set_xscale("log")
ax[1].set_ylabel("hydraulic head / drawdown [m]")

plt.tight_layout()
plt.show()
