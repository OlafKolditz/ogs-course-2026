"""
Leakage source term for the Theis problem (OGS LiquidFlow, 2D) via the
OpenGeoSys Python interface.

Referenced in ogs-leakage-python.prj by
    <python_script>leakage_source_term.py</python_script>
    <source_term>
        <mesh>mesh_leakage</mesh>
        <type>Python</type>
        <source_term_object>leakage</source_term_object>
    </source_term>

OGS calls getFlux(t, coords, primary_vars) at every integration point of the
source-term mesh and integrates  int_Omega N^T * value dOmega  into the RHS,
i.e. the Python ST behaves like the former <type>Volumetric</type> ST.
It must return (value, d value / d primary_vars).

Select the leakage model with MODE:
  "constant"  q_L = Q_CONST                          (= old prj, -1e-7)
  "time"      q_L = A_TIME / sqrt(t)                 (= time-dependent prj)
  "hantush"   q_L = -LEAKANCE * (p - P_REF)          (Hantush-Jacob leaky aquifer,
                                                      flux through aquitard ~ head difference)

Note: "constant" and "time" run with the Picard setup (ogs-leakage-python.prj).
"hantush" depends on p; with Picard the term is lagged and diverges for large
time steps -> use ogs-leakage-python-newton.prj (Newton, which uses dq/dp).
"""

import math

try:  # OGS >= 6.5 (wheel / embedded interpreter)
    import ogs.callbacks as OpenGeoSys
except ModuleNotFoundError:  # older OGS versions
    import OpenGeoSys

# ---------------------------------------------------------------- parameters
MODE = "constant"

Q_CONST = -1.0e-7  # constant leakage rate                [1/s]
A_TIME = -1.0e-5  # amplitude of the 1/sqrt(t) leakage    [1/s^0.5]
LEAKANCE = 1.0e-7  # aquitard leakance K'/b' (per unit p)    [1/(Pa s)]
P_REF = 0.0  # pressure (head) above the aquitard       [Pa]


class LeakageSourceTerm(OpenGeoSys.SourceTerm):
    def getFlux(self, t, coords, primary_vars):
        # x, y, z = coords   # available for spatially varying leakage
        p = primary_vars[0]
        dq_dp = 0.0

        if MODE == "constant":
            q = Q_CONST
        elif MODE == "time":
            q = A_TIME / math.sqrt(t) if t > 0.0 else 0.0
        elif MODE == "hantush":
            q = -LEAKANCE * (p - P_REF)
            dq_dp = -LEAKANCE
        else:
            raise ValueError(f"unknown leakage MODE '{MODE}'")

        # derivative w.r.t. every primary variable (only pressure here);
        # used by OGS for the Jacobian (Newton), ignored by Picard
        return (q, [dq_dp])


# object referenced in the prj file via <source_term_object>
leakage = LeakageSourceTerm()
