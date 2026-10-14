Theis

14.10.2026: theis output in actual directory

(OK1)
Please modify the project file axisym_theis.prj in folder theis1Dradial that the output results are written in the actual folder (i.e. where the project file is).

(OK2)
Please modify the Python code axisym_theis.py that the output results are written in the actual folder (i.e. where the project file is) and the comparison with the analytical solution works.

(OK3)
Please ckeck the theis numerical models in folder 94-CLD/theis. The basic case and the leakage case should write results in the actual folder (where the project file is) and function plot_compare_numerical_solutions.py (in the folder above) should show the comparision of basic and leakage cases.

09.10.2026

Python-Version for leakage

(OK1)
After writing the essay, now we go back to the OGS simulation. The leakage term in ogs-leakage.prj is expressed as <source_term> connected with a <mesh>mesh_leakage</mesh>. Please use a Python-Funktion instead and use the OGS Python interface for modeling the leakage. The OGS project you find in subfolder 94-CLD/theis/leakage.
