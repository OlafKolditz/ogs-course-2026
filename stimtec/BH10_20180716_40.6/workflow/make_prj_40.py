"""
Erzeugt die Projektdateien fuer 40.6 aus den Vorlagen von 51.6
(../../BH10_20180716_51.6/simulation/*.prj): gleiche Modelle, aber
 - Zeitraum 1750..2670 s, dt 5 s
 - pressure_source_term aus injektion_stufen_40.csv
 - permeability_t / storage_t: gleiche Stufenstruktur (Startwerte = Optimum 51.6 Variante 4)
ACHTUNG: ueberschreibt die kalibrierten prj-Dateien in ../simulation (danach Kalibrierung neu starten)
"""
import os, re, sys
import bh40

HERE = bh40.HERE
SRC = os.path.normpath(os.path.join(HERE, "..", "..", "BH10_20180716_51.6", "simulation"))
if len(sys.argv) > 1:
    SRC = sys.argv[1]
DST = os.path.normpath(os.path.join(HERE, "..", "simulation"))
FILES = ["stimtec.prj", "stimtec_bc_right.prj", "stimtec_bc_right_kS_t.prj",
         "stimtec_bc_rightmid_kS_t.prj", "stimtec_bc_rightmid_py.prj"]


def set_expr(c, name, values):
    pat = re.compile(rf"(<name>{name}</name>\s*<type>Function</type>\s*<expression>)(.*?)(</expression>)", re.S)
    c, n = pat.subn(lambda m: m.group(1) + bh40.expression(values) + m.group(3), c)
    assert n == 1, name
    return c


for fn in FILES:
    with open(os.path.join(SRC, fn), encoding="ISO-8859-1") as f:
        c = f.read()
    c = re.sub(r"<t_initial>[^<]*</t_initial>", f"<t_initial>{bh40.T_START}</t_initial>", c)
    c = re.sub(r"<t_end>[^<]*</t_end>", f"<t_end>{bh40.T_END}</t_end>", c)
    c = re.sub(r"<delta_t>[^<]*</delta_t>", f"<delta_t>{bh40.DT:g}</delta_t>", c)
    c = set_expr(c, "pressure_source_term", bh40.Q_STEPS)
    if "permeability_t" in c and "<name>permeability_t</name>" in c:
        c = set_expr(c, "permeability_t", bh40.mapped_values(3.3679e-15, 2.4110e-14))
        c = set_expr(c, "storage_t", bh40.mapped_values(2.4790e-14, 1.4138e-11))
        c = re.sub(r"<!-- k\(t\).*?-->", "<!-- k(t): gleiche Zeitfenster wie pressure_source_term (40.6) -->", c)
        c = re.sub(r"<!-- S\(t\).*?-->", "<!-- S(t): gleiche Zeitfenster wie pressure_source_term (40.6) -->", c)
    with open(os.path.join(DST, fn), "w", encoding="ISO-8859-1") as f:
        f.write(c)
    print("geschrieben:", fn)
