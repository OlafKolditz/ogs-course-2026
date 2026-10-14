"""Erzeugt stimtec_dfn_st.prj aus stimtec_dfn.prj:
Injektion als nodale Quellterme an den Schnittpunkten BH10 x Kluft 1/2 (Knoten 270, 271)
statt Neumann-RB auf BH10.vtu; bc_zmax (3D-Modell) entfernt."""
import re, sys
args = [a for a in sys.argv[1:] if not a.startswith('--')]
src = args[0] if len(args) > 0 else 'stimtec_dfn.prj'
dst = args[1] if len(args) > 1 else 'stimtec_dfn_st.prj'
RIM = '--rim' in sys.argv   # optional: Dirichlet p0 auf den Kluftraendern (bc_rim.vtu)
L_BH10 = 55.6          # q_in ist auf die BH10-Laenge bezogen (kg/s/m, s. Kommentar /55.6)
SPLIT1 = 0.5           # Anteil der Gesamtrate in Kluft 1 (Rest -> Kluft 2)
s = open(src, encoding='ISO-8859-1').read()

# 1) Meshes
s = re.sub(r'<meshes>.*?</meshes>', '''<meshes>
        <mesh>mesh.vtu</mesh>
        <!-- Schnittpunkte BH10 x Kluftnetz (bulk_node_ids 270, 271) -->
        <mesh>BH10_frac1.vtu</mesh>
        <mesh>BH10_frac2.vtu</mesh>
        <mesh>BH10_intersections.vtu</mesh>
    </meshes>''', s, count=1, flags=re.S)

# 2) Output auch an den Schnittpunkten
s = re.sub(r'<!-- <meshes>\s*<mesh>mesh</mesh>\s*<mesh>BH10</mesh>\s*</meshes> -->',
           '<meshes>\n                <mesh>mesh</mesh>\n                <mesh>BH10_intersections</mesh>\n            </meshes>', s)

s = s.replace('<prefix>Stimtec_DFN</prefix>', '<prefix>Stimtec_DFN_{:meshname}</prefix>')

# 3) Quellterm-Parameter aus q_in-Expression ableiten
m = re.search(r'(<parameter>\s*<name>q_in</name>\s*<type>Function</type>\s*<expression>)(.*?)(<!--.*?-->)?\s*(</expression>\s*</parameter>)', s, re.S)
expr = ' '.join(m.group(2).split())
def par(name, f, comment):
    return f'''        <parameter>
            <!-- {comment} -->
            <name>{name}</name>
            <type>Function</type>
            <expression>{f:g}*{L_BH10:g}*{expr}</expression>
        </parameter>
'''
newpars = par('q_frac1', SPLIT1, f'Massenrate Kluft 1 [kg/s] = {SPLIT1:g} * L_BH10 * q_in') + \
          par('q_frac2', 1 - SPLIT1, f'Massenrate Kluft 2 [kg/s] = {1-SPLIT1:g} * L_BH10 * q_in')
end = m.end()
s = s[:end] + '\n' + newpars + s[end:]

# 4) Randbedingungen: bc_zmax und Neumann BH10 auskommentieren, Quellterme hinzufuegen
s = re.sub(r'(<boundary_condition>\s*<type>Dirichlet</type>\s*<mesh>bc_zmax</mesh>.*?</boundary_condition>)',
           lambda k: '<!-- bc_zmax gehoert zum 3D-Modell\n                ' + k.group(1) + ' -->', s, flags=re.S)
s = re.sub(r'(<boundary_condition>\s*<type>Neumann</type>\s*<mesh>BH10</mesh>.*?</boundary_condition>)',
           lambda k: '<!-- ersetzt durch nodale Quellterme\n                ' + k.group(1) + ' -->', s, flags=re.S)
s = s.replace('</boundary_conditions>', '''</boundary_conditions>
            <source_terms>
                <source_term>
                    <mesh>BH10_frac1</mesh>
                    <type>Nodal</type>
                    <parameter>q_frac1</parameter>
                </source_term>
                <source_term>
                    <mesh>BH10_frac2</mesh>
                    <type>Nodal</type>
                    <parameter>q_frac2</parameter>
                </source_term>
            </source_terms>''', 1)
if RIM:
    s = s.replace('<mesh>BH10_intersections.vtu</mesh>', '<mesh>BH10_intersections.vtu</mesh>\n        <!-- Kluftraender (262 Knoten), Dirichlet p0 -->\n        <mesh>bc_rim.vtu</mesh>', 1)
    s = s.replace('</boundary_conditions>', '''    <boundary_condition>
                    <type>Dirichlet</type>
                    <mesh>bc_rim</mesh>
                    <parameter>p0</parameter>
                </boundary_condition>
            </boundary_conditions>''', 1)
open(dst, 'w', encoding='ISO-8859-1').write(s)
print('geschrieben:', dst)
