#!/usr/bin/env python3
"""
plot_mesh.py

Stellt das Mesh einer VTU-Datei (VTK XML UnstructuredGrid, z.B. aus OGS/OpenGeoSys)
grafisch dar - Knoten, Kanten/Elemente und optional Knoten-IDs.

Gedacht u.a. fuer Randmesh-Dateien wie "outer_boundary.vtu" (Linienelemente, die den
Rand eines 2D-Gebiets beschreiben), funktioniert aber generisch auch fuer Dreiecks-,
Viereck- und gemischte 2D-Meshes.

Verwendung:
    python plot_mesh.py                          # nutzt outer_boundary.vtu im aktuellen Verzeichnis
    python plot_mesh.py pfad/zur/datei.vtu
    python plot_mesh.py datei.vtu --node-ids      # Knoten-IDs anzeigen
    python plot_mesh.py datei.vtu --save mesh.png # als PNG speichern statt anzuzeigen
    python plot_mesh.py datei.vtu --no-points     # nur Kanten, keine Knotenmarker

Abhaengigkeiten:
    pip install meshio matplotlib numpy
"""

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

try:
    import meshio
except ImportError:
    sys.exit(
        "Das Paket 'meshio' wird benoetigt. Installation mit:\n"
        "    pip install meshio"
    )

# VTK-Zelltypen, die als Linien (Kanten) dargestellt werden
LINE_TYPES = {"line", "line3"}
# VTK-Zelltypen, die als gefuellte/umrandete Flaechen dargestellt werden
POLY_TYPES = {"triangle", "triangle6", "quad", "quad8", "quad9", "polygon"}


def parse_args():
    p = argparse.ArgumentParser(description="Visualisierung eines VTU-Meshes")
    p.add_argument(
        "vtu_file",
        nargs="?",
        default="outer_boundary.vtu",
        help="Pfad zur .vtu-Datei (Default: outer_boundary.vtu)",
    )
    p.add_argument(
        "--node-ids", action="store_true", help="Knoten-IDs im Plot anzeigen"
    )
    p.add_argument(
        "--no-points", action="store_true", help="Keine Knotenmarker zeichnen"
    )
    p.add_argument(
        "--save",
        metavar="DATEI",
        help="Plot als Bilddatei speichern (z.B. mesh.png) statt interaktiv anzuzeigen",
    )
    p.add_argument(
        "--dpi", type=int, default=150, help="Aufloesung beim Speichern (Default: 150)"
    )
    return p.parse_args()


def load_mesh(path: Path) -> meshio.Mesh:
    if not path.exists():
        sys.exit(f"Datei nicht gefunden: {path}")
    return meshio.read(path)


def plot_mesh(mesh: meshio.Mesh, title: str, show_node_ids: bool, show_points: bool):
    points = mesh.points  # (N, 3)
    x, y = points[:, 0], points[:, 1]

    fig, ax = plt.subplots(figsize=(8, 8))

    n_line_cells = 0
    n_poly_cells = 0

    for cell_block in mesh.cells:
        ctype = cell_block.type
        conn = cell_block.data

        if ctype in LINE_TYPES:
            # jede Zeile in conn ist ein Linienelement [i, j]
            for i, j in conn[:, :2]:
                ax.plot(
                    [x[i], x[j]], [y[i], y[j]],
                    color="tab:blue", linewidth=1.2, zorder=2,
                )
            n_line_cells += len(conn)

        elif ctype in POLY_TYPES:
            for elem in conn:
                poly_x = np.append(x[elem], x[elem[0]])
                poly_y = np.append(y[elem], y[elem[0]])
                ax.fill(poly_x, poly_y, facecolor="tab:blue", alpha=0.08, zorder=1)
                ax.plot(poly_x, poly_y, color="tab:blue", linewidth=0.6, zorder=2)
            n_poly_cells += len(conn)

        else:
            print(f"Hinweis: Zelltyp '{ctype}' wird nicht speziell dargestellt "
                  f"und wird uebersprungen.")

    if show_points:
        ax.scatter(x, y, s=14, color="tab:red", zorder=3, label="Knoten")

    if show_node_ids:
        for idx, (xi, yi) in enumerate(zip(x, y)):
            ax.annotate(str(idx), (xi, yi), fontsize=6, color="dimgray",
                        xytext=(2, 2), textcoords="offset points")

    ax.set_aspect("equal", adjustable="datalim")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    n_elem = n_line_cells + n_poly_cells
    ax.set_title(f"{title}\n{len(points)} Knoten, {n_elem} Elemente")
    if show_points:
        ax.legend(loc="best", fontsize=8)
    ax.grid(True, linestyle=":", linewidth=0.5, alpha=0.6)
    fig.tight_layout()
    return fig


def main():
    args = parse_args()
    path = Path(args.vtu_file)
    mesh = load_mesh(path)

    print(f"Datei:          {path}")
    print(f"Anzahl Knoten:  {len(mesh.points)}")
    for cb in mesh.cells:
        print(f"Zelltyp '{cb.type}': {len(cb.data)} Elemente")

    fig = plot_mesh(
        mesh,
        title=path.name,
        show_node_ids=args.node_ids,
        show_points=not args.no_points,
    )

    if args.save:
        fig.savefig(args.save, dpi=args.dpi)
        print(f"Plot gespeichert unter: {args.save}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
