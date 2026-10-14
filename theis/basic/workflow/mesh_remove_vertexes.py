import gmsh
import sys
import pyvista as pv
import os
import ogstools as ot
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
###########################################################
# Correcting mesh
# Remove vertex elemens
# Load the mesh with the central vertex
mesh = pv.read('circle_graded_mesh.vtu')
# Get mesh information
print(f"Number of points: {mesh.n_points}")
print(f"Number of cells: {mesh.n_cells}")
# This removes the central point cell but keeps the surface/volume elements
non_vertex_mask = mesh.celltypes != 1
# Apply the mask to extract the desired cells
mesh_without_vertex_cells = mesh.extract_cells(non_vertex_mask)
# Save the clean mesh
mesh_without_vertex_cells.save('circle_without_central_vertex_cell.vtu')
print(f"Removed vertex cells. New cell count: {mesh_without_vertex_cells.n_cells}")

# Remove line elements
mesh = pv.read('circle_without_central_vertex_cell.vtu')
print(f"Number of points: {mesh.n_points}")
print(f"Number of cells: {mesh.n_cells}")
# This removes the central point cell but keeps the surface/volume elements
non_line_mask = mesh.celltypes != 3
# Apply the mask to extract the desired cells
mesh_without_line_cells = mesh.extract_cells(non_line_mask)
# Save the clean mesh
mesh_without_line_cells.save('circle_without_line_cells.vtu')
print(f"Removed line cells. New cell count: {mesh_without_line_cells.n_cells}")