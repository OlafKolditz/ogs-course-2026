import pyvista as pv

# Load the circle mesh
mesh = pv.read('circle_without_central_vertex_cell.vtu')

# Extract outer boundary edges
boundary_edges = mesh.extract_feature_edges(
    boundary_edges=True,
    feature_edges=False,
    manifold_edges=False,
    non_manifold_edges=False
)

# Save as VTP (simplest)
boundary_edges.save('circle_boundary.vtp')

# OR save as VTU (after conversion)
boundary_edges.cast_to_unstructured_grid().save('circle_boundary.vtu')

# Visualize (optional)
plotter = pv.Plotter()
plotter.add_mesh(mesh, style='wireframe', color='lightblue', opacity=0.3)
plotter.add_mesh(boundary_edges, color='red', line_width=3)
plotter.show()