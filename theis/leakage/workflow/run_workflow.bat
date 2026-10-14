python meshing_gmsh.py
python mesh_remove_vertexes.py
ExtractBoundary -i circle_without_central_vertex_cell.vtu -o outer_boundary.vtu
rem identifySubdomains -m mesh_domain.vtu -- mesh_leakage.vtu
rem ogs theis-2a-run.prj
rem python plot_numerical_vs_analytical.py
