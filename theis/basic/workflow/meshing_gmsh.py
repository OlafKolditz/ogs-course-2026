import gmsh
import math
import pyvista as pv

gmsh.initialize()
gmsh.model.add("circle_with_graded_mesh")

R = 304.8

# Add points
center = gmsh.model.geo.addPoint(0.0, 0.0, 0.0)
points = []
n = 8
for i in range(n):
    theta = 2*math.pi * i/n
    x = R * math.cos(theta)
    y = R * math.sin(theta)
    points.append(gmsh.model.geo.addPoint(x, y, 0))

# Create geometry
arcs = []
for i in range(n-1):
    arcs.append(gmsh.model.geo.addCircleArc(points[i], center, points[i+1]))
arcs.append(gmsh.model.geo.addCircleArc(points[n-1], center, points[0]))

curve_loop = gmsh.model.geo.addCurveLoop(arcs)
surface = gmsh.model.geo.addPlaneSurface([curve_loop])
gmsh.model.geo.synchronize()

gmsh.model.mesh.embed(0, [center], 2, surface)

# ====== MESH SIZE CONTROL USING FIELDS ======
# Create a Distance field from the center point
tag_dist = gmsh.model.mesh.field.add("Distance")
gmsh.model.mesh.field.setNumbers(tag_dist, "PointsList", [center])

# Create a Threshold field for size transition
tag_threshold = gmsh.model.mesh.field.add("Threshold")
gmsh.model.mesh.field.setNumber(tag_threshold, "InField", tag_dist)
gmsh.model.mesh.field.setNumber(tag_threshold, "SizeMin", 1.0)    # Size at center
gmsh.model.mesh.field.setNumber(tag_threshold, "SizeMax", 20.)     # Size far from center
gmsh.model.mesh.field.setNumber(tag_threshold, "DistMin", 0.0)     # Start transition at center
gmsh.model.mesh.field.setNumber(tag_threshold, "DistMax", 250.)     # End transition at radius 0.5

# Set this as the background mesh field
gmsh.model.mesh.field.setAsBackgroundMesh(tag_threshold)

# Generate mesh
gmsh.model.mesh.generate(2)
gmsh.write("circle_graded_mesh.msh")
gmsh.write("circle_graded_mesh.vtk")

elementTypes, elementTags, nodeTags = gmsh.model.mesh.getElements()
print("Element types present after generation:", elementTypes)


gmsh.fltk.run()
gmsh.finalize()

# Creating 
mesh = pv.read('circle_graded_mesh.vtk')
mesh.save('circle_graded_mesh.vtu')

# See the names of all data arrays in your file
print("Available point data arrays:", mesh.point_data.keys())
print("Available cell data arrays:", mesh.cell_data.keys())

# Create a plotter and add the mesh
plotter = pv.Plotter() # Create a plotting object[citation:2]
plotter.add_mesh(
    mesh,
    show_edges=True,  # Makes cell edges visible for a wireframe-like view[citation:1]
    color='lightblue'
)
# CRITICAL: Set the camera to view the X-Y plane (top-down view)
plotter.view_xy() # Looks along the -Z axis[citation:3]
# Show the plot
plotter.show()
##plotter.show(interactive=False, screenshot='square_with_center_node.png')