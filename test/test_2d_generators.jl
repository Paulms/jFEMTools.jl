include("../src/mesh.jl")
include("../src/mesh_generator.jl")

import Tensors
mesh = rectangle_mesh(TriangleCell, (2,2), Tensors.Vec{2}((0.0,0.0)), Tensors.Vec{2}((1.0,1.0)))

# Plots
#plot solution
coordinates = get_vertices_matrix(mesh);
connectivity = get_conectivity_list(mesh);
using Makie
scene = Scene()
for row in 1:getncells(mesh)
	#read coordinates
	points = node(:poly, Point2f0[coordinates[node,:] for node in connectivity[row]])
	poly!(scene, points, strokewidth = 1, color = :white, strokecolor = :black, show_axis = true, scale_plot = false)
end
display(scene)

# Quad mesh
mesh = rectangle_mesh(RectangleCell, (2,2), Tensors.Vec{2}((0.0,0.0)), Tensors.Vec{2}((1.0,1.0)))

# Mixed mesh
cells = [
	TriangleCell((1,2,4)),
	RectangleCell((2,3,6,5)),
	TriangleCell((5,4,2)),
    RectangleCell((4, 5, 8, 7)),
	RectangleCell((5,6,9,8))
    ];
mixmesh = PolytopeMesh(cells, mesh.nodes)

coordinates = get_vertices_matrix(mixmesh);
connectivity = get_conectivity_list(mixmesh);
scene = Scene()
for row in 1:getncells(mixmesh)
	#read coordinates
	points = node(:poly, Point2f0[coordinates[node,:] for node in connectivity[row]])
	poly!(scene, points, strokewidth = 1, color = :white, strokecolor = :black, show_axis = true, scale_plot = false)
end
display(scene)
