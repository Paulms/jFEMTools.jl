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
