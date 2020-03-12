# VEM computes the virtual element solution of a Poisson problem on a polygonal mesh
# Translation to Julia of MATLAB code from paper:
# O. Suttton, The virtual element method in 50 lines of MATLAB, Numer Algor(2017) 75:1141-1159

using jFEMTools

import Tensors
using SparseArrays
using LinearAlgebra

mesh = rectangle_mesh(RectangleCell, (2,2), Tensors.Vec{2}((0.0,0.0)), Tensors.Vec{2}((1.0,1.0)));

function compute_u(mesh)
	# forcing function
	rhs(x::Tensors.Vec{2}) = 2*π^2*sin(π*x[1])*sin(π*x[2])
	# Boundary condition
	g(x::Tensors.Vec{2}) = sin(π * x[2])*sin(π * x[1]);
	n_dofs = getnvertices(mesh)
	n_polys = 3; # Method has 1 degree of freedom per vertex
	assembler = start_assemble(n_dofs)
	F = zeros(n_dofs); # Forcing vector
	u = zeros(n_dofs); # Degrees of freedom of the virtual element solution
	linear_polynomials = ((0,0), (1,0), (0,1)); # Impose an ordering on the linear polynomials
	for el_id = 1:getncells(mesh)
		vert_ids = getverticesidx(mesh, el_id) # Global IDs of the vertices of this element
		verts = getverticescoords(mesh, el_id) # Coordinates of the vertices of this element
		n_sides = getnedges(mesh.cells[el_id]) # Start computing the geometric information
		area = cell_volume(mesh, el_id)
		centr = cell_centroid(mesh, el_id)
		diameter = cell_diameter(mesh, el_id)
		D = zeros(n_sides, n_polys); D[:, 1] .= 1;
		B = zeros(n_polys, n_sides); B[1, :] .= 1/n_sides;
		for vertex_id = 1:n_sides
			vert = verts[vertex_id]; #This vertex and its neighbours
			prev = verts[mod1(vertex_id - 1, n_sides)];
			nextv = verts[mod1(vertex_id + 1, n_sides)];
			vertex_normal = [nextv[2] - prev[2], prev[1] - nextv[1]]; # Average of edge normals
			for poly_id in 2:n_polys # Only need to loop over non-constant polynomials
				poly_degree = linear_polynomials[poly_id];
				monomial_grad = poly_degree ./ diameter; # Gradient of a linear is constant
				D[vertex_id, poly_id] = dot(vert - centr, poly_degree) / diameter;
				B[poly_id, vertex_id] = 0.5 * dot(monomial_grad, vertex_normal);
			end
		end
		projector = (B*D) \ B; # Compute the local Ritz projector to polynomials
		stabilising_term = (I - D * projector)' * (I - D * projector);
		G = B*D; G[1, :] .= 0;
		local_stiffness = projector' * G * projector + stabilising_term;
		assemble!(assembler,vert_ids,local_stiffness)
		assemble!(F, [dof for dof in vert_ids], rhs(centr)*ones(length(vert_ids)) * area / n_sides)
	end
	K = end_assemble(assembler)
	boundary_nodes = getvertexset(mesh,"boundary")
	boundary_vals = Vector{Float64}(undef, length(boundary_nodes))
	for (i,node) in enumerate(boundary_nodes)
		boundary_vals[i] = g(getvertexcoords(mesh, node))
	end

	internal_dofs = setdiff(1:n_dofs, boundary_nodes)
	F = F - K[:, collect(boundary_nodes)] * boundary_vals; # Apply the boundary condition
	u[internal_dofs] = K[collect(internal_dofs), collect(internal_dofs)] \ F[collect(internal_dofs)]; # Solve
	u[collect(boundary_nodes)] = boundary_vals; # Set the boundary values
	u
end
u = compute_u(mesh)

#plot solution
using Makie
import AbstractPlotting
using ColorSchemes
include("../src/plot_recipes.jl")
scene = Scene(resolution = (400, 200), colormap = ColorSchemes.viridis.colors)
plot!(scene, mesh, color = u)
display(scene)
