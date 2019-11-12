using jFEMTools

import Tensors
using SparseArrays
using LinearAlgebra

mesh = unitSquareMesh(RectangleCell, (2,2));

degree = 2;
dim = 2;
element = VirtualElement(dim,degree);
dofs = DofHandler(mesh, element);
operators = VEMOperators(dofs);
K = get_K(operators);

b = get_constant_load(operators, 1);
boundaries = Boundaries(operators, b);

# Boundary condition
g(x::Tensors.Vec{2}) = (1 - x[1])*x[2]*sin(π * x[1]);
dbc = Dirichlet(boundaries, "boundary", g)
rhs = get_rhs(boundaries)

#Solve
x = K\rhs ;
u = get_nodal_values(x, dofs);
#Plot solution
using Makie
#popdisplay(AbstractPlotting.PlotDisplay())
#AbstractPlotting.inline!(true)
scene = Scene(resolution = (400, 200))
plot!(scene, mesh)
#display(scene)
