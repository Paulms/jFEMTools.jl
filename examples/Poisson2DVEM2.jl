using jFEMTools

import Tensors
using SparseArrays
using LinearAlgebra

mesh = unitSquareMesh(RectangleCell, (1,1));

degree = 1;
dim = 2;
element = VirtualElement(dim,degree);
dofs = DofHandler(mesh, element);
operators = VEMOperators(dofs, element);
# Test
d = sqrt(2)
Bref = 1/4*[1 1 1 1;
            -d d -d d;
            -d -d d d]
Bref ≈ operators.elements[1].B

Dref = 1/4*[4 -d -d;
            4 d -d;
            4 d  d;
            4 -d d]

Dref ≈ operators.elements[1].D

Gref = 1/2*[2 0 0;
            0 1 0;
            0 0 1]

Gref ≈ operators.elements[1].G


degree = 2;
element = VirtualElement(dim,degree);
dofs = DofHandler(mesh, element);
operators = VEMOperators(dofs, element);
# Test
d = sqrt(2)
Bref = 1/12*[0 0 0 0 0 0 0  0 12;
             -d d d -d 0 4*d 0 -4*d 0;
             -d -d d d -4*d 0 4*d 0 0;
             1 1 1 1 0 4 0 4 -12;
             1 -1 1 -1 0 0 0 0 0;
             1 1 1 1 4 0 4 0 -12]
Bref ≈ operators.elements[1].B

Dref = 1/24*[24 -6*d -6*d 3 3 3;
            24 6*d -6*d 3 -3 3;
            24 6*d 6*d 3 3 3;
            24 -6*d 6*d 3 -3 3;
            24 0 -6*d 0 0 3;
            24 6*d 0 3 0 0;
            24 0 6*d 0 0 3;
            24 -6*d 0 3 0 0;
            24 0 0 1 0 1]

Dref ≈ operators.elements[1].D

Gref = 1/24*[24 0 0 1 0 1;
            0 12 0 0 0 0;
            0 0 12 0 0 0;
            0 0 0 2 0 0;
            0 0 0 0 1 0;
            0 0 0 0 0 2]

Gref ≈ operators.elements[1].G


K = get_K(operators)
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
