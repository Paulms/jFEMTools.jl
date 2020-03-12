using jFEMTools
import Tensors

mesh = unitSquareMesh2(RectangleCell, (3,3));

# forcing function
rhs(x::Tensors.Vec{2}) = 2*π^2*sin(π*x[1])*sin(π*x[2])

degree = 2;
dim = 2;
element = PoissonVirtualElement(dim,degree);
u = TrialFunction(element)
dofs = DofHandler(mesh, u);
operators = VEMOperators(dofs, u;load = rhs);

K = assemble_stiffnessMat(operators);
b = assemble_load(operators);

# Boundary condition
g(x::Tensors.Vec{2}) = sin(π * x[2])*sin(π * x[1]);
dbc = Dirichlet(dofs,u,"boundary",g);
apply!(K,b,dbc);

#Solve
x = K\b

#Plot solution
using Makie
import AbstractPlotting
include("src/plot_recipes.jl")
#popdisplay(AbstractPlotting.PlotDisplay())
#AbstractPlotting.inline!(true)
scene = Scene(resolution = (400, 200))
# Plot approximation
vi = jFEMTools.vertexdofs(dofs, u)
plot!(scene, mesh, color = x[vi])

#Plot exact solution
vv = get_vertices_matrix(mesh)
xx = [g(vv) for vv in mesh.vertices];
plot!(scene, mesh, color = xx)
#display(scene)
