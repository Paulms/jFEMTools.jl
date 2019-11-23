using jFEMTools
import Tensors

mesh = unitSquareMesh(RectangleCell, (3,3));

# forcing function
rhs(x::Tensors.Vec{2}) = 2*π^2*sin(π*x[1])*sin(π*x[2])

degree = 2;
dim = 2;
element = VirtualElement(dim,degree);
dofs = DofHandler(mesh, element);
operators = VEMOperators(dofs, element;load = rhs);

K = assemble_K(operators)
b = assemble_load(operators);

# Boundary condition
g(x::Tensors.Vec{2}) = sin(π * x[2])*sin(π * x[1]);
dbc = Dirichlet(dofs,element,"boundary",g);
apply!(K,b,dbc);

#Solve
x = K\b
vi = jFEMTools.vertexdofs(dofs, element)

# Compute exact solution
nv = getnvertices(mesh)
vv = get_vertices_matrix(mesh)
xx = [g(vv.x) for vv in mesh.vertices];

#Plot solution
using Makie
#popdisplay(AbstractPlotting.PlotDisplay())
#AbstractPlotting.inline!(true)
scene = Scene(resolution = (400, 200))
plot!(scene, mesh, color = x[vi])
plot!(scene, mesh, color = xx)
#display(scene)
