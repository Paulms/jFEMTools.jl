using jFEMTools
using Tensors
using SparseArrays

# We start  generating a simple grid with 20x20 quadrilateral elements
# using `generate_grid`. The generator defaults to the unit square,
# so we don't need to specify the corners of the domain.
mesh = unitSquareMesh2(TriangleCell, (3,3));

# ### Initiate function Spaces
dim = 2
P1 = ContinuousLagrange(:Triangle,1)
Wh = FEFunctionSpace(mesh, P1, dim)

# Declare variables
u_h = TrialFunction(Wh)
v_h = TestFunction(Wh)

# ### Degrees of freedom
# Next we need to define a `DofHandler`, which will take care of numbering
# and distribution of degrees of freedom for our approximated fields.
# We create the `DofHandler` and then add a single field called `u`.
dh = DofHandler(mesh,u_h);

# Now that we have distributed all our dofs we can create our tangent matrix,
# using `create_sparsity_pattern`. This function returns a sparse matrix
# with the correct elements stored.
K = create_sparsity_pattern(dh);

# ### Boundary conditions
dbc = Dirichlet(dh, u_h, "boundary", x->0.0)

# ### RHS function
f(x::Vec{dim}) = 2*π^2*sin(π*x[1])*sin(π*x[2])

# ### Assembling the linear system
K_f = Int2D(mesh,dot(grad(u_h),grad(v_h)))
b_f = Int2D(mesh,dot(f,v_h))

K = assemble(K_f, dh)
b = assemble(b_f, dh)

# To account for the boundary conditions we use the `apply!` function.
# This modifies elements in `K` and `b` respectively, such that
# we can get the correct solution vector `u` by using `\`.
apply!(K,b,dbc);
u = K \ b