using jFEMTools
import Tensors: Vec,  ⋅
using SparseArrays
import WriteVTK
const jF = jFEMTools


# We solve a simple 3D Poisson PDE with homogeneous Dirichlet boundary conditions
# using ContinuousLagrange finite elements on a Hexahedral mesh

# Problem: 

# Δu = f  in Ω
# u  = 0  in ∂Ω

# We start  generating a simple grid with 20x20 triangular elements
# using `unitSquareMesh2`. The generator defaults to the unit square,
# so we don't need to specify the corners of the domain.
mesh = hyper_rectagle_mesh2(HexahedronCell,(10,10,10))

# ### Initiate function Spaces
dim = jF.getdim(mesh)
P1 = ContinuousLagrange(:Hexahedron,1)
Wh = FEMFunctionSpace(mesh, P1, 1)

# Declare variables
u_h = TrialFunction(Wh)
#v_h = TestFunction(Wh)

# ### Degrees of freedom
# Next we need to define a `DofHandler`, which will take care of numbering
# and distribution of degrees of freedom for our approximated fields.
# We create the `DofHandler` and then add a single field called `u_h`.
dh = DofHandler(mesh,[u_h])

# Now that we have distributed all our dofs we can create our tangent matrix,
# using `create_sparsity_pattern`. This function returns a sparse matrix
# with the correct elements stored.
K = jF.create_sparsity_pattern(dh);

# ### Boundary conditions
dbc = Dirichlet(dh, u_h, "boundary", 0.0)

# ### RHS function
f(x::Vec{dim}) where {dim} = 2*π^2*sin(π*x[1])*sin(π*x[2])*sin(π*x[3])

# ### Assembling the linear system
# Now we have all the pieces needed to assemble the linear system, $K u = f$.
function doassemble(Wh, K::SparseMatrixCSC, dh::jF.DofHandler)
  # We allocate the element stiffness matrix and element force vector
  # just once before looping over all the cells instead of allocating
  # them every time in the loop.
  n_basefuncs = jF.getnbasefunctions(Wh)
  Ke = zeros(n_basefuncs, n_basefuncs)
  fe = zeros(n_basefuncs)
  b = zeros(jF.ndofs(dh))
  cell_dofs = Vector{Int}(undef, jF.ndofs_per_cell(dh))

  # Next we define the global force `f` and
  # create an assembler. The assembler
  # is just a thin wrapper around `f` and `K` and some extra storage
  # to make the assembling faster.
  assembler = jF.start_assemble(K, b)

  # It is now time to loop over all the cells in our mesh
  @inbounds for (cellcount, cell) in enumerate(CellIterator(mesh))
      # Always remember to reset the function space local data
      fill!(Ke, 0)
      fill!(fe, 0)
      jF.reinit!(Wh, cell)

      # It is now time to loop over all the quadrature points in the cell and
      # assemble the contribution to `Ke` and `fe`. The integration weight
      # can be queried from `cellvalues` by `getdetJdV`.
      for q_point in 1:jF.getnquadpoints(Wh)
          dΩ = jF.getdetJdV(Wh, q_point)
          fh = jF.function_value(f, Wh, cell, q_point)
          # For each quadrature point we loop over all the (local) shape functions.
          # We need the value and gradient of the testfunction `v` and also the gradient
          # of the trial function `u`.
          for i in 1:n_basefuncs
              v  = jF.shape_value(Wh, q_point, i)
              ∇v = jF.shape_gradient(Wh, q_point, i)
              fe[i] += fh*v * dΩ
              for j in 1:n_basefuncs
                  ∇u = jF.shape_gradient(Wh, q_point, j)
                  Ke[i, j] += (∇v ⋅ ∇u) * dΩ
              end
          end
      end
      # The last step in the element loop is to assemble `Ke` and `fe`
      # into the global `K` and `f` with `assemble!`.
      jF.celldofs!(cell_dofs, dh, cell)
      jF.assemble!(assembler, cell_dofs, fe, Ke)
  end
  return K, b
end

# ### Solution of the system
# The last step is to solve the system. First we call `doassemble`
# to obtain the global stiffness matrix `K` and force vector `f`.
K, b = doassemble(Wh, K, dh);

# To account for the boundary conditions we use the `apply!` function.
# This modifies elements in `K` and `f` respectively, such that
# we can get the correct solution vector `u` by using `\`.
apply!(K,b,dbc);
u = K \ b;

# Plot approximation
vi = jFEMTools.vertexdofs(dh, u_h);
vtk_file = vtk_grid("poisson3D", mesh)
vtk_file["u", WriteVTK.VTKPointData()] = u[vi]
outfiles = WriteVTK.vtk_save(vtk_file)
