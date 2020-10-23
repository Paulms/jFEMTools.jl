# jFEMTools

![Lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)
[![Travis](https://travis-ci.org/Paulms/jFEMTools.jl.svg?branch=master)](https://travis-ci.org/Paulms/jFEMTools.jl)
[![AppVeyor](https://ci.appveyor.com/api/projects/status/99gxgpykq2nd7gp8?svg=true)](https://ci.appveyor.com/project/Paulms/jfemtools-jl)
[![Coverage Status](https://coveralls.io/repos/github/Paulms/jFEMTools.jl/badge.svg?branch=master)](https://coveralls.io/github/Paulms/jFEMTools.jl?branch=master)
[![codecov.io](http://codecov.io/github/Paulms/jFEMTools.jl/coverage.svg?branch=master)](http://codecov.io/github/Paulms/jFEMTools.jl?branch=master)

Tools for FEM and VEM code

VEM implementation based on:

- The Hitchhiker's Guide to the Virtual Element Method
L. Beirão da Veiga, F. Brezzi, L. D. Marini, and A. Russo, Mathematical Models and Methods in Applied Sciences 2014 24:08, 1541-1573

FEM implementation based on

- JuaFEM: https://github.com/KristofferC/JuAFEM.jl.git

- Logg, Mardal, Wells, Kirby, FIAT: numerical construction of finite element basis functions.

# FEM
- Continuous Lagrange finite elements of arbitrary order:

```julia
    ContinuousLagrange(:Shape, order)
```

where `:Shape` = :Triangle, :Quad, :Hexahedron, :Tetrahedron

- VTK saving using `WriteVTK.jl` package.

## Example

See `examples` directory for more information. Below, a simple Poisson PDE with homogeneous Dirichlet boundary conditions, in a 3D unitary cube.

```julia
using jFEMTools
import Tensors: Vec,  ⋅
using SparseArrays
import WriteVTK
const jF = jFEMTools

# Problem: 

# Δu = f  in Ω
# u  = 0  in ∂Ω

# We start  generating a simple grid with 10x10x10 Hexahedron elements
# on unitary cube. The generator defaults to the unit hyper_cube.
mesh = hyper_rectagle_mesh2(HexahedronCell,(10,10,10))

# ### Initiate finite element function Space
dim = jF.getdim(mesh)
P1 = ContinuousLagrange(:Hexahedron,1)
Wh = FEMFunctionSpace(mesh, P1, 1)

# Declare variables
u_h = TrialFunction(Wh)
#v_h = TestFunction(Wh)

# ### Degrees of freedom
# We create the `DofHandler` and then add a single field called `u_h`.
dh = DofHandler(mesh,[u_h])

# This function returns a sparse matrix with the correct elements stored.
K = jF.create_sparsity_pattern(dh);

# ### Boundary conditions
dbc = Dirichlet(dh, u_h, "boundary", 0.0)

# ### RHS function
f(x::Vec{dim}) where {dim} = 2*π^2*sin(π*x[1])*sin(π*x[2])*sin(π*x[3])

# ### Assembling the linear system
function doassemble(Wh, K::SparseMatrixCSC, dh::jF.DofHandler)
  # Allocate the element stiffness matrix and element force vector
  # global and local matrices
  n_basefuncs = jF.getnbasefunctions(Wh)
  Ke = zeros(n_basefuncs, n_basefuncs)
  fe = zeros(n_basefuncs)
  b = zeros(jF.ndofs(dh))
  cell_dofs = Vector{Int}(undef, jF.ndofs_per_cell(dh))
  assembler = jF.start_assemble(K, b)

  # It is now time to loop over all the cells in our mesh
  @inbounds for (cellcount, cell) in enumerate(CellIterator(mesh))
      # We recompute local data for each cell
      fill!(Ke, 0)
      fill!(fe, 0)
      jF.reinit!(Wh, cell)

      # Loop over all the quadrature points in the cell and
      # assemble the contribution to `Ke` and `fe`. The integration weight
      # can be queried using `getdetJdV`.
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
K, b = doassemble(Wh, K, dh);

# To account for the boundary conditions we use the `apply!` function.
apply!(K,b,dbc);
u = K \ b;

# Save approximation on vtu file
vi = jFEMTools.vertexdofs(dh, u_h);
vtk_file = vtk_grid("poisson3D", mesh)
vtk_file["u", WriteVTK.VTKPointData()] = u[vi]
outfiles = WriteVTK.vtk_save(vtk_file)
```

# VEM
## Example

We solve a Poisson PDE using Virtual Element Method:

```
using jFEMTools
import Tensors
const jF = jFEMTools

mesh = unitSquareMesh2(RectangleCell, (3,3));

# forcing function
rhs(x::Tensors.Vec{2}) = 2*π^2*sin(π*x[1])*sin(π*x[2])

degree = 2;
dim = jF.getdim(mesh);
element = PoissonVirtualElement(dim,degree);
u = TrialFunction(VEMFunctionSpace(mesh,element))
dofs = DofHandler(mesh, [u]);
operators = VEMOperators(dofs, u;load = rhs);

K = assemble_stiffnessMat(operators);
b = assemble_load(operators);

# Boundary condition
g(x::Tensors.Vec{2}) = sin(π * x[2])*sin(π * x[1]);
dbc = Dirichlet(dofs,u,"boundary",g);
apply!(K,b,dbc);

#Solve
x = K\b

#Plot numerical solution against exact, using Makie
using Makie
import AbstractPlotting
include("src/plot_recipes.jl")
scene = Scene(resolution = (400, 200))
# Plot approximation
vi = jFEMTools.vertexdofs(dofs, u)
plot!(scene, mesh, color = x[vi])

#Plot exact solution
vv = get_vertices_matrix(mesh)
xx = [g(vv) for vv in mesh.vertices];
plot!(scene, mesh, color = xx)
```

