@testset "Test Function Spaces" begin

using jFEMTools
import Tensors: Vec, gradient, dot, norm, Tensor
using LinearAlgebra
jF = jFEMTools

# Load mesh
mesh = unitSquareMesh2(TriangleCell, (2,2));

# ### Define Function Spaces
dim = 2
finiteElement = ContinuousLagrange(:Triangle,1)
Wh = FEMFunctionSpace(mesh, finiteElement, 1)    #Scalar Space
Vh = FEMFunctionSpace(mesh, finiteElement, 2)    #Vector Space

# Basic Test
@test jF.getnlocaldofs(Vh) == 6
@test jF.getnlocaldofs(Wh) == 3
@inbounds for (cellcount, cell) in enumerate(jF.CellIterator(mesh))
    jF.reinit!(Wh, cell)
    jF.reinit!(Vh, cell)
end

# Wh Scalar Space testset
@test jF.getnquadpoints(Wh) == 3
@test jF.getnbasefunctions(Wh) == 3
@test jF.getdetJdV(Wh, 1) ≈ 1/24
@test jF.geometric_value(Wh, 1, 1) ≈ 2/3
@test jF.getdim(Wh) == 2
@test jF.shape_value(Wh,1,1) ≈ 2/3
@test jF.shape_gradient(Wh,1,1) ≈ Vec{2}((0.0,-2.0))
@test jF.shape_divergence(Wh,1,1) ≈ -2.0

# Vh Vector Space testset
@test jF.getnquadpoints(Vh) == 3
@test jF.getnbasefunctions(Vh) == 6
@test jF.getdetJdV(Vh, 1) ≈ 1/24
@test jF.geometric_value(Vh, 1, 1) ≈ 2/3
@test jF.getdim(Vh) == 2
@test jF.shape_value(Vh,1,1) ≈ [2/3, 0.0]
@test jF.shape_gradient(Vh,1,1) ≈ Tensor{2,2}([0.0 -2.0;0.0 0.0])
@test jF.shape_divergence(Vh,1,2) ≈ 2.0

end
