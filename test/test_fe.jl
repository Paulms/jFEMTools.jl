@testset "Test Finite Elements" begin

using jFEMTools
import Tensors: Vec, norm, gradient
jF = jFEMTools

# Test ContinuousLagrange Element
finiteElement = ContinuousLagrange(jF.Triangle,2)
@test jF.gettopology(finiteElement) == Dict(0=>3,1=>3,2=>0)
@test jF.getgeomtopology(finiteElement) == Dict(0=>3,1=>0,2=>0)
x = [Vec{2}((1.0,1.0)),Vec{2}((2.0,1.0)),Vec{2}((2.0,2.0))]
y = [Vec{2}((1.0,1.0)),Vec{2}((2.0,1.0)),Vec{2}((2.0,2.0)),Vec{2}((1.5,1.5)),Vec{2}((1.5,1.0)),Vec{2}((1.0,1.5))]
for i in 1:3
    @test jF.spatial_nodal_coordinate(finiteElement, i, x) == y[i]
end

finiteElement = ContinuousLagrange(jF.Triangle,1)
@test finiteElement.nodal_base_coefs ≈ [0.23570226039551578 0.2357022603955159 0.2357022603955158;
                                 -0.14433756729740646 0.14433756729740646 0.0; -1/12 -1/12 1/6]
@test jF.gettopology(finiteElement) == Dict(0=>3,1=>0,2=>0)
for i in 1:3
    @test jF.spatial_nodal_coordinate(finiteElement, i, x) == x[i]
end

# Test Lagrange Continuous Elements
function ref_value(i::Int, ξ::Vec{2})
    ξ_x = ξ[1]
    ξ_y = ξ[2]
    i == 2 && return ξ_x
    i == 3 && return ξ_y
    i == 1 && return 1. - ξ_x - ξ_y
    throw(ArgumentError("no shape function $i for interpolation $ip"))
end

ξ_ref = Vec{2}((0.5,0.5))
for j = 1:3
    @test abs(ref_value(j,ξ_ref) - jF.value(finiteElement,j,ξ_ref)) < eps(Float64)
    @test norm(gradient(ξ -> ref_value(j, ξ), ξ_ref) - jF.gradient_value(finiteElement,j,ξ_ref), Inf) < eps(Float64)
end

function ref_value1(i::Int, ξ::Vec{1})
    ξ_x = ξ[1]
    i == 1 && return 1-ξ_x
    i == 2 && return ξ_x
    throw(ArgumentError("no shape function $i for interpolation $ip"))
end
ξ_ref = Vec{1}((0.5,))
finiteElement = ContinuousLagrange(jF.Segment,1)
@test jF.gettopology(finiteElement) == Dict(0=>2,1=>0)
for j = 1:2
    @test abs(ref_value1(j,ξ_ref) - jF.value(finiteElement,j,ξ_ref)) < eps(Float64)
    @test norm(gradient(ξ -> ref_value1(j, ξ), ξ_ref) - jF.gradient_value(finiteElement,j,ξ_ref), Inf) < eps(Float64)
end

end
