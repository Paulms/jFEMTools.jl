@testset "Test Finite Elements" begin

using jFEMTools
import Tensors: Vec, norm, gradient
jF = jFEMTools

# Test ContinuousLagrange Element
finiteElement = ContinuousLagrange(jF.Triangle,2)
@test jF.gettopology(finiteElement) == Dict(0=>3,1=>3,2=>0)
@test jF.getgeomtopology(finiteElement) == Dict(0=>3,1=>0,2=>0)

finiteElement = ContinuousLagrange(jF.Triangle,1)
@test finiteElement.polynomial_space.coeffs ≈ [0.23570226039551578 0.2357022603955159 0.2357022603955158;
                                 -0.14433756729740646 0.14433756729740646 0.0; -1/12 -1/12 1/6]
@test jF.gettopology(finiteElement) == Dict(0=>3,1=>0,2=>0)

# Test Lagrange Continuous Elements
########### 2D in Triangles #######################################
function ref_value(i::Int, ξ::Vec{2})
    ξ_x = ξ[1]
    ξ_y = ξ[2]
    i == 2 && return ξ_x
    i == 3 && return ξ_y
    i == 1 && return 1. - ξ_x - ξ_y
end

ξ_ref = Vec{2}((0.5,0.5))
for j = 1:3
    @test abs(ref_value(j,ξ_ref) - jF.value(finiteElement,j,ξ_ref)) < eps(Float64)
    @test norm(gradient(ξ -> ref_value(j, ξ), ξ_ref) - jF.gradient_value(finiteElement,j,ξ_ref), Inf) < eps(Float64)
end

######################## 1D in Segments #######################################3

function ref_value1(i::Int, ξ::Vec{1})
    ξ_x = ξ[1]
    i == 1 && return 1-ξ_x
    i == 2 && return ξ_x
end
ξ_ref = Vec{1}((0.5,))
finiteElement = ContinuousLagrange(jF.Segment,1)
@test jF.gettopology(finiteElement) == Dict(0=>2,1=>0)
for j = 1:2
    @test abs(ref_value1(j,ξ_ref) - jF.value(finiteElement,j,ξ_ref)) < eps(Float64)
    @test norm(gradient(ξ -> ref_value1(j, ξ), ξ_ref) - jF.gradient_value(finiteElement,j,ξ_ref), Inf) < eps(Float64)
end

############## 3D in Tetrahedra #############################3

finiteElement = ContinuousLagrange(jF.Tetrahedron, 1)
function ref_value3(i::Int, ξ::Vec{3})
    ξ_x = ξ[1]
    ξ_y = ξ[2]
    ξ_z = ξ[3]
    i == 1 && return 1.0 - ξ_x - ξ_y - ξ_z
    i == 2 && return ξ_x
    i == 3 && return ξ_y
    i == 4 && return ξ_z
end

ξ_ref = Vec{3}((0.5,0.5,0.5))
for j = 1:4
    @test abs(ref_value3(j,ξ_ref) - jF.value(finiteElement,j,ξ_ref)) < eps(Float64)
end

########################### 2D in Rectangles #################

finiteElement = ContinuousLagrange(jF.Rectangle, 1)

ξ_ref = Vec{2}((0.5,0.5))
for j = 1:4
    @test abs(0.25 - jF.value(finiteElement,j,ξ_ref)) < eps(Float64)
end

########################### 3D in Hexahedra #################

finiteElement = ContinuousLagrange(jF.Hexahedron, 1)

ξ_ref = Vec{3}((0.5,0.5,0.5))
for j = 1:4
    @test abs(0.125 - jF.value(finiteElement,j,ξ_ref)) < eps(Float64)
end

end
