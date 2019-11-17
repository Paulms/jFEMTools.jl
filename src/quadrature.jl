
abstract type AbstractShape end

struct Segment <: AbstractShape end
struct RefSimplex <: AbstractShape end

function reference_coordinates(::Type{RefSimplex}, ::Type{Val{1}})
    return [Tensors.Vec{1, Float64}((0.0,)),Tensors.Vec{1, Float64}((1.0,))]
end

function reference_coordinates(::Type{RefSimplex}, ::Type{Val{2}})
    [Tensors.Vec{2, Float64}((0.0, 0.0)),
     Tensors.Vec{2, Float64}((1.0, 0.0)),
     Tensors.Vec{2, Float64}((0.0, 1.0))]
end

abstract type AbstractQuadratureRule end

struct QuadratureRule{dim,shape,T} <: AbstractQuadratureRule
    weights::Vector{T}
    points::Vector{Tensors.Vec{dim,T}}
end
