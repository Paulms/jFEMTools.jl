
abstract type AbstractQuadratureRule end

struct QuadratureRule{shape,N,T} <: AbstractQuadratureRule
    weights::Vector{T}
    points::Vector{Tensors.Vec{N,T}}
end

getweights(quad::QuadratureRule) = quad.weights
getpoints(quad::QuadratureRule) = quad.points

"""
function integrate(f::function,qr::QuadratureRule)
integrate function f on reference shape of quadrature qr
"""
function integrate(f::Function,qr::QuadratureRule)
    p = getpoints(qr); w = getweights(qr)
    int_val = zero(typeof(f(p[1])))
    for (i,x) in enumerate(p)
        int_val += f(x)*w[i]
    end
    int_val
end

#quadratures
struct GrundmannMoeller <: AbstractQuadratureRule end
struct Strang <: AbstractQuadratureRule end
struct DefaultQuad <: AbstractQuadratureRule end
struct GaussLegendre <: AbstractQuadratureRule end
struct Gauss <: AbstractQuadratureRule end

function (::Type{QuadratureRule{Triangle}})(quad_type::DefaultQuad, order::Int)
    if order <= 6
        return QuadratureRule{Triangle}(Strang(),order)
    elseif order > 6 && isodd(order)
        s = Int((order-1)/2)
        return QuadratureRule{Triangle}(GrundmannMoeller(),s)
    else
        throw(ArgumentError("Quadrature rule of order $order not available"))
    end
end

function (::Type{QuadratureRule{Tetrahedron}})(quad_type::DefaultQuad, order::Int)
    if order <= 4
        return QuadratureRule{Tetrahedron}(Gauss(),order)
    elseif order > 6 && isodd(order)
        s = Int((order-1)/2)
        return QuadratureRule{Tetrahedron}(GrundmannMoeller(),s)
    else
        throw(ArgumentError("Quadrature rule of order $order not available"))
    end
end

function (::Type{QuadratureRule{HyperCube{dim}}})(quad_type::DefaultQuad, order::Int) where {dim}
    return QuadratureRule{HyperCube{dim}}(GaussLegendre(),order)
end

# Get GaussLegendre weigths and points from FastGaussQuadrature
function (::Type{QuadratureRule{Segment}})(quad_type::GaussLegendre, order::Int)
    points, weights = FastGaussQuadrature.gausslegendre(order)
    # Shift interval from (-1,1) to (0,1)
    weights *= 0.5
    points = points .+ 1.0; points /= 2.0
    return QuadratureRule{Segment,Float64}(weights, [Tensors.Tensor{1,1}([x]) for x in points])
end

function (::Type{QuadratureRule{Segment}})(quad_type::DefaultQuad, order::Int)
    return QuadratureRule{Segment}(GaussLegendre(), order)
end