struct Monomials{dim,order, T}
    indices::Union{Vector{NTuple{dim,Int}},Vector{Int}}
    centroid::Tensors.Vec{dim,T}
    diameter::T
end
function Monomials(dim,degree, centroid, diameter)
    if dim == 1
        return Monomials1D(degree, centroid, diameter)
    elseif dim == 2
        return Monomials2D(degree, centroid, diameter)
    else
        throw("Monomials of dim $dim are not implemented!")
    end
end
function Monomials1D(degree::Int, centroid::Tensors.Vec{1,T}, diameter::T) where {T}
    indices = 0:degree
    Monomials{1,degree, T}(collect(indices), centroid, diameter)
end
function Monomials2D(degree::Int, centroid::Tensors.Vec{2,T}, diameter::T) where {T}
    indices = get_monomial2DIndices(degree)
    Monomials{2,degree, T}(indices, centroid, diameter)
end

function get_monomial2DIndices(degree)
    if degree == 0
        return [(0,0)]
    else
        tmp = Vector{NTuple{2,Int}}()
        a = (degree,0)
        for i in 1:(binomial(degree+1,1))
            push!(tmp,a)
            a = (a[1]-1,a[2]+1)
        end
        return [get_monomial2DIndices(degree-1)...,tmp...]
    end
end


getnbasefunctions(::Monomials{2,order}) where {order} = Int((order+1)*(order+2)/2)
getnbasefunctions(::Monomials{1,order}) where {order} = order+1


# Evaluate monomials
function value(ip::Monomials{dim,order}, k::Int, ξ::Tensors.Vec{dim,T}) where {order, T, dim}
    if k > getnbasefunctions(ip);throw(ArgumentError("no shape function $k for interpolation $ip"));end
    prod((1/ip.diameter*(ξ-ip.centroid)).^ip.indices[k])
end

function gradient_value(ip::Monomials{dim,order}, k::Int, ξ::Tensors.Vec{dim,T}) where {dim,order,T}
    if k >getnbasefunctions(ip);throw(ArgumentError("no shape function $k for interpolation $ip"));end
    Tensors.gradient(ξ -> value(ip, k, ξ), ξ)
end

function laplace_value(ip::Monomials{dim,order}, k::Int, ξ::Tensors.Vec{dim,T}) where {dim,order,T}
    if k >getnbasefunctions(ip);throw(ArgumentError("no shape function $k for interpolation $ip"));end
    Tensors.laplace(ξ -> value(ip, k, ξ), ξ)
end


struct RotMonomials{dim,order, T}
    indices::Union{Vector{NTuple{dim,Int}},Vector{Int}}
    centroid::Tensors.Vec{dim,T}
    diameter::T
end
function RotMonomials(dim,order, centroid, diameter)
    if dim == 2
        return RotMonomials2D(order, centroid, diameter)
    else
        throw("Monomials of dim $dim are not implemented!")
    end
end
function RotMonomials2D(order::Int, centroid::Tensors.Vec{2,T}, diameter::T) where {T}
    indices = get_monomial2DIndices(order-1)
    RotMonomials{2,order, T}(indices, centroid, diameter)
end


getnbasefunctions(::RotMonomials{2,order}) where {order} = Int((order)*(order+1)/2)


# Evaluate monomials

rot(ξ::Tensors.Vec{2}) = Tensors.Vec((-ξ[2],ξ[1])) 

function value(ip::RotMonomials{dim,order}, k::Int, ξ::Tensors.Vec{dim,T}) where {order, T, dim}
    if k > getnbasefunctions(ip);throw(ArgumentError("no shape function $k for interpolation $ip"));end
    rot(ξ)*prod((1/ip.diameter*(ξ-ip.centroid)).^ip.indices[k])
end

function gradient_value(ip::RotMonomials{dim,order}, k::Int, ξ::Tensors.Vec{dim,T}) where {dim,order,T}
    if k >getnbasefunctions(ip);throw(ArgumentError("no shape function $k for interpolation $ip"));end
    Tensors.gradient(ξ -> value(ip, k, ξ), ξ)
end

function laplace_value(ip::RotMonomials{dim,order}, k::Int, ξ::Tensors.Vec{dim,T}) where {dim,order,T}
    if k >getnbasefunctions(ip);throw(ArgumentError("no shape function $k for interpolation $ip"));end
    Tensors.laplace(ξ -> value(ip, k, ξ), ξ)
end

