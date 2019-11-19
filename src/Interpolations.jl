struct Monomials{dim,order, T}
    indices::Vector{Tuple{Int,Int}}
    centroid::Tensors.Vec{dim,T}
    diameter::T
end
function Monomials(dim,degree, centroid, diameter)
    if dim == 2
        return Monomials2D(degree, centroid, diameter)
    else
        throw("Monomials of dim $dim are not implemented!")
    end
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

function value(ip::Monomials{2,order}, k::Int, ξ::Tensors.Vec{2,T}) where {order, T}
    if k > getnbasefunctions(ip);throw(ArgumentError("no shape function $k for interpolation $ip"));end
    prod((1/ip.diameter*(ξ-ip.centroid)).^ip.indices[k])
end

function gradient_value(ip::Monomials{2,order}, k::Int, ξ::Tensors.Vec{dim,T}) where {dim,order,T}
    if k >getnbasefunctions(ip);throw(ArgumentError("no shape function $k for interpolation $ip"));end
    Tensors.gradient(ξ -> value(ip, k, ξ), ξ)
end

function laplace_value(ip::Monomials{2,order}, k::Int, ξ::Tensors.Vec{dim,T}) where {dim,order,T}
    if k >getnbasefunctions(ip);throw(ArgumentError("no shape function $k for interpolation $ip"));end
    Tensors.laplace(ξ -> value(ip, k, ξ), ξ)
end
