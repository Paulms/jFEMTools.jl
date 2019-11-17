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

############
# Dubiner #
###########
# struct Dubiner{dim,shape,order} <: Interpolation{dim,shape,order} end
#
# isnodal(ip::Dubiner) = false
#
# getnbasefunctions(::Dubiner{2,RefTetrahedron,order}) where {order} = Int((order+1)*(order+2)/2)
# #nvertexdofs(::Dubiner{2,RefTetrahedron,order}) = 1
#
# #vertices(::Dubiner{2,RefTetrahedron,order}) where {order} = (1,2,3)
# #faces(::Dubiner{2,RefTetrahedron,order}) where {order} = ((1,2), (2,3), (3,1))
#
# """
# value(ip::Dubiner{2,RefTetrahedron,order}, j::Int, ξ::AbstactVector) where {order}
# Compute value of dubiner basis `j` at point ξ
# on the reference triangle ((0,0),(1,0),(0,1))
# """
# function value(ip::Dubiner{2,RefTetrahedron,order}, j::Int, ξ::Vec{2,T}) where {order, T}
#     r = ξ[1]
#     s = ξ[2]
#     return dubiner_basis(r,s,j)
# end
#
# function dubiner_basis(x,y,j::Integer)
#     #Compute degrees
#     t=-3/2+(1/2)*sqrt(1+8*j)
#     n=((ceil(t)+1)*(ceil(t)+2))/2-j
#     m=ceil(t)-n
#     #Compute Dubiner_nm(ξ, η)
#     n = Int(n); m = Int(m)
#     #check domain
#     @assert ((y>=0)&&(x>=0))&&(1>=x+y) "point not in domain"
#     # Map to reference square
#     ξ=abs(x) < eps() ? -one(typeof(x)) : 2*x/(1-y)-1
#     η=2*y-1
#     k=2*n+1
#     #Compute Dubiner_nm(ξ, η)
#     P=jacobi(ξ,n,0,0)*jacobi(η, m,k,0)*((1-η)/2)^n
#     #normalize
#     N=sqrt(2/((2*n+1)*(m+n+1)))
#     return (2*P)/N
# end
#
# function gradient_value(ip::Dubiner{2,RefTetrahedron,order}, k::Int, ξ::Vec{2,T}) where {order,T}
#     if k >getnbasefunctions(ip);throw(ArgumentError("no shape function $k for interpolation $ip"));end
#     gradient(ξ -> value(ip, k, ξ), ξ)
# end
