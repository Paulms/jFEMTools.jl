struct Monomials{dim,order}
    indices::Vector{Tuple{Int,Int}}
end
function Monomials(dim,degree)
    if dim == 2
        return Monomials2D(degree)
    else
        throw("Monomials of dim $dim are not implemented!")
    end
end
function Monomials2D(degree::Int)
    indices = Vector{Tuple{Int,Int}}()
    for (i,j) in zip(repeat(0:2,inner=3),repeat(0:2,outer=3))
        if i + j < degree+1
            push!(indices,(i,j))
        end
    end
    Monomials{2,degree}(indices)
end

getnbasefunctions(::Monomials{2,order}) where {order} = Int((order+1)*(order+2)/2)

function value(ip::Monomials{2,order}, k::Int, ξ::Tensors.Vec{2,T}) where {order, T}
    if k > getnbasefunctions(ip);throw(ArgumentError("no shape function $k for interpolation $ip"));end
    sum(ξ.^ip.indices[k])
end

function gradient_value(ip::Monomials{2,order}, k::Int, ξ::Tensors.Vec{dim,T}) where {dim,order,T}
    if k >getnbasefunctions(ip);throw(ArgumentError("no shape function $k for interpolation $ip"));end
    gradient(ξ -> value(ip, k, ξ), ξ)
end

function laplace_value(ip::Monomials{2,order}, k::Int, ξ::Tensors.Vec{dim,T}) where {dim,order,T}
    if k >getnbasefunctions(ip);throw(ArgumentError("no shape function $k for interpolation $ip"));end
    laplace(ξ -> value(ip, k, ξ), ξ)
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
