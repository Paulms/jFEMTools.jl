abstract type Interpolation{shape,order} end

"""
Return the dimension of an `Interpolation`
"""
@inline getdim(ip::Interpolation{shape}) where {shape} = getdim(shape)

"""
Return the reference shape of an `Interpolation`
"""
@inline getrefshape(ip::Interpolation{shape}) where {shape} = shape

"""
Return the polynomial order of the `Interpolation`
"""
@inline getorder(ip::Interpolation{shape,order}) where {shape,order} = order

"""
Return the topology Dict of the `Interpolation`
"""
@inline gettopology(ip::Interpolation) = Dict{Int,Int}()

"""
Compute the value of the shape functions at a point ξ for a given interpolation
"""
function value(ip::Interpolation, ξ::Tensors.Vec{dim,T}) where {dim,T}
    [value(ip, i, ξ) for i in 1:getnbasefunctions(ip)]
end

function value(ip::Interpolation{shape}, j::Int, ξ::Tuple) where {shape}
    value(ip,j,Tensors.Vec{getdim(shape)}(ξ))
end

"""
Compute the gradients of the shape functions at a point ξ for a given interpolation
"""
function gradient_value(ip::Interpolation, ξ::Tensors.Vec{dim,T}) where {dim,T}
    [gradient_value(ip, i, ξ) for i in 1:getnbasefunctions(ip)]
end

function gradient_value(ip::Interpolation{shape}, j::Int, ξ::Tuple) where {shape}
    gradient_value(ip,j,Tensors.Vec{getdim(shape)}(ξ))
end

# Default to nodal interpolator for geom
function get_default_geom_interpolator(shape::Shape)
    return Lagrange(shape,1)
end

############
# Dubiner #
###########
struct Dubiner{shape,order} <: Interpolation{shape,order} end
isnodal(ip::Dubiner) = false

getnbasefunctions(::Dubiner{Triangle,order}) where {order} = Int((order+1)*(order+2)/2)

"""
value(ip::Dubiner{2,RefTetrahedron,order}, j::Int, ξ::AbstactVector) where {order}
Compute value of dubiner basis `j` at point ξ
on the reference triangle ((0,0),(1,0),(0,1))
"""
function value(ip::Dubiner{Triangle,order}, j::Int, ξ::Tensors.Vec{2,T}) where {order, T}
    r = ξ[1]
    s = ξ[2]
    if j == 0
        return zero(T)
    else 
        return dubiner_basis(r,s,j)
    end
end

"""
gradient_value(ip::Dubiner{2,RefTetrahedron,order}, j::Int, ξ::AbstactVector) where {order}
Compute gradient of dubiner basis `j` at point ξ
on the reference triangle ((0,0),(1,0),(0,1))
"""
function gradient_value(ip::Dubiner{Triangle,order}, j::Int, ξ::Tensors.Vec{2,T}) where {order, T}
    if j >getnbasefunctions(ip);throw(ArgumentError("no shape function $j for interpolation $ip"));end
    Tensors.gradient(ξ -> value(ip, j, ξ), ξ)
end


"""
    jacobi(x, p::Integer, α, β)
Evaluate the Legendre polynomial with parameters `α`, `β` of degree `p` at `x`
using the three term recursion [Karniadakis and Sherwin, Spectral/hp Element
Methods for CFD, Appendix A].
Author: H. Ranocha (see PolynomialBases.jl)
"""
function jacobi(x, p::Integer, α, β)
    T = typeof( (2+α+β)*x / 2 )
    a = one(T)
    b = ((2+α+β)*x + α - β) / 2
    if p <= 0
        return a
    elseif p == 1
        return b
    end

    for n in 2:p
        a1 = 2n*(n+α+β)*(2n-2+α+β)
        a2 = (2n-1+α+β)*(α+β)*(α-β)
        a3 = (2n-2+α+β)*(2n-1+α+β)*(2n+α+β)
        a4 = 2*(n-1+α)*(n-1+β)*(2n+α+β)
        a, b = b, ( (a2+a3*x)*b - a4*a ) / a1
    end
    b
end

"""
dubiner(x,y,n::Integer,m::Integer)
Compute the dubiner polynomial of degree `n`,`m` at point (x,y)
on the reference triangle ((0,0),(1,0),(0,1))
"""
function dubiner(x,y,n::Integer,m::Integer)
    #check domain
    @assert ((y>=0)&&(x>=0))&&(1>=x+y) "point not in domain"
    # Map to reference square
    ξ=abs(x) < eps() ? -one(typeof(x)) : 2*x/(1-y)-1
    η=2*y-1
    k=2*n+1
    #Compute Dubiner_nm(ξ, η)
    P=jacobi(ξ,n,0,0)*jacobi(η, m,k,0)*((1-η)/2)^n
    #normalize
    N=sqrt(2/((2*n+1)*(m+n+1)))
    return (2*P)/N
end

"""
dubiner_basis(x,y,j::Integer)
Evaluate the dubiner basis `j` at point (x,y)
on the reference triangle ((0,0),(1,0),(0,1))
"""
function dubiner_basis(x,y,j::Integer)
    #Compute degrees
    t=-3/2+(1/2)*sqrt(1+8*j)
    n=((ceil(t)+1)*(ceil(t)+2))/2-j
    m=ceil(t)-n
    #Compute Dubiner_nm(ξ, η)
    return dubiner(x,y,Int(n),Int(m))
end

####################
# Lagrange
####################
struct Lagrange{shape,order} <: Interpolation{shape,order}
    nodal_base_coefs::Matrix{Float64}
    topology::Dict{Int,Int}
end

isnodal(ip::Lagrange) = true

function getdefaultdualbasis(shape::Type{s},order::Int) where {s<:Shape}
    if getdim(shape()) == 1
        return Legendre{shape,order}()
    elseif getdim(shape()) == 2
        return Dubiner{shape,order}()
    else
        throw("Not dual basis available for dimension $dim")
    end
end

function Lagrange(shape,order::Int)
    nodal_points, topology = get_nodal_points(shape, order)
    nodals=[x->x(point) for point in nodal_points]
    ip_prime = getdefaultdualbasis(shape,order)
    nbasefuncs = getnbasefunctions(ip_prime)
    prime_base = [x->value(ip_prime, j, x) for j in 1:nbasefuncs]
    V = reshape([nodals[i](prime_base[j]) for j = 1:nbasefuncs for i=1:nbasefuncs],(nbasefuncs,nbasefuncs))
    nodal_base_coefs = inv(V)
    Lagrange{shape, order}(nodal_base_coefs, topology)
end

@inline getnbasefunctions(::Lagrange{Segment,order}) where {order} = order + 1
@inline getnbasefunctions(::Lagrange{Triangle,order}) where {order} = Int((order+1)*(order+2)/2)
@inline gettopology(ip::Lagrange) = ip.topology

"""
value(ip::Lagrange{shape,order}, k::Int, ξ::Tensors.Vec{dim,T}) where {shape,order, dim,T}
Compute value of Lagrange basis `j` at point ξ
on the reference shape
"""
function value(ip::Lagrange{shape,order}, k::Int, ξ::Tensors.Vec{dim,T}) where {shape,order, dim,T}
    if k > getnbasefunctions(ip);throw(ArgumentError("no shape function $k for interpolation $ip"));end
    Tensors.dot(ip.nodal_base_coefs[:,k], value(getdefaultdualbasis(shape,order), ξ))
end

"""
gradient_value(ip::Lagrange{shape,order}, k::Int, ξ::Tensors.Vec{dim,T}) where {dim,shape,order,T}
Compute value of Lagrange basis `j` gradient at point ξ
on the reference shape
"""
function gradient_value(ip::Lagrange{shape,order}, k::Int, ξ::Tensors.Vec{dim,T}) where {dim,shape,order,T}
    if k >getnbasefunctions(ip);throw(ArgumentError("no shape function $k for interpolation $ip"));end
    Tensors.gradient(ξ -> value(ip, k, ξ), ξ)
end

function _get_nodal_transformation_matrix(fe::Lagrange{shape,order}) where {shape,order}
    # Matrix to get spacial coordinates
    nodal_points, topology = get_nodal_points(shape, order)
    T = eltype(nodal_points[1])
    geom_interpol = get_default_geom_interpolator(shape)
    qrs = QuadratureRule{shape,getdim(shape),T}(fill(T(NaN), length(nodal_points)), nodal_points) # weights will not be used
    n_qpoints = length(getweights(qrs))
    n_geom_basefuncs = getnbasefunctions(geom_interpol)
    M =    fill(zero(T)           * T(NaN), n_geom_basefuncs, n_qpoints)
    for (qp, ξ) in enumerate(qrs.points)
        for i in 1:n_geom_basefuncs
            M[i, qp] = value(geom_interpol, i, ξ)
        end
    end
    M
end

function spatial_nodal_coordinate(fe::Lagrange, M::Matrix{T}, n_point::Int, x::AbstractVector{Tensors.Vec{dim,T}}) where {dim,T}
    n_base_funcs = size(M,1)
    @assert length(x) == n_base_funcs
    vec = zero(Tensors.Vec{dim,T})
    @inbounds for i in 1:n_base_funcs
        vec += M[i, n_point] * x[i]
    end
    return vec
end

####################
# Legendre
####################
struct Legendre{shape,order} <: Interpolation{shape,order} end
isnodal(Legendre) = false

@inline getnbasefunctions(::Legendre{Segment,order}) where {order} = order + 1

"""
value(ip::Legendre{Segment,order}, j::Int, ξ::AbstactVector) where {order}
Compute value of Legendre basis `j` at point ξ
on the reference line (0,1)
"""
function value(ip::Legendre{Segment,order}, k::Int, ξ::Tensors.Vec{1,T}) where {order, T}
    if k > getnbasefunctions(ip);throw(ArgumentError("no shape function $k for interpolation $ip"));end
    return sqrt((2*(k-1)+1))*jacobi(2*ξ[1]-1,k-1,0.0,0.0)
end

"""
gradient_value(ip::Legendre{Segment,order}, j::Int, ξ::AbstactVector) where {order}
Compute value of Legendre basis `j` derivative at point ξ
on the reference line (0,1)
"""
function gradient_value(ip::Legendre{Segment,order}, k::Int, ξ::Tensors.Vec{1,T}) where {order, T}
    if k >getnbasefunctions(ip);throw(ArgumentError("no shape function $k for interpolation $ip"));end
    Tensors.gradient(ξ -> value(ip, k, ξ), ξ)
end
