getshape(::FiniteElement{dim,shape}) where {dim,shape} = shape

function spatial_nodal_coordinate(mesh::AbstractPolytopalMesh{dim,T}, ci::Int,fe::FiniteElement{dim},n_point::Int) where {dim,T}
    x = getverticescoords(mesh,getcell(mesh,ci))
    n_base_funcs = getngeombasefunctions(fe)
    @assert length(x) == n_base_funcs
    vec = zero(Tensors.Vec{dim,T})
    @inbounds for i in 1:n_base_funcs
        vec += fe.M[i, n_point] * x[i]
    end
    return vec
end

"""
gradient_value(ip::FiniteElement{dim,shape<:Shape,order}, k::Int, ξ::Vec{dim,T}) where {dim,shape,order, T}
Compute value of Finite Element basis `j` derivative at point ξ
on the reference shape `shape`
"""
function gradient_value(ip::FiniteElement{dim}, k::Int, ξ::Tensors.Vec{dim,T}) where {dim, T}
    Tensors.gradient(ξ -> value(ip, k, ξ), ξ)
end
