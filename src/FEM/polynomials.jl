getpolynomialbasisdim(order::Int, ::Type{Triangle}) = Int((order+1)*(order+2)/2)
getpolynomialbasisdim(order::Int, ::Type{Tetrahedron}) = Int((order+1)*(order+2)*(order+3)/6)
getpolynomialbasisdim(order::Int, ::Type{HyperCube{dim}}) where dim = (order+1)^dim

struct PolynomialSpace{dim,shape,order}
    polynomial_basis :: Interpolation{shape,order}
    degree::Int
    coeffs::Matrix{Float64}
end

value(npb::PolynomialSpace{dim}, k::Int, ξ::Tensors.Vec{dim,T}) where {dim, T} = 
    Tensors.dot(npb.coeffs[:,k], value(npb.polynomial_basis, ξ))

getnbasefunctions(nps::PolynomialSpace) = getnbasefunctions(nps.polynomial_basis)