mutable struct FEFSJacobians{dim,T,MM}
    detJ::T
    Jinv::Tensors.Tensor{2,dim,T,MM}
end

struct FEMFunctionSpace{dim,T<:Real,FE<:FiniteElement,mtype <: AbstractPolytopalMesh,MM} <: AbstractFEMFunctionSpace{dim,T,FE}
    N::Matrix{T}
    dNdx::Matrix{Vec{dim,T}}
    dNdξ::Matrix{Vec{dim,T}}
    M::Matrix{T}
    dMdξ::Matrix{Vec{dim,T}}
    qr_weights::Vector{T}
    fe::FE
    mesh::mtype
    jacobians::FEFSJacobians{dim,T,MM}
    components::Int
end

getelement(fs::FEMFunctionSpace) = fs.fe

function Base.show(io::IO, fs::FEMFunctionSpace)
    if fs.components == 1
        println(io, "Discrete scalar function space")
    else
        println(io, "Discrete vectorial function space (components = $(fs.components))")
    end
    println(io, "Finite Element: \n", fs.fe)
    println(io, "Domain mesh: \n", fs.mesh)
    println(io, "Number of local dofs: ", getnlocaldofs(fs))
  end

function FEMFunctionSpace(mesh::AbstractPolytopalMesh, felem::FiniteElement{dim,shape,order,gorder}, components::Int;
    quad_degree::Int = order+1) where {dim, shape, order,gorder}
    quad_rule = QuadratureRule{shape}(DefaultQuad(), quad_degree)
    _scalar_fs(Float64, mesh, quad_rule, felem, components)
end

function _scalar_fs(::Type{T}, mesh::AbstractPolytopalMesh{dim,T}, quad_rule::QuadratureRule{shape},
    felem::FiniteElement{dim,shape,order,1}, components::Int) where {dim,T,shape<:Shape,order}
    n_qpoints = length(getpoints(quad_rule))
    n_cells = getncells(mesh)

    # Function interpolation
    n_func_basefuncs = getnbasefunctions(felem)
    N    = fill(zero(T)          * T(NaN), n_func_basefuncs, n_qpoints)
    dNdx = fill(zero(Vec{dim,T}) * T(NaN), n_func_basefuncs, n_qpoints)
    dNdξ = fill(zero(Vec{dim,T}) * T(NaN), n_func_basefuncs, n_qpoints)

    # Geometry interpolation
    n_geom_basefuncs = getngeombasefunctions(felem)
    M    = fill(zero(T)          * T(NaN), n_geom_basefuncs, n_qpoints)
    dMdξ = fill(zero(Vec{dim,T}) * T(NaN), n_geom_basefuncs, n_qpoints)

    for (qp, ξ) in enumerate(quad_rule.points)
        for i in 1:n_func_basefuncs
            dNdξ[i, qp], N[i, qp]  = gradient(ξ -> value(felem, i, ξ), ξ, :all)
        end
        for i in 1:n_geom_basefuncs
            dMdξ[i, qp], M[i, qp] = gradient(ξ -> geom_value(felem, i, ξ), ξ, :all)
        end
    end

    detJ = T(NaN)
    Jinv = zero(Tensors.Tensor{2,dim,T,2*dim}) 

    MM = Tensors.n_components(Tensors.get_base(typeof(Jinv)))

    jacobians = FEFSJacobians{dim,T,MM}(detJ,Jinv)
    FEMFunctionSpace{dim,T,typeof(felem),typeof(mesh),MM}(N, dNdx, dNdξ,
    M, dMdξ, getweights(quad_rule), felem, mesh,jacobians,components)
end

function reinit!(fs::FEMFunctionSpace{dim}, x::AbstractVector{Vec{dim,T}}) where {dim,T}
    n_geom_basefuncs = getngeombasefunctions(fs.fe)
    n_func_basefuncs = getnbasefunctions(fs.fe)
    @assert length(x) == n_geom_basefuncs
    fecv_J = zero(Tensors.Tensor{2,dim})
    for j in 1:n_geom_basefuncs
        fecv_J += x[j] ⊗ fs.dMdξ[j, 1]
    end
    detJ = det(fecv_J)
    detJ > 0.0 || throw(ArgumentError("det(J) is not positive: det(J) = $(detJ)"))
    fs.jacobians.detJ = detJ
    fs.jacobians.Jinv = inv(fecv_J)

    for j in 1:n_func_basefuncs, i in 1:length(fs.qr_weights)
        fs.dNdx[j, i] = fs.dNdξ[j, i] ⋅ fs.jacobians.Jinv
    end
end

########### Data Functions
getngeobasefunctions(fs::AbstractFEMFunctionSpace) = size(fs.M, 1)
getnquadpoints(fs::AbstractFEMFunctionSpace) = length(fs.qr_weights)
getnbasefunctions(fs::AbstractFEMFunctionSpace) = size(fs.N,1)*fs.components
getdetJdV(fs::FEMFunctionSpace, q_point::Int) = fs.jacobians.detJ*fs.qr_weights[q_point]
geometric_value(fs::AbstractFEMFunctionSpace, q_point::Int, base_func::Int) = fs.M[base_func, q_point]
getdim(::AbstractFEMFunctionSpace{dim}) where {dim} = dim
reference_coordinate(fs::AbstractFEMFunctionSpace{dim,T},cell::Int, mesh, x::Vec{dim,T}) where {dim,T} = fs.jacpbians.Jinv⋅(x-mesh.nodes[mesh.cells[cell].nodes[1]].x)
getfiniteelement(fs::AbstractFEMFunctionSpace) = fs.fe
getnlocaldofs(fs::AbstractFEMFunctionSpace) = getnbasefunctions(fs)
getmesh(fs::AbstractFEMFunctionSpace) = fs.mesh
getncomponents(fs::AbstractFEMFunctionSpace) = fs.components

function shape_value(fs::AbstractFEMFunctionSpace{dim,T}, q_point::Int, base_func::Int) where {dim,T}
    if fs.components == 1 
        fs.N[base_func, q_point]
    else
        N_comp = zeros(T, fs.components)
        n = size(fs.N,1)
        N_comp[div(base_func,n+1)+1] = fs.N[mod1(base_func,n),q_point]
        return N_comp
    end
end

function shape_gradient(fs::AbstractFEMFunctionSpace{dim,T}, q_point::Int, base_func::Int)  where {dim,T}
    if fs.components == 1
        fs.dNdx[base_func, q_point]
    else
        dN_comp = zeros(T, fs.components, fs.components)
        n = size(fs.N,1)
        dN_comp[div(base_func,n+1)+1, :] = fs.dNdξ[mod1(base_func,n), q_point]
        return Tensors.Tensor{2,dim,T}((dN_comp...,)) ⋅ fs.jacobians.Jinv
    end
end

function shape_divergence(fs::AbstractFEMFunctionSpace, q_point::Int, base_func::Int)
    if fs.components == 1
        sum(fs.dNdx[base_func, q_point])
    else
        tr(shape_gradient(fs, q_point, base_func))
    end
end