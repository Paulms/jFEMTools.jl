# Apply Dirichlet boundary condition
struct Dirichlet{T}  #<: Constraint
    prescribed_dofs::Vector{Int}
    values::Vector{T}
end

function Dirichlet(dh::DofHandler{2}, u::TrialFunction, edgeset::String,f::Union{Function,Real})
    element = u.element
    edgeset = getedgeset(dh.mesh, edgeset)
    prescribed_dofs = Vector{Int}()
    values = Vector{Float64}()
    for edge_idx in edgeset
        cell_idx = edge_idx[1]  #cell idx
        cell = getcells(dh.mesh)[cell_idx]
        cell_topology = gettopology(dh.mesh,cell, element)
        edge_lidx = edge_idx[2]  # edge idx
        l_dof = Int[]
        offset::Int = dh.cell_dofs_offset[cell_idx] - 1 + field_offset(dh, u, cell_idx)
        for j = 1:2 #Add vertex dofs
            local_offset = reference_edge_vertices(dh.mesh,cell)[edge_lidx][j]
            if !(dh.cell_dofs[offset+local_offset] ∈ prescribed_dofs)
                push!(prescribed_dofs, dh.cell_dofs[offset+local_offset])
                push!(l_dof, local_offset)
            end
        end
        dofs_per_edge = Int(cell_topology[1]/getnedges(dh.mesh,cell))
        for j in 1:dofs_per_edge  #Add edge dofs
            local_offset::Int = cell_topology[0] + dofs_per_edge*(edge_lidx-1) + j
            if !(dh.cell_dofs[offset+local_offset] ∈ prescribed_dofs)
                push!(prescribed_dofs, dh.cell_dofs[offset+local_offset])
                push!(l_dof, local_offset)
            end
        end
        _push_values!(values, cell_idx, dh.mesh, l_dof, element, f)
    end
    #now put all in order
    p = sortperm(prescribed_dofs)
    return Dirichlet(prescribed_dofs[p], values[p])
end

function _push_values!(values::Vector, cell::Int, mesh,l_dof::Vector{Int}, felem::AbstractElement, f::Function)
    for i in l_dof
        vals = f(spatial_nodal_coordinate(mesh,cell,felem,i))
        push!(values,vals)
    end
end

function _push_values!(values::Vector, cell::Int, mesh,l_dof::Vector{Int}, felem::AbstractElement, vals::Real)
    for i in l_dof
        push!(values,vals)
    end
end

@enum(ApplyStrategy, APPLY_TRANSPOSE, APPLY_INPLACE)

function apply!(KK::Union{SparseMatrixCSC,Symmetric}, f::AbstractVector, dirichlet::Dirichlet;
                strategy::ApplyStrategy=APPLY_TRANSPOSE)
    K = isa(KK, Symmetric) ? KK.data : KK
    @assert length(f) == 0 || length(f) == size(K, 1)
    @boundscheck checkbounds(K, dirichlet.prescribed_dofs, dirichlet.prescribed_dofs)
    @boundscheck length(f) == 0 || checkbounds(f, dirichlet.prescribed_dofs)

    m = meandiag(K) # Use the mean of the diagonal here to not ruin things for iterative solver
    @inbounds for i in 1:length(dirichlet.values)
        d = dirichlet.prescribed_dofs[i]
        v = dirichlet.values[i]

        if v != 0
            for j in nzrange(K, d)
                f[K.rowval[j]] -= v * K.nzval[j]
            end
        end
    end
    zero_out_columns!(K, dirichlet.prescribed_dofs)
    if strategy == APPLY_TRANSPOSE
        K′ = copy(K)
        transpose!(K′, K)
        zero_out_columns!(K′, dirichlet.prescribed_dofs)
        transpose!(K, K′)
    elseif strategy == APPLY_INPLACE
        K[dirichlet.prescribed_dofs, :] = 0
    else
        error("Unknown apply strategy")
    end
    @inbounds for i in 1:length(dirichlet.values)
        d = dirichlet.prescribed_dofs[i]
        v = dirichlet.values[i]
        K[d, d] = m
        if length(f) != 0
            f[d] = v * m
        end
    end
end

# columns need to be stored entries, this is not checked
function zero_out_columns!(K, dofs::Vector{Int}) # can be removed in 0.7 with #24711 merged
    #@debug assert(issorted(dofs))
    for col in dofs
        r = nzrange(K, col)
        K.nzval[r] .= 0.0
    end
end

function meandiag(K::AbstractMatrix)
    z = zero(eltype(K))
    for i in 1:size(K, 1)
        z += abs(K[i, i])
    end
    return z / size(K, 1)
end
