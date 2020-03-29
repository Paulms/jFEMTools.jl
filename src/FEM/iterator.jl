# this file defines iterators used for looping over a grid
mutable struct CellData
    current_cellid::Int
    cell
end

struct CellIterator{dim,T}
    mesh::AbstractPolytopalMesh{dim,T}
    celldata::CellData
    function CellIterator{dim,T}(mesh::AbstractPolytopalMesh{dim,T}) where {dim,T}
        cellid = 0
        cell = nothing
        return new{dim,T}(mesh, CellData(cellid,cell))
    end
end

CellIterator(mesh::AbstractPolytopalMesh{dim,T}) where {dim,T} =
    CellIterator{dim,T}(mesh)

# iterator interface
function Base.iterate(ci::CellIterator, state = 1)
    if state > getncells(ci.mesh)
        return nothing
    else
        return (reinit!(ci, state), state+1)
    end
end
Base.length(ci::CellIterator)  = getncells(ci.mesh)

Base.IteratorSize(::Type{T})   where {T<:CellIterator} = Base.HasLength() # this is default in Base
Base.IteratorEltype(::Type{T}) where {T<:CellIterator} = Base.HasEltype() # this is default in Base
Base.eltype(::Type{T})         where {T<:CellIterator} = T

# utility
@inline cellid(ci::CellIterator) = ci.celldata.current_cellid

function reinit!(ci::CellIterator{dim}, i::Int) where {dim}
    ci.celldata.current_cellid = i
    ci.celldata.cell = getcell(ci.mesh,i)
    return ci
end

@inline reinit!(fs::AbstractFEMFunctionSpace{dim,T,FE}, ci::CellIterator{dim,T}) where {dim,T,FE} = reinit!(fs, getverticescoords(ci.mesh, ci.celldata.cell))
function_value(f::Function, fs::AbstractFEMFunctionSpace{dim,T,FE}, ci::CellIterator,q_point::Int) where {dim,T,FE} = function_value(f,fs, ci.celldata.current_cellid,q_point)
celldofs!(v::Vector, dh::DofHandler, ci::CellIterator) = celldofs!(v, dh, ci.celldata.current_cellid)