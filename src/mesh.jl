import Tensors
abstract type AbstractPolytopeMesh end
abstract type AbstractCell{dim,V,F} end

"""
A `FaceIndex` wraps an (Int, Int) and defines a face by pointing to a (cell, face).
"""
struct FaceIndex
    cellidx::Int
    idx::Int
end

"""
A `EdgeIndex` wraps an (Int, Int) and defines an edge by pointing to a (cell, edge).
"""
struct EdgeIndex
    cellidx::Int
    idx::Int
end

"""
A `VertexIndex` wraps an (Int, Int) and defines an vertex by pointing to a (cell, vertex).
"""
struct NodeIndex
    cellidx::Int
    idx::Int
end

# Nodes
struct Node{dim,T}
    x::Tensors.Vec{dim, T}
end

@inline get_coordinates(node::Node) = node.x

#--------------- Cells
struct Cell{dim, N, M, L}
    nodes::NTuple{N,Int}
end

#Common cell types
const TriangleCell = Cell{2,3,3,1}
@inline get_cell_name(::TriangleCell) = "Triangle"
@inline reference_edge_nodes(::Type{TriangleCell}) = ((2,3),(3,1),(1,2))

const RectangleCell = Cell{2,4,4,1}
@inline get_cell_name(::RectangleCell) = "Rectangle"
@inline reference_edge_nodes(::Type{RectangleCell}) = ((1,2),(2,3),(3,4),(4,1))

# ----------------- Mesh
struct PolytopeMesh{dim,T,C} <: AbstractPolytopeMesh
    cells::Vector{C}
    nodes::Vector{Node{dim,T}}
    # Sets
    cellsets::Dict{String,Set{Int}}
    facesets::Dict{String,Set{FaceIndex}}
    edgesets::Dict{String,Set{EdgeIndex}}
    nodesets::Dict{String,Set{NodeIndex}}
end

function PolytopeMesh(cells,
              nodes;
              cellsets::Dict{String,Set{Int}}=Dict{String,Set{Int}}(),
              facesets::Dict{String,Set{FaceIndex}}=Dict{String,Set{FaceIndex}}(),
              edgesets::Dict{String,Set{EdgeIndex}}=Dict{String,Set{EdgeIndex}}(),
              nodesets::Dict{String,Set{NodeIndex}}=Dict{String,Set{NodeIndex}}()) where {dim}
    return PolytopeMesh(cells, nodes, cellsets, facesets, edgesets, nodesets)
end

# API

@inline getncells(mesh::PolytopeMesh) = length(mesh.cells)

function get_vertices_matrix(mesh::PolytopeMesh{dim,T,C}) where {dim,T,C}
    nodes_m = Matrix{T}(undef,length(mesh.nodes),dim)
    for (k,node) in enumerate(mesh.nodes)
        nodes_m[k,:] = node.x
    end
    nodes_m
end
function get_conectivity_list(mesh::PolytopeMesh{dim,T,C}) where {dim,T,C}
    cells_m = Vector()
    for k = 1:getncells(mesh)
        push!(cells_m,mesh.cells[k].nodes)
    end
    cells_m
end
