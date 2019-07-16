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
struct VertexIndex
    cellidx::Int
    idx::Int
end

# Vertices
struct Vertex{dim,T}
    x::Tensors.Vec{dim, T}
end

@inline get_coordinates(vertex::Vertex) = vertex.x

#--------------- Cells
struct Cell{dim, N, M, L}
    vertices::NTuple{N,Int}
end

#Common cell types
const TriangleCell = Cell{2,3,3,1}
@inline get_cell_name(::TriangleCell) = "Triangle"
@inline reference_edge_vertices(::Type{TriangleCell}) = ((2,3),(3,1),(1,2))

const RectangleCell = Cell{2,4,4,1}
@inline get_cell_name(::RectangleCell) = "Rectangle"
@inline reference_edge_vertices(::Type{RectangleCell}) = ((1,2),(2,3),(3,4),(4,1))

# API
@inline getnvertices(cell::Cell{dim,N}) where {dim,N} = N
@inline getnedges(cell::Cell{dim,N,M}) where {dim,N,M} = M
@inline getnfaces(cell::Cell{dim,N,M,P}) where {dim,N,M,P} = P

# ----------------- Mesh
struct PolytopeMesh{dim,T,C} <: AbstractPolytopeMesh
    cells::Vector{C}
    vertices::Vector{Vertex{dim,T}}
    # Sets
    cellsets::Dict{String,Set{Int}}
    facesets::Dict{String,Set{FaceIndex}}
    edgesets::Dict{String,Set{EdgeIndex}}
    vertexsets::Dict{String,Set{VertexIndex}}
end

function PolytopeMesh(cells,
              vertices;
              cellsets::Dict{String,Set{Int}}=Dict{String,Set{Int}}(),
              facesets::Dict{String,Set{FaceIndex}}=Dict{String,Set{FaceIndex}}(),
              edgesets::Dict{String,Set{EdgeIndex}}=Dict{String,Set{EdgeIndex}}(),
              vertexsets::Dict{String,Set{VertexIndex}}=Dict{String,Set{VertexIndex}}()) where {dim}
    return PolytopeMesh(cells, vertices, cellsets, facesets, edgesets, vertexsets)
end

# API

@inline getncells(mesh::PolytopeMesh) = length(mesh.cells)
@inline getnvertices(mesh::PolytopeMesh) = length(mesh.vertices)
@inline getverticesidx(mesh, cell_idx) = mesh.cells[cell_idx].vertices
@inline getvertexset(mesh::PolytopeMesh, set::String) = mesh.vertexsets[set]

"""
function getcoords(mesh, vertex::VertexIndex)
Return a Tensor.Vec with the coordinates of `vertex`
"""
function getcoords(mesh::PolytopeMesh, vertex::VertexIndex)
    K = mesh.cells[vertex.cellidx]
    return mesh.vertices[K.vertices[vertex.idx]].x
end

@inline getvertexcoords(mesh::PolytopeMesh, vertex_idx::Int) = mesh.vertices[vertex_idx].x

function mapToGlobalIdx(mesh,vertexset::Set{VertexIndex})
    globalVertexSet = Set{Int}()
    for vertex in vertexset
        K = mesh.cells[vertex.cellidx]
        push!(globalVertexSet,K.vertices[vertex.idx])
    end
    return globalVertexSet
end

"""
    getverticescoords(mesh::PolytopeMesh, cell_idx)
Return a vector with the coordinates of the vertices of cell number `cell`.
"""
@inline function getverticescoords(mesh::PolytopeMesh{dim,T}, cell_idx::Int) where {dim,T}
    N = getnvertices(mesh.cells[cell_idx])
    coords = Vector{Tensors.Vec{dim,T}}(undef, N)
    for (i,j) in enumerate(mesh.cells[cell_idx].vertices)
        coords[i] = mesh.vertices[j].x
    end
    return coords
end

function get_vertices_matrix(mesh::PolytopeMesh{dim,T,C}) where {dim,T,C}
    vertices_m = Matrix{T}(undef,length(mesh.vertices),dim)
    for (k,vertex) in enumerate(mesh.vertices)
        vertices_m[k,:] = vertex.x
    end
    vertices_m
end
function get_conectivity_list(mesh::PolytopeMesh{dim,T,C}) where {dim,T,C}
    cells_m = Vector()
    for k = 1:getncells(mesh)
        push!(cells_m,mesh.cells[k].vertices)
    end
    cells_m
end

function cell_volume(mesh::PolytopeMesh{2}, cell_idx::Int)
    N = getnvertices(mesh.cells[cell_idx])
    verts = getverticescoords(mesh,cell_idx)
    return 0.5*abs(sum(verts[j][1]*verts[mod1(j+1,N)][2]-verts[mod1(j+1,N)][1]*verts[j][2] for j ∈ 1:N))
end

function cell_centroid(mesh::PolytopeMesh{2}, cell_idx::Int)
    verts = getverticescoords(mesh,cell_idx)
    Ve = cell_volume(mesh, cell_idx)
    N = getnvertices(mesh.cells[cell_idx])
    xc = 1/(6*Ve)*sum((verts[j][1]*verts[mod1(j+1,N)][2]-verts[mod1(j+1,N)][1]*verts[j][2])*(verts[j][1]+verts[mod1(j+1,N)][1]) for j ∈ 1:N)
    yc = 1/(6*Ve)*sum((verts[j][1]*verts[mod1(j+1,N)][2]-verts[mod1(j+1,N)][1]*verts[j][2])*(verts[j][2]+verts[mod1(j+1,N)][2]) for j ∈ 1:N)
    return Tensors.Vec{2}((xc,yc))
end

function cell_diameter(mesh::PolytopeMesh{dim,T}, cell_idx::Int) where {dim,T}
    K = mesh.cells[cell_idx]
    verts = getverticescoords(mesh,cell_idx)
    h = zero(T)
     for k in reference_edge_vertices(typeof(K))
        mσ = norm(verts[k[2]] - verts[k[1]])
        h = max(h, mσ)
    end
    h
end
