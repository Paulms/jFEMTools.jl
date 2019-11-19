abstract type AbstractPolytopalMesh end
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

@inline reference_edge_vertices(::Type{Cell{2, N}})  where {N} = ((i,mod1(i+1,N)) for i in 1:N)

# API
@inline getnvertices(cell::Cell{dim,N}) where {dim,N} = N
@inline getnedges(cell::Cell{dim,N,M}) where {dim,N,M} = M
@inline getnfaces(cell::Cell{dim,N,M,P}) where {dim,N,M,P} = P
gettopology(cell::Cell{1,N,M,P}) where {N,M,P} = Dict(0=>N,1=>M)
gettopology(cell::Cell{2,N,M,P}) where {N,M,P} = Dict(0=>N,1=>M,2=>P)
gettopology(cell::Cell{3,N,M,P}) where {N,M,P} = Dict(0=>N,1=>M,2=>P,3=>1)

# ----------------- Mesh
struct PolytopalMesh{dim,T,C} <: AbstractPolytopalMesh
    cells::Vector{C}
    vertices::Vector{Vertex{dim,T}}
    # Sets
    cellsets::Dict{String,Set{Int}}
    facesets::Dict{String,Set{FaceIndex}}
    edgesets::Dict{String,Set{EdgeIndex}}
    vertexsets::Dict{String,Set{Int}}
end

function PolytopalMesh(cells,
              vertices;
              cellsets::Dict{String,Set{Int}}=Dict{String,Set{Int}}(),
              facesets::Dict{String,Set{FaceIndex}}=Dict{String,Set{FaceIndex}}(),
              edgesets::Dict{String,Set{EdgeIndex}}=Dict{String,Set{EdgeIndex}}(),
              vertexsets::Dict{String,Set{Int}}=Dict{String,Set{Int}}()) where {dim}
    return PolytopalMesh(cells, vertices, cellsets, facesets, edgesets, vertexsets)
end

# API

@inline getncells(mesh::PolytopalMesh) = length(mesh.cells)
@inline getnvertices(mesh::PolytopalMesh) = length(mesh.vertices)
@inline getverticesidx(mesh, cell_idx) = mesh.cells[cell_idx].vertices
@inline getvertexset(mesh::PolytopalMesh, set::String) = mesh.vertexsets[set]
@inline getedgeset(mesh::PolytopalMesh, set::String) = mesh.edgesets[set]
getcells(mesh::PolytopalMesh) = mesh.cells
function topology_elements(mesh,cellidx,element::Int)
    if element == 0
        return mesh.cells[cellidx].vertices
    elseif element == 1
        return Tuple(EdgeIndex(cellidx,i) for i in 1:getnedges(mesh.cells[cellidx]))
    else
        throw("Topology element of order $element not available for cell type")
    end
end

"""
function getcoords(mesh, vertex_idx::Int)
Return a Tensor.Vec with the coordinates of vertex with index `vertex_idx`
"""
@inline getvertexcoords(mesh::PolytopalMesh, vertex_idx::Int) = mesh.vertices[vertex_idx].x

"""
    getverticescoords(mesh::PolytopalMesh, cell_idx)
Return a vector with the coordinates of the vertices of cell number `cell`.
"""
function getverticescoords(mesh::PolytopalMesh{dim,T}, cell_idx::Int) where {dim,T}
    N = getnvertices(mesh.cells[cell_idx])
    coords = Vector{Tensors.Vec{dim,T}}(undef, N)
    for (i,j) in enumerate(mesh.cells[cell_idx].vertices)
        coords[i] = mesh.vertices[j].x
    end
    return coords
end

function getvertexcoords(mesh::PolytopalMesh{dim,T}, cell_idx::Int, vidx::Int) where {dim,T}
    return mesh.vertices[mesh.cells[cell_idx].vertices[vidx]].x
end

function getverticescoords(mesh::PolytopalMesh{dim,T}, edge_idx::EdgeIndex) where {dim,T}
    cell = getcells(mesh)[edge_idx.cellidx]
    ref_edge = reference_edge_vertices(typeof(cell))[edge_idx.idx]
    coords = Vector{Tensors.Vec{dim,T}}(undef, 2)
    for i in 1:2
        coords[i] = mesh.vertices[cell.vertices[ref_edge[i]]].x
    end
    return coords
end

function getverticesindices(mesh::PolytopalMesh{dim,T}, edge_idx::EdgeIndex) where {dim,T}
    cell = getcells(mesh)[edge_idx.cellidx]
    ref_edge = reference_edge_vertices(typeof(cell))[edge_idx.idx]
    return [cell.vertices[ref_edge[i]] for i in 1:2]
end

function get_Normal(mesh::PolytopalMesh{dim,T}, edge_idx::EdgeIndex) where {dim,T}
    coords = getverticescoords(mesh, edge_idx)
    v1 =  coords[2] - coords[1]
    n1 = Tensors.Vec{2}((v1[2], -v1[1]))
    return n1/norm(n1)
end

function get_vertices_matrix(mesh::PolytopalMesh{dim,T,C}) where {dim,T,C}
    vertices_m = Matrix{T}(undef,length(mesh.vertices),dim)
    for (k,vertex) in enumerate(mesh.vertices)
        vertices_m[k,:] = vertex.x
    end
    vertices_m
end
function get_conectivity_list(mesh::PolytopalMesh{dim,T,C}) where {dim,T,C}
    cells_m = Vector()
    for k = 1:getncells(mesh)
        push!(cells_m,mesh.cells[k].vertices)
    end
    cells_m
end

function cell_volume(mesh::PolytopalMesh{2}, cell_idx::Int)
    N = getnvertices(mesh.cells[cell_idx])
    verts = getverticescoords(mesh,cell_idx)
    return 0.5*abs(sum(verts[j][1]*verts[mod1(j+1,N)][2]-verts[mod1(j+1,N)][1]*verts[j][2] for j ∈ 1:N))
end

function cell_centroid(mesh::PolytopalMesh{2}, cell_idx::Int)
    verts = getverticescoords(mesh,cell_idx)
    vertices = [StaticArrays.SVector(x[1],x[2]) for x in verts]
    chull = PlanarConvexHulls.ConvexHull{PlanarConvexHulls.CCW}(vertices)
    return Tensors.Vec{2}(PlanarConvexHulls.centroid(chull))
end

function cell_diameter(mesh::PolytopalMesh{dim,T}, cell_idx::Int) where {dim,T}
    verts = getverticescoords(mesh,cell_idx)
    vertices = [StaticArrays.SVector(x[1],x[2]) for x in verts]
    chull = PlanarConvexHulls.ConvexHull{PlanarConvexHulls.CCW}(vertices)
    h = 0.0
    for vert in chull.vertices
        h = max(h,maximum([norm(vert-vert2) for vert2 in chull.vertices]))
    end
    h
end
