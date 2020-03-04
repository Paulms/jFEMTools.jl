struct MeshEntity{dim}
  index::Int 
end 
const MeshVertex = MeshEntity{0}
const MeshEdge = MeshEntity{1}
const MeshFace = MeshEntity{2}

function get_entity_name(entity,dim)
  if entity == 0
    return "vertex"
  elseif entity > 0 && entity < dim
    if entity == 1
      return "edge"
    elseif entity == 2
      return "face"
    end
  elseif entity == dim
    return "cell"
  else
    return "dim $entity"
  end
end

#M = dim + 1
struct MeshConectivity
  indices::NTuple
  offsets::NTuple
end

#dim1 = dim+1 (better solutions?)
struct PolytopalMesh2{dim,dim1,T} #<: AbstractPolytopalMesh
  entities::NTuple{dim1, Int}
  vertices::Vector{NTuple{dim, T}}
  geometry::Dict{NTuple{2,Int},MeshConectivity}
  # Sets
  entitysets::Dict{Int,Dict{String,Set{Int}}}
end

#Common API
getncells(mesh::PolytopalMesh2) = mesh.entities[end]
getnvertices(mesh::PolytopalMesh2) = mesh.entities[1]
getentityset(mesh::PolytopalMesh2, entity::Int, set::String) = mesh.entitysets[entity][set]
getvertexset(mesh::PolytopalMesh2, set::String) = mesh.entitysets[0][set]
getedgeset(mesh::PolytopalMesh2, set::String) = mesh.entitysets[1][set]
function getcells(mesh::PolytopalMesh2{dim}) where {dim}
  return [MeshEntity{dim}(idx) for idx in 1:getncells(mesh)]
end

@inline getnentities(mesh::PolytopalMesh2,d) = mesh.entities[d+1]
function getentities(mesh::PolytopalMesh2, d)
  return [MeshEntity{d}(idx) for idx in 1:getnentities(d)]
end

function topology_elements(mesh::PolytopalMesh2{dim},cellidx,element::Int) where {dim}
  connectivity = get_connectivity!(mesh,dim,element)
  return [MeshEntity{element}(idx) for idx in get_entity_indices(mesh,cellidx)]
end

getvertexcoords(mesh::PolytopalMesh2, vertex_idx::Int) = mesh.vertices[vertex_idx]

function getverticescoords(mesh::PolytopalMesh2, entity::MeshEntity)
  return [mesh.vertices[i] for i in getverticesidx(mesh, entity)]
end

function getvertexcoords(mesh::PolytopalMesh2, entity::MeshEntity, vidx::Int)
  return mesh.vertices[getverticesidx(mesh, entity)[vidx]]
end

#Internal
@inline get_entity_indices(connectivity, idx) = connectivity.indices[connectivity.offsets[idx]:connectivity.offsets[idx+1]-1]
function getverticesidx(mesh::PolytopalMesh2{dim}, entity::MeshEntity{d}) where {dim,d}
  connectivity = get_connectivity!(mesh,d,0)
  return get_entity_indices(connectivity,entity.index)
end

function getentities(mesh::PolytopalMesh2, d2, d1, entity)
  connectivity = get_connectivity!(mesh,d2,d1)
  return [MeshEntity{d1}(i) for i in get_entity_indices(connectivity,entity.index)]
end

function _pack_connectivity(indices_mat)
  offsets = [1]
  k = 1
  indices = []
  for x in indices_mat
    k = k + size(x,1)
    push!(offsets,k)
    push!(indices,x...)
  end
  MeshConectivity(Tuple(indices), Tuple(offsets))
end

# Compute d1 -> d2 from d2 -> d1
function _transpose(mesh, d2, d1)
  indices = [Int[] for _ in 1:getnentities(mesh,d1)]
  for entity_j in getentities(mesh,d2)
    for entity_i in getentities(mesh,d2,d1,entity_j)
      push!(indices[entity_i.index], entity_j.index)
    end
  end
  _pack_connectivity(indices)
end

# Compute d1 -> d2 from d1 -> d3 and d3 -> d2
function _intersection(mesh,d1,d2,d3)
  indices = [Int[] for _ in 1:getnentities(mesh,d1)]
  for entity_i in getentities(mesh,d1)
    for entity_k in getentities(mesh,d1,d3,entity_i)
      for entity_j in getentities(mesh,d3,d2,entity_k)
        if (d1 == d2 && entity_i.index != entity_j.index) || 
          (d1 > d2 && all(x in getentities(mesh,d1,0,entity_i) for x in getentities(mesh,d2,0,entity_j)))
          if !(entity_j.index in indices[entity_i.index]) 
            push!(indices[entity_i.index], entity_j.index)
          end
        end
      end
    end
  end
  _pack_connectivity(indices)
end

function get_connectivity!(mesh::PolytopalMesh2{D},d1::Int,d2::Int) where {D}
  if haskey(mesh.geometry,(d1,d2))
    get(mesh.geometry,(d1,d2),nothing)
  else
    if d1 < d2
      V = get_connectivity!(mesh,d2,d1)
      V = _transpose(mesh,d2,d1)
    else
      d3 = (d1 == 0 && d2 == 0) ? D : 0
      V = get_connectivity!(mesh,d1,d3)
      V = get_connectivity!(mesh,d3,d2)
      V = _intersection(mesh,d1,d2,d3)
    end
    push!(mesh.geometry,(d1,d2)=>V)
    return V
  end
end

function local_vertices(mesh, cell::MeshEntity, d)
  if haskey(mesh.geometry,(d,0))
    V = getverticesidx(mesh, cell)
  elseif d==1

  else
  end
end

function build_d(mesh, d)
  k = 0
  for cell in getcells(mesh)
    local_vertices(mesh,cell,d)
  end
end

################ Auxiliary functions
function get_Normal(mesh::PolytopalMesh2{dim,T}, edge_idx::EdgeIndex) where {dim,T}
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



geometry = Dict((2,0) => MeshConectivity((1,2,4,2,3,4),(1,4,7)))
mesh1 = PolytopalMesh2((4,5,2),[(0.0,0.0),(1.,0.),(0.,1.),(1.,1.)],geometry,Dict{Int,Dict{String,Set{Int}}}())

import Tensors
using jFEMTools
mesh2 = rectangle_mesh(TriangleCell, (1,1), Tensors.Vec{2}((0.0,0.0)), Tensors.Vec{2}((1.0,1.0)))

using BenchmarkTools
@btime PolytopalMesh2((4,5,2),[(0.0,0.0),(1.,0.),(0.,1.),(1.,1.)],geometry,Dict{Int,Dict{String,Set{Int}}}())
@btime rectangle_mesh(TriangleCell, (1,1), Tensors.Vec{2}((0.0,0.0)), Tensors.Vec{2}((1.0,1.0)))