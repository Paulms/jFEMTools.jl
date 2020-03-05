# Local indices
"""
A `FaceIndex` wraps an (Int, Int) and defines a face by pointing to a (cell, face).
"""
struct FaceIndex
    cellidx::Int
    idx::Int
end

"""
A `FacetIndex` wraps an (Int, Int) and defines a facet by pointing to a (cell, facet).
"""
struct FacetIndex
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

function get_Normal(mesh, edge_idx::EdgeIndex)
  coords = getverticescoords(mesh, edge_idx)
  v1 = coords[2] - coords[1]
  n1 = Tensors.Vec{2}((v1[2], -v1[1]))
  return n1/norm(n1)
end

function get_vertices_matrix(mesh)
  vertices_m = Matrix{eltype(mesh.vertices[1])}(undef,getnvertices(mesh),getdim(mesh))
  for k in 1:getnvertices(mesh)
      vertices_m[k,:] = getvertexcoords(mesh,k)
  end
  vertices_m
end

function cell_volume2(mesh, cell_idx::Int)
  d = getdim(mesh) 
  if  d == 2
    N = getncellvertices(mesh, cell_idx)
    verts = getverticescoords(mesh, cell_idx)
    return 0.5*abs(sum(verts[j][1]*verts[mod1(j+1,N)][2]-verts[mod1(j+1,N)][1]*verts[j][2] for j ∈ 1:N))
  else
    error("cell volume not available for cells of dim $d")
  end
end