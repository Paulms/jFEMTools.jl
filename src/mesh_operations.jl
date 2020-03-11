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

function Base.getindex(e::EdgeIndex, n::Int64)
  if n == 1
    return e.cellidx
  elseif n ==2
    return e.idx
  else
    throw(BoundsError(e,n))
  end
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

function cell_volume(mesh, cell_idx::Int)
  d = getdim(mesh) 
  if  d == 2
    N = getncellvertices(mesh, cell_idx)
    verts = getverticescoords(mesh, cell_idx)
    return 0.5*abs(sum(verts[j][1]*verts[mod1(j+1,N)][2]-verts[mod1(j+1,N)][1]*verts[j][2] for j ∈ 1:N))
  else
    error("cell volume not available for cells of dim $d")
  end
end

function cell_centroid(mesh, cell_idx::Int)
  d = getdim(mesh) 
  if  d == 2
    verts = getverticescoords(mesh, cell_idx)
    vertices = [StaticArrays.SVector(x[1],x[2]) for x in verts]
    chull = PlanarConvexHulls.ConvexHull{PlanarConvexHulls.CCW}(vertices)
    return Tensors.Vec{2}(PlanarConvexHulls.centroid(chull))
  else
    error("cell centroid not available for cells of dim $d")
  end
end

function cell_diameter(mesh, cell_idx::Int)
  d = getdim(mesh) 
  if  d == 2
    verts = getverticescoords(mesh, cell_idx)
    vertices = [StaticArrays.SVector(x[1],x[2]) for x in verts]
    chull = PlanarConvexHulls.ConvexHull{PlanarConvexHulls.CCW}(vertices)
    h = 0.0
    for vert in chull.vertices
        h = max(h,maximum([norm(vert-vert2) for vert2 in chull.vertices]))
    end
    h
  else
    error("cell diameter not available for cells of dim $d")
  end
end