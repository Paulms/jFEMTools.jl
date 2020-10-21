# Generate 2D vertices in a rectangular Lattice
function _generate_2d_vertices!(vertices::Vector{Tensors.Vec{2,T}}, nx, ny, LL, LR, UR, UL) where {T}
  for i in 0:ny-1
    ratio_bounds = i / (ny-1)

    x0 = LL[1] * (1 - ratio_bounds) + ratio_bounds * UL[1]
    x1 = LR[1] * (1 - ratio_bounds) + ratio_bounds * UR[1]

    y0 = LL[2] * (1 - ratio_bounds) + ratio_bounds * UL[2]
    y1 = LR[2] * (1 - ratio_bounds) + ratio_bounds * UR[2]

    for j in 0:nx-1
        ratio = j / (nx-1)
        x = x0 * (1 - ratio) + ratio * x1
        y = y0 * (1 - ratio) + ratio * y1
        push!(vertices, Tensors.Vec{2}((x, y)))
    end
end
end

# Shortcuts
unitSquareMesh2(cellType, nel::NTuple{2,Int}) = rectangle_mesh2(cellType, nel, Tensors.Vec{2}((0.0,0.0)), Tensors.Vec{2}((1.0,1.0)))

#########################
# Triangle Cells 2D   #
#########################
"""
rectangle_mesh2(::Type{TriangleCell}, nel::NTuple{2,Int}, LL::Vec{2,T}, UR::Vec{2,T})
Generate a rectangular mesh with triangular cells, where `LL` is the low left vertex
and `UR` is the upper right one. `nel` is a tuple with the number of partions to be
used in each dimension.
"""
function rectangle_mesh2(::Type{TriangleCell}, nel::NTuple{2,Int}, LL::Tensors.Vec{2,T}, UR::Tensors.Vec{2,T}) where {T}
  LR = Tensors.Vec{2}((UR[1],LL[2]))
  UL = Tensors.Vec{2}((LL[1],UR[2]))
  nel_x = nel[1]; nel_y = nel[2]; nel_tot = 2*nel_x*nel_y
  n_vertices_x = nel_x + 1; n_vertices_y = nel_y + 1
  n_vertices = n_vertices_x * n_vertices_y

  # Generate vertices
  vertices = Tensors.Vec{2,Float64}[]
  _generate_2d_vertices!(vertices, n_vertices_x, n_vertices_y, LL, LR, UR, UL)
  typeof(vertices)

  # Generate cells
  vertex_array = reshape(collect(1:n_vertices), (n_vertices_x, n_vertices_y))
  offsets = [1]
  k = 1
  indices = []
  for j in 1:nel_y, i in 1:nel_x
      k = k + 3
      push!(offsets,k)
      push!(indices,(vertex_array[i,j], vertex_array[i+1,j], vertex_array[i,j+1])...) # ◺
      k = k + 3
      push!(offsets,k)
      push!(indices, (vertex_array[i+1,j], vertex_array[i+1,j+1], vertex_array[i,j+1])...) # ◹
  end

  geometry = Dict((2,0) => MeshConectivity(Tuple(indices), Tuple(offsets)))
  entities = (n_vertices, nel_x+nel_y+3*nel_x*nel_y, nel_tot)

  # Compute boundaries
  cell_array = reshape(collect(1:nel_tot),(2, nel_x, nel_y))
  boundary = Tuple{Int,Int}[[(cl, 1) for cl in cell_array[1,:,1]];
                               [(cl, 1) for cl in cell_array[2,end,:]];
                               [(cl, 2) for cl in cell_array[2,:,end]];
                               [(cl, 3) for cl in cell_array[1,1,:]]]
  offset = 0
  edgesets = Dict{String,Set{NTuple{2,Int}}}()
  edgesets["bottom"] = Set{Tuple{Int,Int}}(boundary[(1:length(cell_array[1,:,1]))   .+ offset]); offset += length(cell_array[1,:,1])
  edgesets["right"]  = Set{Tuple{Int,Int}}(boundary[(1:length(cell_array[2,end,:])) .+ offset]); offset += length(cell_array[2,end,:])
  edgesets["top"]    = Set{Tuple{Int,Int}}(boundary[(1:length(cell_array[2,:,end])) .+ offset]); offset += length(cell_array[2,:,end])
  edgesets["left"]   = Set{Tuple{Int,Int}}(boundary[(1:length(cell_array[1,1,:]))   .+ offset]); offset += length(cell_array[1,1,:])
  edgesets["boundary"] = union(edgesets["bottom"],edgesets["right"],edgesets["top"],edgesets["left"])

  return PolytopalMesh2(entities,vertices,geometry,Dict(1=>edgesets))
end

#########################
# Rectangle Cells 2D   #
#########################
"""
rectangle_mesh2(::Type{RectangleCell}, nel::NTuple{2,Int}, LL::Vec{2,T}, UR::Vec{2,T}) where {T}
Generate a rectangular mesh with triangular cells, where `LL` is the low left vertex
and `UR` is the upper right one. `nel` is a tuple with the number of partions to be
used in each dimension.
"""
function rectangle_mesh2(::Type{RectangleCell}, nel::NTuple{2,Int}, LL::Tensors.Vec{2,T}, UR::Tensors.Vec{2,T}) where {T}
  LR = Tensors.Vec{2}((UR[1],LL[2]))
  UL = Tensors.Vec{2}((LL[1],UR[2]))
  nel_x = nel[1]; nel_y = nel[2]; nel_tot = nel_x*nel_y
  n_vertices_x = nel_x + 1; n_vertices_y = nel_y + 1
  n_vertices = n_vertices_x * n_vertices_y

  # Generate vertices
  vertices = Tensors.Vec{2,Float64}[]
  _generate_2d_vertices!(vertices, n_vertices_x, n_vertices_y, LL, LR, UR, UL)

  # Generate cells
  vertex_array = reshape(collect(1:n_vertices), (n_vertices_x, n_vertices_y))
  offsets = [1]
  k = 1
  indices = []
  for j in 1:nel_y, i in 1:nel_x
    k = k + 4
    push!(offsets,k)
    push!(indices,(vertex_array[i,j], vertex_array[i+1,j], vertex_array[i+1,j+1], vertex_array[i,j+1])...)
  end

  geometry = Dict((2,0) => MeshConectivity(Tuple(indices), Tuple(offsets)))
  entities = (n_vertices, nel_x+nel_y+2*nel_x*nel_y, nel_tot)

  # Cell faces
  cell_array = reshape(collect(1:nel_tot),(nel_x, nel_y))
  boundary = Tuple{Int,Int}[[(cl, 1) for cl in cell_array[:,1]];
                            [(cl, 2) for cl in cell_array[end,:]];
                            [(cl, 3) for cl in cell_array[:,end]];
                            [(cl, 4) for cl in cell_array[1,:]]]
  # Cell face sets
  offset = 0
  edgesets = Dict{String, Set{Tuple{Int,Int}}}()
  edgesets["bottom"] = Set{Tuple{Int,Int}}(boundary[(1:length(cell_array[:,1]))   .+ offset]); offset += length(cell_array[:,1])
  edgesets["right"]  = Set{Tuple{Int,Int}}(boundary[(1:length(cell_array[end,:])) .+ offset]); offset += length(cell_array[end,:])
  edgesets["top"]    = Set{Tuple{Int,Int}}(boundary[(1:length(cell_array[:,end])) .+ offset]); offset += length(cell_array[:,end])
  edgesets["left"]   = Set{Tuple{Int,Int}}(boundary[(1:length(cell_array[1,:]))   .+ offset]); offset += length(cell_array[1,:])
  edgesets["boundary"] = union(edgesets["bottom"],edgesets["right"],edgesets["top"],edgesets["left"])
  return PolytopalMesh2(entities,vertices,geometry,Dict(1=>edgesets))
end

# #########################
# # Hexagonal Cells 2D
# #############################

function _gen_hexagon2(centroid,LL,UR,hex_width,hex_heigth)
  hex_verts = [
      Tensors.Vec{2}((centroid[1]-hex_width/2,centroid[2]-hex_heigth/4)),
      Tensors.Vec{2}((centroid[1],centroid[2]-hex_heigth/2)),
      Tensors.Vec{2}((centroid[1]+hex_width/2,centroid[2]-hex_heigth/4)),
      Tensors.Vec{2}((centroid[1]+hex_width/2,centroid[2]+hex_heigth/4)),
      Tensors.Vec{2}((centroid[1],centroid[2]+hex_heigth/2)),
      Tensors.Vec{2}((centroid[1]-hex_width/2,centroid[2]+hex_heigth/4)),
  ]
  filter_verts = filter(x->(x[1]>=LL[1] && x[1]<=UR[1]
                        && x[2]>=LL[2] && x[2] <= UR[2]),hex_verts)
  if size(filter_verts,1) == 2 #corners case
      x = centroid[1] + (centroid[1] == LL[1] ? hex_width/2 : -hex_width/2)
      new_vertex = Tensors.Vec{2}((x, centroid[2]))
      if centroid == LL
          filter_verts = (centroid,new_vertex,filter_verts...)
      elseif centroid[2] == LL[2] && centroid[1] == UR[1]
          filter_verts = (new_vertex,centroid,filter_verts...)
      elseif centroid[2] == UR[2] && centroid[1] == LL[1]
          filter_verts = (filter_verts...,new_vertex,centroid)
      elseif centroid == UR
          filter_verts = (filter_verts...,centroid,new_vertex)
      else
          throw("error on hexagon at corner $centroid")
      end
  elseif size(filter_verts,1) == 3 #top and bottom pentagon/triangle case
      nv1 = Tensors.Vec{2}((centroid[1]-hex_width/2, centroid[2]))
      nv2 = Tensors.Vec{2}((centroid[1]+hex_width/2, centroid[2]))
      if centroid[2] == LL[2]
          filter_verts = (nv1,nv2,filter_verts...)
      elseif centroid[2] == UR[2]
          filter_verts = (nv1,filter_verts...,nv2)
      elseif centroid[2] > UR[2]
          filter_verts = (filter_verts...,)
      else
          throw("error at bottom centroid $centroid")
      end
  end
  return filter_verts
end

function rectangle_mesh2(::Type{HexagonCell}, nel::NTuple{2,Int}, LL::Tensors.Vec{2,T}, UR::Tensors.Vec{2,T}) where {T}
  LR = Tensors.Vec{2}((UR[1],LL[2]))
  UL = Tensors.Vec{2}((LL[1],UR[2]))
  nel_x = nel[1]; nel_y = nel[2]

  # Generate vertices
  vertices = Tensors.Vec{2,Float64}[]
  hex_width = 1.0
  hex_heigth = 1.0
  w_scale = (LR[1] - LL[1])/nel_x
  h_scale = (UL[2] - LL[1])/nel_y
  LLn = Tensors.Vec{2}((0.0,0.0))
  URn = Tensors.Vec{2}((nel_x,nel_y))

  map_func(x::Tensors.Vec{2,T}) where T = Tensors.Vec{2}((x[1]*w_scale+LL[1],x[2]*h_scale+LL[2]))

  centroids = Tensors.Vec{2,T}[]
  _generate_2d_hex_centroids!(centroids,LL, nel_y, nel_x, hex_width, hex_heigth)
  used_vertices = Dict{Tensors.Vec{2,T},Int}()
  #Add cells
  offsets = [1]
  k = 1
  indices = []
  nextvert = 1 # next free vertex to use
  for c_i in 1:size(centroids,1)
          cell_verts = Int[]
          for vert in _gen_hexagon2(centroids[c_i],LLn,URn,hex_width,hex_heigth)
              nextvert = _push_cell_vertex!(vert,cell_verts, used_vertices,vertices,nextvert, map_func)
          end
          k = k + size(cell_verts,1)
          push!(offsets,k)
          push!(indices,cell_verts...)
  end
  nel_tot = size(centroids,1)
  n_vertices = size(vertices,1)
  geometry = Dict((2,0) => MeshConectivity(Tuple(indices), Tuple(offsets)))

  # Cell edges
  n_edges = 0
  edgesets = Dict{String,Set{NTuple{2,Int}}}()
  t_edge = mod1(nel_y,3)==1 ? 4 : 3
  t_idxs = (nel_tot-nel_x + (mod1(nel_y,3)==3 ? 0 : 1)):nel_tot
  l_bounds = Tuple{Int,Int}[]
  r_bounds = Tuple{Int,Int}[]
  c_idx = 1
  for i in 1:nel_y
      if mod1(i,3) == 1
          push!(l_bounds,(c_idx,4))
          push!(r_bounds,(c_idx+nel_x,2))
          ed_idx = (i == nel_y ? 5 : 6)
          push!(l_bounds,(c_idx+nel_x+1,ed_idx))
          push!(r_bounds,(c_idx+2*nel_x,3))
          c_idx = c_idx + nel_x*2+1
          n_edges = n_edges + 3 + 6*nel_x
      elseif mod1(i,3) == 2
          push!(l_bounds,(c_idx,4))
          push!(r_bounds,(c_idx+nel_x,2))
          c_idx = c_idx + nel_x + 1
          n_edges = n_edges + 1 + 6*nel_x
      else
          push!(l_bounds,(c_idx,6))
          push!(r_bounds,(c_idx+nel_x-1,3))
          if i == nel_y
              push!(l_bounds,(c_idx+nel_x,4))
              push!(r_bounds,(c_idx+2*nel_x,2))
          end
          c_idx = c_idx + nel_x
          n_edges = n_edges + 3 + 5*nel_x
      end
  end
  edgesets["bottom"] = Set{Tuple{Int,Int}}([(c_i,1) for c_i in 1:(nel_x+1)])
  edgesets["top"]    = Set{Tuple{Int,Int}}([(i,t_edge) for i in t_idxs])
  edgesets["right"]  = Set{Tuple{Int,Int}}([(x[1],x[2]) for x in r_bounds])
  edgesets["left"]   = Set{Tuple{Int,Int}}([(x[1],x[2]) for x in l_bounds])
  edgesets["boundary"] = union(edgesets["bottom"],edgesets["right"],edgesets["top"],edgesets["left"])

  entities = (n_vertices, n_edges, nel_tot)

  return PolytopalMesh2(entities,vertices,geometry,Dict(1=>edgesets))
end

##########################
# 3D Data
##########################
function _build_hexahedron_geometry!(geometry,indices, offsets, nel, nedges,nfaces)
    indices1 = [Int[] for _ in 1:nel]
    indices2 = [Int[] for _ in 1:nfaces]
    indices3 = [Int[] for _ in 1:nel]
    indices4 = [Int[] for _ in 1:nedges]
    facesDict = Dict{Set{Int},Int}()
    edgesDict = Dict{Set{Int},Int}()
    nextidx = 1 #first face index
    next_e_idx = 1 #first face index
    for cell in 1:nel
        nodes = indices[offsets[cell]:offsets[cell+1]-1]
        # Compute faces geometries
        Vi = [nodes[1:4],nodes[5:8],nodes[[2,3,6,7]],nodes[[3,4,7,8]],nodes[[1,4,5,8]],nodes[[1,2,5,6]]]
        for vi in Vi
            token = Base.ht_keyindex2!(facesDict, Set(vi))
            if token > 0 # reuse edge index
            reuse_idx = facesDict.vals[token]
            push!(indices1[cell], reuse_idx)
            else
            push!(indices1[cell], nextidx)
            indices2[nextidx] = vi
            Base._setindex!(facesDict, nextidx, Set(vi), -token)
            nextidx += 1
            end
        end
        # Compute edges geometries
        Vi = [nodes[1:2],nodes[2:3],nodes[3:4],nodes[[1,4]],nodes[5:6],nodes[6:7],
              nodes[7:8],nodes[[5,8]],nodes[[2,6]],nodes[[3,7]],nodes[[4,8]],nodes[[1,5]]]
        for vi in Vi
            token = Base.ht_keyindex2!(edgesDict, Set(vi))
            if token > 0 # reuse edge index
            reuse_idx = edgesDict.vals[token]
            push!(indices3[cell], reuse_idx) 
            else
            push!(indices3[cell], next_e_idx)
            indices4[next_e_idx] = vi
            Base._setindex!(edgesDict, next_e_idx, Set(vi), -token)
            next_e_idx += 1
            end
        end
    end
    V = _pack_connectivity(indices1)
    push!(geometry,(3,2)=>V)
    V = _pack_connectivity(indices2)
    push!(geometry,(2,0)=>V)
    V = _pack_connectivity(indices3)
    push!(geometry,(3,1)=>V)
    V = _pack_connectivity(indices4)
    push!(geometry,(1,0)=>V)
  end

#####################################
# Hexahedron 3D
#####################################
function hyper_rectagle_mesh2(::Type{HexahedronCell}, nel::NTuple{3,Int}, left::Vec{3,T}=Vec{3}((-1.0,-1.0,-1.0)), right::Vec{3,T}=Vec{3}((1.0,1.0,1.0))) where {T}
    nel_x = nel[1]; nel_y = nel[2]; nel_z = nel[3]; nel_tot = nel_x*nel_y*nel_z
    n_nodes_x = nel_x + 1; n_nodes_y = nel_y + 1; n_nodes_z = nel_z + 1
    n_nodes = n_nodes_x * n_nodes_y * n_nodes_z

    # Generate vertices
    coords_x = range(left[1], stop=right[1], length=n_nodes_x)
    coords_y = range(left[2], stop=right[2], length=n_nodes_y)
    coords_z = range(left[3], stop=right[3], length=n_nodes_z)
    nodes = Tensors.Vec{3,T}[]
    for k in 1:n_nodes_z, j in 1:n_nodes_y, i in 1:n_nodes_x
        push!(nodes, Tensors.Vec{3}((coords_x[i], coords_y[j], coords_z[k])))
    end

    # Generate cells
    node_array = reshape(collect(1:n_nodes), (n_nodes_x, n_nodes_y, n_nodes_z))
    offsets = [1]
    l = 1
    indices = Int[]
    for k in 1:nel_z, j in 1:nel_y, i in 1:nel_x
        l = l + 8
        push!(offsets, l)
        push!(indices, (node_array[i,j,k], node_array[i+1,j,k], node_array[i+1,j+1,k], node_array[i,j+1,k],
                                 node_array[i,j,k+1], node_array[i+1,j,k+1], node_array[i+1,j+1,k+1], node_array[i,j+1,k+1])...)
    end

    n_edges = (nel_x+nel_y+2*nel_x*nel_y)*(nel_z+1) + (nel_y+1)*nel_z*(nel_x+1)
    n_faces = nel_x*nel_y*(nel_z+1)+nel_y*nel_z*(nel_x+1)+nel_x*nel_z*(nel_y+1)
    entities = (n_nodes, n_edges, n_faces,nel_tot)

    geometry = Dict((3,0) => MeshConectivity(Tuple(indices), Tuple(offsets)))
    _build_hexahedron_geometry!(geometry, indices, offsets,nel_tot, n_edges, n_faces)

    # Cell faces
    cell_array = reshape(collect(1:nel_tot),(nel_x, nel_y, nel_z))
    boundary = Tuple{Int,Int}[[(cl, 1) for cl in cell_array[:,:,1][:]];
                              [(cl, 2) for cl in cell_array[:,1,:][:]];
                              [(cl, 3) for cl in cell_array[end,:,:][:]];
                              [(cl, 4) for cl in cell_array[:,end,:][:]];
                              [(cl, 5) for cl in cell_array[1,:,:][:]];
                              [(cl, 6) for cl in cell_array[:,:,end][:]]]

    # Cell edge sets
    offset = 0
    edgesets = Dict{String,Set{Tuple{Int,Int}}}()
    edgesets["bottom"] = Set{Tuple{Int,Int}}(boundary[(1:length(cell_array[:,:,1][:]))   .+ offset]); offset += length(cell_array[:,:,1][:])
    edgesets["front"]  = Set{Tuple{Int,Int}}(boundary[(1:length(cell_array[:,1,:][:]))   .+ offset]); offset += length(cell_array[:,1,:][:])
    edgesets["right"]  = Set{Tuple{Int,Int}}(boundary[(1:length(cell_array[end,:,:][:])) .+ offset]); offset += length(cell_array[end,:,:][:])
    edgesets["back"]   = Set{Tuple{Int,Int}}(boundary[(1:length(cell_array[:,end,:][:])) .+ offset]); offset += length(cell_array[:,end,:][:])
    edgesets["left"]   = Set{Tuple{Int,Int}}(boundary[(1:length(cell_array[1,:,:][:]))   .+ offset]); offset += length(cell_array[1,:,:][:])
    edgesets["top"]    = Set{Tuple{Int,Int}}(boundary[(1:length(cell_array[:,:,end][:])) .+ offset]); offset += length(cell_array[:,:,end][:])
    edgesets["boundary"] = union(edgesets["bottom"],edgesets["right"],edgesets["top"],edgesets["left"],edgesets["back"],edgesets["front"])

    return PolytopalMesh2(entities,nodes,geometry,Dict(1=>edgesets))
end

# Tetrahedron
function hyper_rectagle_mesh2(::Type{TetrahedronCell}, nel::NTuple{3,Int}, left::Vec{3,T}=Vec{3}((-1.0,-1.0,-1.0)), right::Vec{3,T}=Vec{3}((1.0,1.0,1.0))) where {T}
    nodes_per_dim = nel .+ 1

    cells_per_cube = 6
    n_nodes = prod(nodes_per_dim)
    nel_tot = cells_per_cube * prod(nel)

    n_nodes_x, n_nodes_y, n_nodes_z = nodes_per_dim
    nel_x, nel_y, nel_z = nel

    # Generate nodes
    coords_x = range(left[1], stop=right[1], length=n_nodes_x)
    coords_y = range(left[2], stop=right[2], length=n_nodes_y)
    coords_z = range(left[3], stop=right[3], length=n_nodes_z)
    numbering = reshape(1:n_nodes, nodes_per_dim)

    # Pre-allocate the nodes & cells
    nodes = Vector{Tensors.Vec{3,T}}(undef, n_nodes)

    # Generate nodes
    node_idx = 1
    @inbounds for k in 1:n_nodes_z, j in 1:n_nodes_y, i in 1:n_nodes_x
        nodes[node_idx] = Tensors.Vec{3}((coords_x[i], coords_y[j], coords_z[k]))
        node_idx += 1
    end

    # Generate cells, case 1 from: http://www.baumanneduard.ch/Splitting%20a%20cube%20in%20tetrahedras2.htm
    # cube = (1, 2, 3, 4, 5, 6, 7, 8)
    # left = (1, 4, 5, 8), right = (2, 3, 6, 7)
    # front = (1, 2, 5, 6), back = (3, 4, 7, 8)
    # bottom = (1, 2, 3, 4), top = (5, 6, 7, 8)
    offsets = [1]
    l = 1
    indices = []
    @inbounds for k in 1:nel_z, j in 1:nel_y, i in 1:nel_x
        cell = (
            numbering[i  , j  , k],
            numbering[i+1, j  , k],
            numbering[i+1, j+1, k],
            numbering[i  , j+1, k],
            numbering[i  , j  , k+1],
            numbering[i+1, j  , k+1],
            numbering[i+1, j+1, k+1],
            numbering[i  , j+1, k+1]
        )
        l += 4; push!(offsets, l)
        push!(indices, (cell[1], cell[2], cell[4], cell[8])...)
        l += 4; push!(offsets, l)
        push!(indices, (cell[1], cell[5], cell[2], cell[8])...)
        l += 4; push!(offsets, l)
        push!(indices, (cell[2], cell[3], cell[4], cell[8])...)
        l += 4; push!(offsets, l)
        push!(indices, (cell[2], cell[7], cell[3], cell[8])...)
        l += 4; push!(offsets, l)
        push!(indices, (cell[2], cell[5], cell[6], cell[8])...)
        l += 4; push!(offsets, l)
        push!(indices, (cell[2], cell[6], cell[7], cell[8])...)
    end

    geometry = Dict((3,0) => MeshConectivity(Tuple(indices), Tuple(offsets)))
    n_edges = 2*(nel_x+nel_y+nel_x*nel_y)*(nel_z+1) + (2*nel_y+1)*nel_z*(nel_x+1) + nel_y*nel_x*nel_z
    n_faces = 2*(nel_x*nel_y*(nel_z+1)+nel_y*nel_z*(nel_x+1)+nel_x*nel_z*(nel_y+1)) + 6*nel_x*nel_y*nel_z
    entities = (n_nodes, n_edges, n_faces,nel_tot)

    # Order the cells as c_nxyz[n, x, y, z] such that we can look up boundary cells
    c_nxyz = reshape(1:nel_tot, (cells_per_cube, nel...))

    @views le = [map(x -> (x,4), c_nxyz[1, 1, :, :][:])   ; map(x -> (x,2), c_nxyz[2, 1, :, :][:])]
    @views ri = [map(x -> (x,1), c_nxyz[4, end, :, :][:]) ; map(x -> (x,1), c_nxyz[6, end, :, :][:])]
    @views fr = [map(x -> (x,1), c_nxyz[2, :, 1, :][:])   ; map(x -> (x,1), c_nxyz[5, :, 1, :][:])]
    @views ba = [map(x -> (x,3), c_nxyz[3, :, end, :][:]) ; map(x -> (x,3), c_nxyz[4, :, end, :][:])]
    @views bo = [map(x -> (x,1), c_nxyz[1, :, :, 1][:])   ; map(x -> (x,1), c_nxyz[3, :, :, 1][:])]
    @views to = [map(x -> (x,3), c_nxyz[5, :, :, end][:]) ; map(x -> (x,3), c_nxyz[6, :, :, end][:])]

    edgesets = Dict(
        "left" => Set(le),
        "right" => Set(ri),
        "front" => Set(fr),
        "back" => Set(ba),
        "bottom" => Set(bo),
        "top" => Set(to),
        "boundary" => Set(union(le,ri,fr,ba,bo,to))
    )

    return PolytopalMesh2(entities,nodes,geometry,Dict(1=>edgesets))

end