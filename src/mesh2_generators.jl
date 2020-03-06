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
#unitSquareMesh(cellType, nel::NTuple{2,Int}) = rectangle_mesh(cellType, nel, Tensors.Vec{2}((0.0,0.0)), Tensors.Vec{2}((1.0,1.0)))

#########################
# Triangle Cells 2D   #
#########################
"""
rectangle_mesh(::Type{RefTetrahedron}, ::Type{Val{2}}, nel::NTuple{2,Int}, LL::Vec{2,T}, UR::Vec{2,T})
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
  cells = NTuple[]
  for j in 1:nel_y, i in 1:nel_x
      push!(cells, (vertex_array[i,j], vertex_array[i+1,j], vertex_array[i,j+1])) # ◺
      push!(cells, (vertex_array[i+1,j], vertex_array[i+1,j+1], vertex_array[i,j+1])) # ◹
  end

  geometry = Dict((2,0) => _pack_connectivity(cells))
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
# @inline _build_cell(::Type{RectangleCell}, el_vertices, el_faces) = RectangleCell(el_vertices,(el_faces[1],el_faces[2],el_faces[3],el_faces[4]))
# """
# rectangle_mesh(::Type{RectangleCell}, nel::NTuple{2,Int}, LL::Vec{2,T}, UR::Vec{2,T}) where {T}
# Generate a rectangular mesh with triangular cells, where `LL` is the low left vertex
# and `UR` is the upper right one. `nel` is a tuple with the number of partions to be
# used in each dimension.
# """
# function rectangle_mesh(::Type{RectangleCell}, nel::NTuple{2,Int}, LL::Tensors.Vec{2,T}, UR::Tensors.Vec{2,T}) where {T}
# LR = Tensors.Vec{2}((UR[1],LL[2]))
# UL = Tensors.Vec{2}((LL[1],UR[2]))
# nel_x = nel[1]; nel_y = nel[2]; nel_tot = nel_x*nel_y
# n_vertices_x = nel_x + 1; n_vertices_y = nel_y + 1
# n_vertices = n_vertices_x * n_vertices_y

# # Generate vertices
# vertices = Vertex{2,T}[]
# _generate_2d_vertices!(vertices, n_vertices_x, n_vertices_y, LL, LR, UR, UL)

# # Generate cells
# vertex_array = reshape(collect(1:n_vertices), (n_vertices_x, n_vertices_y))
# cells = RectangleCell[]
# for j in 1:nel_y, i in 1:nel_x
#     push!(cells, RectangleCell((vertex_array[i,j], vertex_array[i+1,j], vertex_array[i+1,j+1], vertex_array[i,j+1])))
# end

# # Cell faces
# cell_array = reshape(collect(1:nel_tot),(nel_x, nel_y))
# boundary = EdgeIndex[[EdgeIndex(cl, 1) for cl in cell_array[:,1]];
#                   [EdgeIndex(cl, 2) for cl in cell_array[end,:]];
#                   [EdgeIndex(cl, 3) for cl in cell_array[:,end]];
#                   [EdgeIndex(cl, 4) for cl in cell_array[1,:]]]
# # Cell face sets
# offset = 0
# edgesets = Dict{String, Set{EdgeIndex}}()
# edgesets["bottom"] = Set{EdgeIndex}(boundary[(1:length(cell_array[:,1]))   .+ offset]); offset += length(cell_array[:,1])
# edgesets["right"]  = Set{EdgeIndex}(boundary[(1:length(cell_array[end,:])) .+ offset]); offset += length(cell_array[end,:])
# edgesets["top"]    = Set{EdgeIndex}(boundary[(1:length(cell_array[:,end])) .+ offset]); offset += length(cell_array[:,end])
# edgesets["left"]   = Set{EdgeIndex}(boundary[(1:length(cell_array[1,:]))   .+ offset]); offset += length(cell_array[1,:])
# edgesets["boundary"] = union(edgesets["bottom"],edgesets["right"],edgesets["top"],edgesets["left"])
# vertexsets = Dict{String,Set{Int}}()
# for set in edgesets
#     vertexsets[set.first] = _get_vertexset_from_edges(cells,set.second)
# end
# return PolytopalMesh(cells, vertices; edgesets = edgesets, vertexsets = vertexsets)
# end

# #########################
# # Hexagonal Cells 2D
# #############################

# function _generate_2d_hex_centroids(T,LL, n_centroid_rows, n_centroid_cols, hex_width, hex_heigth)
# centroids = Tensors.Vec{2,T}[]
# x_coord = LL[1]; y_coord = LL[2]; sign = 1;
# for j in 0:(n_centroid_rows-1)
#     for i in 0:(sign > 0 ? n_centroid_cols : n_centroid_cols-1)
#         centroid = Tensors.Vec{2,T}((x_coord,y_coord))
#         push!(centroids, centroid)
#         x_coord = x_coord + hex_width
#     end
#     y_coord = y_coord + hex_heigth*3/4
#     x_coord = LL[1] + max(0,sign*hex_width/2)
#     sign =  -sign;
# end
# return centroids
# end

# function _gen_hexagon(centroid,LL,UR,hex_width,hex_heigth)
# hex_verts = [
#     Vertex(Tensors.Vec{2}((centroid[1]-hex_width/2,centroid[2]-hex_heigth/4))),
#     Vertex(Tensors.Vec{2}((centroid[1],centroid[2]-hex_heigth/2))),
#     Vertex(Tensors.Vec{2}((centroid[1]+hex_width/2,centroid[2]-hex_heigth/4))),
#     Vertex(Tensors.Vec{2}((centroid[1]+hex_width/2,centroid[2]+hex_heigth/4))),
#     Vertex(Tensors.Vec{2}((centroid[1],centroid[2]+hex_heigth/2))),
#     Vertex(Tensors.Vec{2}((centroid[1]-hex_width/2,centroid[2]+hex_heigth/4))),
# ]
# filter_verts = filter(x->(x.x[1]>=LL[1] && x.x[1]<=UR[1]
#                       && x.x[2]>=LL[2] && x.x[2] <= UR[2]),hex_verts)
# if size(filter_verts,1) == 2 #corners case
#     x = centroid[1] + (centroid[1] == LL[1] ? hex_width/2 : -hex_width/2)
#     new_vertex = Vertex(Tensors.Vec{2}((x, centroid[2])))
#     if centroid == LL
#         filter_verts = (Vertex(centroid),new_vertex,filter_verts...)
#     elseif centroid[2] == LL[2] && centroid[1] == UR[1]
#         filter_verts = (new_vertex,Vertex(centroid),filter_verts...)
#     elseif centroid[2] == UR[2] && centroid[1] == LL[1]
#         filter_verts = (filter_verts...,new_vertex,Vertex(centroid))
#     elseif centroid == UR
#         filter_verts = (filter_verts...,Vertex(centroid),new_vertex)
#     else
#         throw("error on hexagon at corner $centroid")
#     end
# elseif size(filter_verts,1) == 3 #top and bottom pentagon case
#     nv1 = Vertex(Tensors.Vec{2}((centroid[1]-hex_width/2, centroid[2])))
#     nv2 = Vertex(Tensors.Vec{2}((centroid[1]+hex_width/2, centroid[2])))
#     if centroid[2] == LL[2]
#         filter_verts = (nv1,nv2,filter_verts...)
#     elseif centroid[2] == UR[2]
#         filter_verts = (nv1,filter_verts...,nv2)
#     else
#         throw("error at bottom centroid $centroid")
#     end
# end
# return filter_verts
# end

# # Use map_func to avoid float precision differences between the same vertex
# function _push_cell_vertex!(vert,cell_verts, used_vertices,vertices,nextvert, map_func)
#     token = ht_keyindex2!(used_vertices, vert)
#     if token > 0 # reuse dofs
#         reuse_vert = used_vertices.vals[token]
#         push!(cell_verts, reuse_vert)
#     else # token <= 0, use new vertex
#         Base._setindex!(used_vertices, nextvert, vert, -token)
#         push!(cell_verts, nextvert)
#         push!(vertices, map_func(vert))
#         nextvert += 1
#     end
#     return nextvert
# end

# function rectangle_mesh(::Type{HexagonCell}, nel::NTuple{2,Int}, LL::Tensors.Vec{2,T}, UR::Tensors.Vec{2,T}) where {T}
# LR = Tensors.Vec{2}((UR[1],LL[2]))
# UL = Tensors.Vec{2}((LL[1],UR[2]))
# nel_x = nel[1]; nel_y = isodd(nel[2]) ? nel[2] : nel[2]+1
# nel_tot = 2*nel_x*nel_y +nel_y - nel_x

# # Generate vertices
# vertices = Vertex{2,T}[]
# hex_width = 1.0
# hex_heigth = 1.0
# w_scale = (LR[1] - LL[1])/nel_x
# h_scale = (UL[2] - LL[1])/nel_y
# LLn = Tensors.Vec{2}((0.0,0.0))
# URn = Tensors.Vec{2}((nel_x,nel_y))

# map_func(x::Vertex) = Vertex(Tensors.Vec{2}((x.x[1]*w_scale+LL[1],x.x[2]*h_scale+LL[2])))

# n_centroid_cols = nel_x; n_centroid_rows = nel_y+2
# centroids = _generate_2d_hex_centroids(T,LL, n_centroid_rows, n_centroid_cols, hex_width, hex_heigth)
# cells = Cell[]
# used_vertices = Dict{Vertex,Int}()
# #Add cells
# nextvert = 1 # next free vertex to use
# c_i = 1
# sign = 1
# for j in 0:(n_centroid_rows-1)
#     for i in 0:(sign > 0 ? n_centroid_cols : n_centroid_cols-1)
#         cell_verts = Int[]
#         for vert in _gen_hexagon(centroids[c_i],LLn,URn,hex_width,hex_heigth)
#             nextvert = _push_cell_vertex!(vert,cell_verts, used_vertices,vertices,nextvert, map_func)
#         end
#         n = size(cell_verts,1)
#         push!(cells,Cell{2,n,n,1}(Tuple(cell_verts)))
#         c_i +=1
#     end
#     sign =  -sign;
# end

# # Cell edges
# edgesets = Dict{String,Set{EdgeIndex}}()
# ncells = size(cells,1)
# l_idxs = zip(repeat([6,4],div(nel_y,2)+1),cumsum(repeat([nel_x+1,nel_x],2)).+1)
# ft_idx = sum(repeat([nel_x+1,nel_x],2)).+1
# t_idxs = zip([3,repeat([4],nel_x-1)...,3],ft_idx:ft_idx+nel_x)
# r_idxs = zip(repeat([3,2],div(nel_y,2)+1),cumsum(repeat([nel_x,nel_x+1],2)).+(nel_x+1))
# edgesets["bottom"] = Set{EdgeIndex}([EdgeIndex(c_i,1) for c_i in 1:(nel_x+1)])
# edgesets["right"]  = Set{EdgeIndex}([EdgeIndex(1+nel_x,2),[EdgeIndex(j,i) for (i,j) in r_idxs]...])
# edgesets["top"]    = Set{EdgeIndex}([EdgeIndex(j,i) for (i,j) in t_idxs])
# edgesets["left"]   = Set{EdgeIndex}([EdgeIndex(1,4),[EdgeIndex(j,i) for (i,j) in l_idxs]...])
# edgesets["boundary"] = union(edgesets["bottom"],edgesets["right"],edgesets["top"],edgesets["left"])

# vertexsets = Dict{String,Set{Int}}()
# for set in edgesets
#     vertexsets[set.first] = _get_vertexset_from_edges(cells,set.second)
# end
# return PolytopalMesh(cells, vertices; edgesets = edgesets, vertexsets = vertexsets)
# end
