# Generate 2D vertices in a rectangular Lattice
function _generate_2d_vertices!(vertices, nx, ny, LL, LR, UR, UL)
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
          push!(vertices, Vertex(Tensors.Vec{2}((x, y))))
      end
  end
  end
  
  # Check edge orientation consistency
  function _check_vertex_data(vertices, n1,n2,n3)
  a = vertices[n2].x-vertices[n1].x
  b = vertices[n3].x-vertices[n1].x
  if (a[1]*b[2]-a[2]*b[1]) < 0
      #swap vertices 2 and 3
      return (n1,n3,n2)
  end
  return (n1,n2,n3)
  end

@inline _mapToGlobalIdx(cells,cellidx,localvertexidx) = cells[cellidx].vertices[localvertexidx]

function _get_vertexset_from_edges(cells,edgeset)
  vertices = Set{Int}()
  for edge in edgeset
      push!(vertices, _mapToGlobalIdx(cells, edge.cellidx, reference_edge_vertices(cells[edge.cellidx])[edge.idx][1]))
      push!(vertices, _mapToGlobalIdx(cells, edge.cellidx, reference_edge_vertices(cells[edge.cellidx])[edge.idx][2]))
  end
  return vertices
end

# Shortcuts
unitSquareMesh(cellType, nel::NTuple{2,Int}) = rectangle_mesh(cellType, nel, Tensors.Vec{2}((0.0,0.0)), Tensors.Vec{2}((1.0,1.0)))

#########################
# Triangle Cells 2D   #
#########################
"""
rectangle_mesh(::Type{RefTetrahedron}, ::Type{Val{2}}, nel::NTuple{2,Int}, LL::Vec{2,T}, UR::Vec{2,T})
Generate a rectangular mesh with triangular cells, where `LL` is the low left vertex
and `UR` is the upper right one. `nel` is a tuple with the number of partions to be
used in each dimension.
"""
function rectangle_mesh(::Type{TriangleCell}, nel::NTuple{2,Int}, LL::Tensors.Vec{2,T}, UR::Tensors.Vec{2,T}) where {T}
    LR = Tensors.Vec{2}((UR[1],LL[2]))
    UL = Tensors.Vec{2}((LL[1],UR[2]))
    nel_x = nel[1]; nel_y = nel[2]; nel_tot = 2*nel_x*nel_y
    n_vertices_x = nel_x + 1; n_vertices_y = nel_y + 1
    n_vertices = n_vertices_x * n_vertices_y

    # Generate vertices
    vertices = Vertex{2,T}[]
    _generate_2d_vertices!(vertices, n_vertices_x, n_vertices_y, LL, LR, UR, UL)

    # Generate cells
    vertex_array = reshape(collect(1:n_vertices), (n_vertices_x, n_vertices_y))
    cells = TriangleCell[]
    for j in 1:nel_y, i in 1:nel_x
        push!(cells, TriangleCell((vertex_array[i,j], vertex_array[i+1,j], vertex_array[i,j+1]))) # ◺
        push!(cells, TriangleCell((vertex_array[i+1,j], vertex_array[i+1,j+1], vertex_array[i,j+1]))) # ◹
    end

    # Cell edges
    cell_array = reshape(collect(1:nel_tot),(2, nel_x, nel_y))
    boundary = EdgeIndex[[EdgeIndex(cl, 3) for cl in cell_array[1,:,1]];
                           [EdgeIndex(cl, 3) for cl in cell_array[2,end,:]];
                           [EdgeIndex(cl, 1) for cl in cell_array[2,:,end]];
                           [EdgeIndex(cl, 2) for cl in cell_array[1,1,:]]]
    offset = 0
    edgesets = Dict{String,Set{EdgeIndex}}()
    edgesets["bottom"] = Set{EdgeIndex}(boundary[(1:length(cell_array[1,:,1]))   .+ offset]); offset += length(cell_array[1,:,1])
    edgesets["right"]  = Set{EdgeIndex}(boundary[(1:length(cell_array[2,end,:])) .+ offset]); offset += length(cell_array[2,end,:])
    edgesets["top"]    = Set{EdgeIndex}(boundary[(1:length(cell_array[2,:,end])) .+ offset]); offset += length(cell_array[2,:,end])
    edgesets["left"]   = Set{EdgeIndex}(boundary[(1:length(cell_array[1,1,:]))   .+ offset]); offset += length(cell_array[1,1,:])
    edgesets["boundary"] = union(edgesets["bottom"],edgesets["right"],edgesets["top"],edgesets["left"])
    vertexsets = Dict{String,Set{Int}}()
    for set in edgesets
        vertexsets[set.first] = _get_vertexset_from_edges(cells,set.second)
    end
    return PolytopalMesh(cells, vertices; edgesets = edgesets, vertexsets = vertexsets)
end

#########################
# Rectangle Cells 2D   #
#########################
@inline _build_cell(::Type{RectangleCell}, el_vertices, el_faces) = RectangleCell(el_vertices,(el_faces[1],el_faces[2],el_faces[3],el_faces[4]))
"""
rectangle_mesh(::Type{RectangleCell}, nel::NTuple{2,Int}, LL::Vec{2,T}, UR::Vec{2,T}) where {T}
Generate a rectangular mesh with triangular cells, where `LL` is the low left vertex
and `UR` is the upper right one. `nel` is a tuple with the number of partions to be
used in each dimension.
"""
function rectangle_mesh(::Type{RectangleCell}, nel::NTuple{2,Int}, LL::Tensors.Vec{2,T}, UR::Tensors.Vec{2,T}) where {T}
    LR = Tensors.Vec{2}((UR[1],LL[2]))
    UL = Tensors.Vec{2}((LL[1],UR[2]))
    nel_x = nel[1]; nel_y = nel[2]; nel_tot = nel_x*nel_y
    n_vertices_x = nel_x + 1; n_vertices_y = nel_y + 1
    n_vertices = n_vertices_x * n_vertices_y

    # Generate vertices
    vertices = Vertex{2,T}[]
    _generate_2d_vertices!(vertices, n_vertices_x, n_vertices_y, LL, LR, UR, UL)

    # Generate cells
    vertex_array = reshape(collect(1:n_vertices), (n_vertices_x, n_vertices_y))
    cells = RectangleCell[]
    for j in 1:nel_y, i in 1:nel_x
        push!(cells, RectangleCell((vertex_array[i,j], vertex_array[i+1,j], vertex_array[i+1,j+1], vertex_array[i,j+1])))
    end

    # Cell faces
    cell_array = reshape(collect(1:nel_tot),(nel_x, nel_y))
    boundary = EdgeIndex[[EdgeIndex(cl, 1) for cl in cell_array[:,1]];
                      [EdgeIndex(cl, 2) for cl in cell_array[end,:]];
                      [EdgeIndex(cl, 3) for cl in cell_array[:,end]];
                      [EdgeIndex(cl, 4) for cl in cell_array[1,:]]]
    # Cell face sets
    offset = 0
    edgesets = Dict{String, Set{EdgeIndex}}()
    edgesets["bottom"] = Set{EdgeIndex}(boundary[(1:length(cell_array[:,1]))   .+ offset]); offset += length(cell_array[:,1])
    edgesets["right"]  = Set{EdgeIndex}(boundary[(1:length(cell_array[end,:])) .+ offset]); offset += length(cell_array[end,:])
    edgesets["top"]    = Set{EdgeIndex}(boundary[(1:length(cell_array[:,end])) .+ offset]); offset += length(cell_array[:,end])
    edgesets["left"]   = Set{EdgeIndex}(boundary[(1:length(cell_array[1,:]))   .+ offset]); offset += length(cell_array[1,:])
    edgesets["boundary"] = union(edgesets["bottom"],edgesets["right"],edgesets["top"],edgesets["left"])
    vertexsets = Dict{String,Set{Int}}()
    for set in edgesets
        vertexsets[set.first] = _get_vertexset_from_edges(cells,set.second)
    end
    return PolytopalMesh(cells, vertices; edgesets = edgesets, vertexsets = vertexsets)
end

#########################
# Hexagonal Cells 2D
#############################

function _gen_hexagon(centroid,LL,UR,hex_width,hex_heigth)
    hex_verts = [
        Vertex(Tensors.Vec{2}((centroid[1]-hex_width/2,centroid[2]-hex_heigth/4))),
        Vertex(Tensors.Vec{2}((centroid[1],centroid[2]-hex_heigth/2))),
        Vertex(Tensors.Vec{2}((centroid[1]+hex_width/2,centroid[2]-hex_heigth/4))),
        Vertex(Tensors.Vec{2}((centroid[1]+hex_width/2,centroid[2]+hex_heigth/4))),
        Vertex(Tensors.Vec{2}((centroid[1],centroid[2]+hex_heigth/2))),
        Vertex(Tensors.Vec{2}((centroid[1]-hex_width/2,centroid[2]+hex_heigth/4))),
    ]
    filter_verts = filter(x->(x.x[1]>=LL[1] && x.x[1]<=UR[1]
                          && x.x[2]>=LL[2] && x.x[2] <= UR[2]),hex_verts)
    if size(filter_verts,1) == 2 #corners case
        x = centroid[1] + (centroid[1] == LL[1] ? hex_width/2 : -hex_width/2)
        new_vertex = Vertex(Tensors.Vec{2}((x, centroid[2])))
        if centroid == LL
            filter_verts = (Vertex(centroid),new_vertex,filter_verts...)
        elseif centroid[2] == LL[2] && centroid[1] == UR[1]
            filter_verts = (new_vertex,Vertex(centroid),filter_verts...)
        elseif centroid[2] == UR[2] && centroid[1] == LL[1]
            filter_verts = (filter_verts...,new_vertex,Vertex(centroid))
        elseif centroid == UR
            filter_verts = (filter_verts...,Vertex(centroid),new_vertex)
        else
            throw("error on hexagon at corner $centroid")
        end
    elseif size(filter_verts,1) == 3 #top and bottom pentagon/triangle case
        nv1 = Vertex(Tensors.Vec{2}((centroid[1]-hex_width/2, centroid[2])))
        nv2 = Vertex(Tensors.Vec{2}((centroid[1]+hex_width/2, centroid[2])))
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

function rectangle_mesh(::Type{HexagonCell}, nel::NTuple{2,Int}, LL::Tensors.Vec{2,T}, UR::Tensors.Vec{2,T}) where {T}
    LR = Tensors.Vec{2}((UR[1],LL[2]))
    UL = Tensors.Vec{2}((LL[1],UR[2]))
    nel_x = nel[1]; nel_y = nel[2]
    nel_tot = 2*nel_x*nel_y +nel_y - nel_x

    # Generate vertices
    vertices = Vertex{2,T}[]
    hex_width = 1.0
    hex_heigth = 1.0
    w_scale = (LR[1] - LL[1])/nel_x
    h_scale = (UL[2] - LL[1])/nel_y
    LLn = Tensors.Vec{2}((0.0,0.0))
    URn = Tensors.Vec{2}((nel_x,nel_y))

    map_func(x::Vertex) = Vertex(Tensors.Vec{2}((x.x[1]*w_scale+LL[1],x.x[2]*h_scale+LL[2])))

    centroids = Tensors.Vec{2,T}[]
    _generate_2d_hex_centroids!(centroids,LL, nel_y, nel_x, hex_width, hex_heigth)
    cells = Cell[]
    used_vertices = Dict{Vertex,Int}()
    #Add cells
    nextvert = 1 # next free vertex to use
    for c_i in 1:size(centroids,1)
        cell_verts = Int[]
        for vert in _gen_hexagon(centroids[c_i],LLn,URn,hex_width,hex_heigth)
            nextvert = _push_cell_vertex!(vert,cell_verts, used_vertices,vertices,nextvert, map_func)
        end
        n = size(cell_verts,1)
        push!(cells,Cell{2,n,n,1}(Tuple(cell_verts)))
    end

    # Cell edges
    edgesets = Dict{String,Set{EdgeIndex}}()
    ncells = size(cells,1) 
    t_edge = mod1(nel_y,3)==1 ? 4 : 3
    t_idxs = (ncells-nel_x + (mod1(nel_y,3)==3 ? 0 : 1)):ncells
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
        elseif mod1(i,3) == 2
            push!(l_bounds,(c_idx,4))
            push!(r_bounds,(c_idx+nel_x,2))
            c_idx = c_idx + nel_x + 1
        else
            push!(l_bounds,(c_idx,6))
            push!(r_bounds,(c_idx+nel_x-1,3))
            if i == nel_y
                push!(l_bounds,(c_idx+nel_x,4))
                push!(r_bounds,(c_idx+2*nel_x,2))
            end
            c_idx = c_idx + nel_x
        end
    end
    edgesets["bottom"] = Set{EdgeIndex}([EdgeIndex(c_i,1) for c_i in 1:(nel_x+1)])
    edgesets["top"]    = Set{EdgeIndex}([EdgeIndex(i,t_edge) for i in t_idxs])
    edgesets["right"]  = Set{EdgeIndex}([EdgeIndex(x[1],x[2]) for x in r_bounds])
    edgesets["left"]   = Set{EdgeIndex}([EdgeIndex(x[1],x[2]) for x in l_bounds])
    edgesets["boundary"] = union(edgesets["bottom"],edgesets["right"],edgesets["top"],edgesets["left"])

    vertexsets = Dict{String,Set{Int}}()
    for set in edgesets
        vertexsets[set.first] = _get_vertexset_from_edges(cells,set.second)
    end
    return PolytopalMesh(cells, vertices; edgesets = edgesets, vertexsets = vertexsets)
end
