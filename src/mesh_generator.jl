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

function _get_vertexset_from_edges(cells,edgeset,CellType)
    vertices = Set{Int}()
    for edge in edgeset
        push!(vertices, _mapToGlobalIdx(cells, edge.cellidx, reference_edge_vertices(CellType)[edge.idx][1]))
        push!(vertices, _mapToGlobalIdx(cells, edge.cellidx, reference_edge_vertices(CellType)[edge.idx][2]))
    end
    return vertices
end

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
        vertexsets[set.first] = _get_vertexset_from_edges(cells,set.second, TriangleCell)
    end
    return PolytopeMesh(cells, vertices; edgesets = edgesets, vertexsets = vertexsets)
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
        vertexsets[set.first] = _get_vertexset_from_edges(cells,set.second, RectangleCell)
    end
    return PolytopeMesh(cells, vertices; edgesets = edgesets, vertexsets = vertexsets)
end
