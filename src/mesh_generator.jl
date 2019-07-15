function _generate_2d_nodes!(nodes, nx, ny, LL, LR, UR, UL)
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
            push!(nodes, Node(Tensors.Vec{2}((x, y))))
        end
    end
end

# Check edge orientation consistency
function _check_node_data(nodes, n1,n2,n3)
    a = nodes[n2].x-nodes[n1].x
    b = nodes[n3].x-nodes[n1].x
    if (a[1]*b[2]-a[2]*b[1]) < 0
        #swap vertices 2 and 3
        return (n1,n3,n2)
    end
    return (n1,n2,n3)
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
    n_nodes_x = nel_x + 1; n_nodes_y = nel_y + 1
    n_nodes = n_nodes_x * n_nodes_y

    # Generate nodes
    nodes = Node{2,T}[]
    _generate_2d_nodes!(nodes, n_nodes_x, n_nodes_y, LL, LR, UR, UL)

    # Generate cells
    node_array = reshape(collect(1:n_nodes), (n_nodes_x, n_nodes_y))
    cells = TriangleCell[]
    for j in 1:nel_y, i in 1:nel_x
        push!(cells, TriangleCell((node_array[i,j], node_array[i+1,j], node_array[i,j+1]))) # ◺
        push!(cells, TriangleCell((node_array[i+1,j], node_array[i+1,j+1], node_array[i,j+1]))) # ◹
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

    return PolytopeMesh(cells, nodes; edgesets = edgesets)
end

#########################
# Rectangle Cells 2D   #
#########################
@inline _build_cell(::Type{RectangleCell}, el_nodes, el_faces) = RectangleCell(el_nodes,(el_faces[1],el_faces[2],el_faces[3],el_faces[4]))
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
    n_nodes_x = nel_x + 1; n_nodes_y = nel_y + 1
    n_nodes = n_nodes_x * n_nodes_y

    # Generate nodes
    nodes = Node{2,T}[]
    _generate_2d_nodes!(nodes, n_nodes_x, n_nodes_y, LL, LR, UR, UL)

    # Generate cells
    node_array = reshape(collect(1:n_nodes), (n_nodes_x, n_nodes_y))
    cells = RectangleCell[]
    for j in 1:nel_y, i in 1:nel_x
        push!(cells, RectangleCell((node_array[i,j], node_array[i+1,j], node_array[i+1,j+1], node_array[i,j+1])))
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

    return PolytopeMesh(cells, nodes; edgesets = edgesets)
end
