function _generate_2d_hex_centroids(T,LL, n_centroid_rows, n_centroid_cols, hex_width, hex_heigth)
    centroids = Tensors.Vec{2,T}[]
    x_coord = LL[1]; y_coord = LL[2]; sign = 1;
    for j in 0:(n_centroid_rows-1)
        for i in 0:(sign > 0 ? n_centroid_cols : n_centroid_cols-1)
            centroid = Tensors.Vec{2,T}((x_coord,y_coord))
            push!(centroids, centroid)
            x_coord = x_coord + hex_width
        end
        y_coord = y_coord + hex_heigth
        x_coord = LL[1] + max(0,sign*hex_width/2)
        sign =  -sign;
    end
    return centroids
end

function _gen_up_pentagon(centroid,hex_width,hex_heigth)
    return [
        Vertex(Tensors.Vec{2}((centroid[1]-hex_width/2,centroid[2]))),
        Vertex(Tensors.Vec{2}((centroid[1]+hex_width/2,centroid[2]))),
        Vertex(Tensors.Vec{2}((centroid[1]+hex_width/2,centroid[2]+hex_heigth/4))),
        Vertex(Tensors.Vec{2}((centroid[1],centroid[2]+hex_heigth/2))),
        Vertex(Tensors.Vec{2}((centroid[1]-hex_width/2,centroid[2]+hex_heigth/4))),
    ]
end

function _gen_down_pentagon(centroid,hex_width,hex_heigth)
    return [
        Vertex(Tensors.Vec{2}((centroid[1]-hex_width/2,centroid[2]))),
        Vertex(Tensors.Vec{2}((centroid[1]-hex_width/2,centroid[2]-hex_heigth/4))),
        Vertex(Tensors.Vec{2}((centroid[1],centroid[2]-hex_heigth/2))),
        Vertex(Tensors.Vec{2}((centroid[1]+hex_width/2,centroid[2]-hex_heigth/4))),
        Vertex(Tensors.Vec{2}((centroid[1]+hex_width/2,centroid[2]))),
    ]
end

function _gen_hexagon(centroid,LL,UR,hex_width,hex_heigth)
    hex_verts = [
        Vertex(Tensors.Vec{2}((centroid[1]-hex_width/2,centroid[2]-hex_heigth/4))),
        Vertex(Tensors.Vec{2}((centroid[1],centroid[2]-hex_heigth/2))),
        Vertex(Tensors.Vec{2}((centroid[1]+hex_width/2,centroid[2]-hex_heigth/4))),
        Vertex(Tensors.Vec{2}((centroid[1]+hex_width/2,centroid[2]+hex_heigth/4))),
        Vertex(Tensors.Vec{2}((centroid[1],centroid[2]+hex_heigth/2))),
        Vertex(Tensors.Vec{2}((centroid[1]-hex_width/2,centroid[2]+hex_heigth/4))),
    ]
    filter_verts = filter(x->(x.x[1]>=LL[1] && x.x[1]<=UR[1]),hex_verts)
    return size(filter_verts,1) < 4 ? (Vertex(centroid), filter_verts...) : filter_verts
end

function _push_cell_vertex!(vert,cell_verts, used_vertices,vertices,nextvert)
        token = ht_keyindex2!(used_vertices, vert)
        if token > 0 # reuse dofs
            reuse_vert = used_vertices.vals[token]
            push!(cell_verts, reuse_vert)
        else # token <= 0, use new vertex
            Base._setindex!(used_vertices, nextvert, vert, -token)
            push!(cell_verts, nextvert)
            push!(vertices, vert)
            nextvert += 1
        end
        return nextvert
end

function rectangle_mesh(::Type{HexagonCell}, nel::NTuple{2,Int}, LL::Tensors.Vec{2,T}, UR::Tensors.Vec{2,T}) where {T}
    LR = Tensors.Vec{2}((UR[1],LL[2]))
    UL = Tensors.Vec{2}((LL[1],UR[2]))
    nel_x = nel[1]; nel_y = isodd(nel[2]) ? nel[2] : nel[2]+1
    nel_tot = 2*nel_x*nel_y +nel_y - nel_x

    # Generate vertices
    vertices = Vertex{2,T}[]
    hex_width = (LR[1] - LL[1])/nel_x
    hex_heigth = (UL[2] - LL[1])/nel_y

    n_centroid_cols = nel_x; n_centroid_rows = nel_y+2
    centroids = _generate_2d_hex_centroids(T,LL, n_centroid_rows, n_centroid_cols, hex_width, hex_heigth)

    cells = Cell[]
    used_vertices = Dict{Vertex,Int}()
    #Add first row

    nextvert = 1 # next free vertex to use
    cell_verts = [
        Vertex(centroids[1]),
        Vertex(Tensors.Vec{2}((centroids[1][1]+hex_width/2,centroids[1][2]))),
        Vertex(Tensors.Vec{2}((centroids[1][1]+hex_width/2,centroids[1][2]+hex_heigth/4))),
        Vertex(Tensors.Vec{2}((centroids[1][1],centroids[1][2]+hex_heigth/2)))
    ]
    cell_verts_idx = Int[]
    for vert in cell_verts
        nextvert = _push_cell_vertex!(vert,cell_verts_idx, used_vertices,vertices,nextvert)
    end
    push!(cells,Cell{2,4,4,1}(Tuple(cell_verts_idx)))
    c_i = 2
    for i in 1:(n_centroid_cols-1)
        cell_verts = Int[]
        for vert in _gen_up_pentagon(centroids[c_i],hex_width,hex_heigth)
            nextvert = _push_cell_vertex!(vert,cell_verts, used_vertices,vertices,nextvert)
        end
        push!(cells,Cell{2,5,5,1}(Tuple(cell_verts)))
        c_i +=1
    end
    cell_verts = [
        Vertex(Tensors.Vec{2}((centroids[c_i][1]-hex_width/2,centroids[c_i][2]))),
        Vertex(centroids[c_i]),
        Vertex(Tensors.Vec{2}((centroids[c_i][1],centroids[c_i][2]+hex_heigth/2))),
        Vertex(Tensors.Vec{2}((centroids[c_i][1]-hex_width/2,centroids[c_i][2]+hex_heigth/4)))
    ]
    cell_verts_idx = Int[]
    for vert in cell_verts
        nextvert = _push_cell_vertex!(vert,cell_verts_idx, used_vertices,vertices,nextvert)
    end
    push!(cells,Cell{2,4,4,1}(Tuple(cell_verts_idx)))
    c_i+1

    # Add middle rows
    sign = -1
    for j in 1:(n_centroid_rows-2)
        for i in 0:(sign > 0 ? n_centroid_cols : n_centroid_cols-1)
            cell_verts = Int[]
            for vert in _gen_hexagon(centroids[c_i],LL,UR,hex_width,hex_heigth)
                nextvert = _push_cell_vertex!(vert,cell_verts, used_vertices,vertices,nextvert)
            end
            n = size(cell_verts,1)
            push!(cells,Cell{2,n,n,1}(Tuple(cell_verts)))
            c_i +=1
        end
        sign =  -sign;
    end

    #Add final row
    cell_verts = [
        Vertex(Tensors.Vec{2}((centroids[c_i][1],centroids[c_i][2]-hex_heigth/2))),
        Vertex(Tensors.Vec{2}((centroids[c_i][1]+hex_width/2,centroids[c_i][2]-hex_heigth/4))),
        Vertex(Tensors.Vec{2}((centroids[c_i][1]+hex_width/2,centroids[c_i][2]))),
        Vertex(centroids[c_i])
    ]
    cell_verts_idx = Int[]
    for vert in cell_verts
        nextvert = _push_cell_vertex!(vert,cell_verts_idx, used_vertices,vertices,nextvert)
    end
    push!(cells,Cell{2,4,4,1}(Tuple(cell_verts_idx)))
    c_i +=1
    for i in 1:(n_centroid_cols-1)
        cell_verts = Int[]
        for vert in _gen_down_pentagon(centroids[c_i],hex_width,hex_heigth)
            nextvert = _push_cell_vertex!(vert,cell_verts, used_vertices,vertices,nextvert)
        end
        push!(cells,Cell{2,5,5,1}(Tuple(cell_verts)))
        c_i +=1
    end
    cell_verts = [
        Vertex(Tensors.Vec{2}((centroids[c_i][1]-hex_width/2,centroids[c_i][2]-hex_heigth/4))),
        Vertex(Tensors.Vec{2}((centroids[c_i][1],centroids[c_i][2]-hex_heigth/2))),
        Vertex(centroids[c_i]),
        Vertex(Tensors.Vec{2}((centroids[c_i][1]-hex_width/2,centroids[c_i][2])))
    ]
    cell_verts_idx = Int[]
    for vert in cell_verts
        nextvert = _push_cell_vertex!(vert,cell_verts_idx, used_vertices,vertices,nextvert)
    end
    push!(cells,Cell{2,4,4,1}(Tuple(cell_verts_idx)))


    # Cell edges

    # vertexsets = Dict{String,Set{Int}}()
    # for set in edgesets
    #     vertexsets[set.first] = _get_vertexset_from_edges(cells,set.second, TriangleCell)
    # end
    return PolytopalMesh(cells, vertices)#; edgesets = edgesets, vertexsets = vertexsets)
end

import Tensors
import Base.ht_keyindex2!
include("../src/mesh.jl")
nel = (3,3)
LL = Tensors.Vec{2}((0.0,0.0))
UR = Tensors.Vec{2}((6.0,6.0))

mesh = rectangle_mesh(HexagonCell,(3,3),LL,UR)
