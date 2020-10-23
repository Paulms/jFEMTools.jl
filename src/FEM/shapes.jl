abstract type Shape{dim} end

getdim(::Shape{dim}) where {dim} = dim

#########################
# Triangles   #
#########################
struct Simplex{dim} <: Shape{dim} end
const Triangle = Simplex{2}
@inline get_num_faces(::Type{Triangle}) = 3
@inline get_num_vertices(::Type{Triangle}) = 3

function reference_coordinates(::Type{Triangle})
    [Tensors.Vec{2, Float64}((0.0, 0.0)),
    Tensors.Vec{2, Float64}((1.0, 0.0)),
    Tensors.Vec{2, Float64}((0.0, 1.0))]
end
function reference_edges(::Type{Triangle})
    [[Tensors.Vec{2, Float64}((1.0, 0.0)),Tensors.Vec{2, Float64}((0.0, 1.0))],
     [Tensors.Vec{2, Float64}((0.0, 1.0)),Tensors.Vec{2, Float64}((0.0, 0.0))],
     [Tensors.Vec{2, Float64}((0.0, 0.0)),Tensors.Vec{2, Float64}((1.0, 0.0))]]
end
@inline reference_edge_nodes(::Type{Triangle}) = ((1,2),(2,3),(3,1))

function gettopology(::Type{Triangle})
    return Dict(0=>3,1=>3,2=>1)
end

##########################
# Tetrahedron
######################

const Tetrahedron = Simplex{3}

get_num_faces(::Type{Tetrahedron}) = 3
get_num_vertices(::Type{Tetrahedron}) = 3

function reference_coordinates(::Type{Tetrahedron})
    [Tensors.Vec{3, Float64}((0.0, 0.0,0.0)),
     Tensors.Vec{3, Float64}((1.0, 0.0,0.0)),
     Tensors.Vec{3, Float64}((0.0,1.0, 0.0)),
     Tensors.Vec{3, Float64}((0.0,0.0, 1.0))]
end
function reference_edges(::Type{Tetrahedron})
    coords = reference_coordinates(Tetrahedron)
    [[coords[i[1]], coords[i[2]]] for i in reference_edge_nodes(Tetrahedron)]
end
function reference_faces(::Type{Tetrahedron})
    coords = reference_coordinates(Tetrahedron)
    [[coords[i[1]], coords[i[2]], coords[i[3]]] for i in reference_face_nodes(Tetrahedron)]
end
reference_edge_nodes(::Type{Tetrahedron}) = ((1,2),(2,3),(3,4),(4,1),(1,3),(2,4))
reference_face_nodes(::Type{Tetrahedron}) = ((1,2,3),(2,3,4),(1,3,4),(1,2,4))


function gettopology(::Type{Tetrahedron})
    return Dict(0=>4,1=>6,2=>4,3=>1)
end

"""
get points for a nodal basis of order `order` on a `dim`
    dimensional simplex
"""
function get_nodal_points(::Type{Simplex{dim}}, order) where {dim}
    points = Vector{Tensors.Vec{dim,Float64}}()
    vertices = reference_coordinates(Simplex{dim})
    topology = Dict{Int, Int}()
    append!(points, vertices)
    push!(topology, 0=>length(points))
    [append!(points, _interior_points(verts, order)) for verts in reference_edges(Simplex{dim})]
    push!(topology, 1=>length(points)-topology[0])
    if dim == 2
        append!(points, _interior_points(vertices, order))
        push!(topology, 2=>length(points)-topology[0]-topology[1])
    else
        [append!(points, _interior_points(verts, order)) for verts in reference_faces(Simplex{dim})]
        push!(topology, 2=>length(points)-topology[0]-topology[1])
        append!(points, _interior_points(vertices, order))
        push!(topology, 3=>length(points)-topology[0]-topology[1]-topology[2])
    end
    points, topology
end

function _interior_points(verts, order)
    n = length(verts)
    ls = [(verts[i] - verts[1])/order for i in 2:n]
    m = length(ls)
    grid_indices =  []
    if m == 1
        grid_indices = [[i] for i in 1:order-1]
    elseif m == 2 && order > 2
        grid_indices = [[i,j] for i in 1:order-1 for j in 1:order-i-1]
    elseif m == 3 && order > 3
        grid_indices = [[i,j,k] for i in 1:order-3 for j in 1:order-i-2 for k in 1:order-1-j-i]
    end
    pts = Vector{typeof(verts[1])}()
    for indices in grid_indices
        res = verts[1]
        for (i,ii) in enumerate(indices)
            res += (ii) * ls[m - i+1]
        end
        push!(pts,res)
    end
    pts
end

#########################
# Rectangles  #
#########################
struct HyperCube{dim} <: Shape{dim} end
const Rectangle = HyperCube{2}
@inline get_num_faces(::Type{Rectangle}) = 4
@inline get_num_vertices(::Type{Rectangle}) = 4

function reference_coordinates(::Type{Rectangle})
    [Tensors.Vec{2, Float64}((0.0, 0.0)),
    Tensors.Vec{2, Float64}((1.0, 0.0)),
    Tensors.Vec{2, Float64}((0.0, 1.0)),
    Tensors.Vec{2, Float64}((1.0, 1.0))]
end
function reference_edges(::Type{Rectangle})
    [[Tensors.Vec{2, Float64}((1.0, 0.0)),Tensors.Vec{2, Float64}((1.0, 1.0))],
     [Tensors.Vec{2, Float64}((1.0, 1.0)),Tensors.Vec{2, Float64}((0.0, 1.0))],
     [Tensors.Vec{2, Float64}((0.0, 1.0)),Tensors.Vec{2, Float64}((0.0, 0.0))],
     [Tensors.Vec{2, Float64}((0.0, 0.0)),Tensors.Vec{2, Float64}((1.0, 0.0))]]
end
@inline reference_edge_nodes(::Type{Rectangle}) = ((1,2),(2,3),(3,4),(4,1))

function gettopology(::Type{Rectangle})
    return Dict(0=>4,1=>4,2=>1)
end

##########################
# Hexahedron
######################

const Hexahedron = HyperCube{3}

get_num_faces(::Type{Hexahedron}) = 6
get_num_vertices(::Type{Hexahedron}) = 8

function reference_coordinates(::Type{Hexahedron})
    [Tensors.Vec{3}([0.0, 0.0, 0.0]),
    Tensors.Vec{3}([1.0, 0.0, 0.0]),
    Tensors.Vec{3}([1.0, 1.0, 0.0]),
    Tensors.Vec{3}([0.0, 1.0, 0.0]),
    Tensors.Vec{3}([0.0, 0.0, 1.0]),
    Tensors.Vec{3}([1.0, 0.0, 1.0]),
    Tensors.Vec{3}([1.0, 1.0, 1.0]),
    Tensors.Vec{3}([0.0, 1.0, 1.0])]
end

function reference_edges(::Type{Hexahedron} )
    coords = reference_coordinates(Hexahedron)
    [[coords[1],coords[2]],[coords[2],coords[3]],
    [coords[3],coords[4]],[coords[4],coords[1]],
    [coords[5],coords[6]],[coords[6],coords[7]],
    [coords[7],coords[8]],[coords[8],coords[5]],
    [coords[1],coords[4]],[coords[2],coords[3]],
    [coords[6],coords[7]],[coords[5],coords[8]]]
end

function reference_faces(::Type{Hexahedron})
    coords = reference_coordinates(Hexahedron)
    [[coords[1],coords[2],coords[3],coords[4]],
    [coords[5],coords[6],coords[7],coords[8]],
    [coords[2],coords[3],coords[7],coords[6]],
    [coords[1],coords[5],coords[8],coords[4]],
    [coords[1],coords[2],coords[5],coords[6]],
    [coords[3],coords[4],coords[8],coords[7]]]
end

reference_edge_nodes(::Type{Hexahedron}) = ((1,2),(2,3),(3,4),(4,1),(5,6),(6,7),(7,8),(8,5),(1,4),(2,3),(6,7),(5,8))
reference_face_nodes(::Type{Hexahedron}) = ((1,2,3,4),(6,5,8,7),(2,6,7,3),(1,4,8,5),(3,7,8,4),(1,5,6,2))

function gettopology(::Type{Hexahedron})
    return Dict(0=>8,1=>12,2=>6,3=>1)
end

"""
get points for a nodal basis of order `order` on a `dim`
    dimensional hypercube
"""
function get_nodal_points(::Type{HyperCube{dim}}, order) where {dim}
    points = Vector{Tensors.Vec{dim,Float64}}()
    vertices = reference_coordinates(HyperCube{dim})
    topology = Dict{Int, Int}()
    append!(points, vertices)
    push!(topology, 0=>length(points))
    [append!(points, _interior_points(verts, order)) for verts in reference_edges(HyperCube{dim})]
    push!(topology, 1=>length(points)-topology[0])
    if dim == 2
        append!(points, _cube_interior_points(vertices, order, dim))
        push!(topology, 2=>length(points)-topology[0]-topology[1])
    else
        [append!(points, _cube_interior_points(verts, order, dim)) for verts in reference_faces(HyperCube{dim})]
        push!(topology, 2=>length(points)-topology[0]-topology[1])
        append!(points, _cube_interior_points(vertices, order, dim))
        push!(topology, 3=>length(points)-topology[0]-topology[1]-topology[2])
    end
    points, topology
end

function _cube_interior_points(verts, order, dim)
    n = length(verts)
    ls = [(verts[i] - verts[1])/order for i in 2:n]
    m = length(ls)
    lx = ls[1]
    
    if dim == 2
        for point in ls
            if !any([x ≈ 0 for x in point])
                lx = point
            end
        end
        return [eltype(verts)((verts[1][1]+i*lx[1],verts[1][2]+j*lx[2])) for i in 1:order-1 for j in 1:order-1]
    elseif dim == 3
        if m == 3
            for point in ls
                if sum(point .≈ 0.0) == 1
                    lx = point
                end
            end
            if lx[1] ≈ 0.0
                return [eltype(verts)((verts[1][1],verts[1][2]+j*lx[2], verts[1][3]+k*lx[3])) for j in 1:order-1 for k in 1:order-1]
            elseif lx[2] ≈ 0.0
                return [eltype(verts)((verts[1][1]+i*lx[1],verts[1][2], verts[1][3]+k*lx[3])) for i in 1:order-1 for k in 1:order-1]
            elseif lx[3] ≈ 0.0
                return [eltype(verts)((verts[1][1]+i*lx[1],verts[1][2]+j*lx[2], verts[1][3])) for i in 1:order-1 for j in 1:order-1]
            end
        elseif m == 7 && order > 1
            for point in ls
                if !any([x ≈ 0 for x in point])
                    lx = point
                end
            end
            return [eltype(verts)((verts[1][1]+i*lx[1],verts[1][2]+j*lx[2], verts[1][3]+k*lx[3])) for i in 1:order-1 for j in 1:order-1 for k in 1:order-1]
        end
    else
        error("interior nodes for hypercube of dimension $dim not supported")
    end
    return []
end

############
# Segment
###########
const Segment = Simplex{1}
@inline get_num_faces(::Type{Segment}) = 1
@inline get_num_vertices(::Type{Segment}) = 2

function reference_coordinates(::Type{Segment})
    return [Tensors.Vec{1, Float64}((0.0,)),Tensors.Vec{1, Float64}((1.0,))]
end

function get_nodal_points(::Type{Segment}, order)
    points = Vector{Tensors.Vec{1,Float64}}()
    vertices = reference_coordinates(Segment)
    topology = Dict{Int, Int}()
    append!(points, vertices)
    push!(topology, 0=>length(points))
    append!(points, _interior_points(vertices, order))
    push!(topology, 1=>length(points)-topology[0])
    points, topology
end


### Utils
function map_shape_symbols(symbol)
    if symbol == :Triangle
        Triangle()
    elseif symbol == :Segment
        Segment()
    elseif symbol == :Hexahedron
        Hexahedron()
    elseif symbol == :Quad
        Rectangle()
    elseif symbol == :Tetrahedron
        Tetrahedron()
    else
        error("Shape not available")
    end
end