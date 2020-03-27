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

"""
get points for a nodal basis of order `order` on a `dim`
    dimensional simplex
"""
function get_nodal_points(::Type{Triangle}, order)
    points = Vector{Tensors.Vec{2,Float64}}()
    vertices = reference_coordinates(Triangle)
    topology = Dict{Int, Int}()
    append!(points, vertices)
    push!(topology, 0=>length(points))
    [append!(points, _interior_points(verts, order)) for verts in reference_edges(Triangle)]
    push!(topology, 1=>length(points)-topology[0])
    append!(points, _interior_points(vertices, order))
    push!(topology, 2=>length(points)-topology[0]-topology[1])
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
    else
        error("Shape not available")
    end
end