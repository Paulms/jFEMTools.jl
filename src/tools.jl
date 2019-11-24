map_unit_to_segment(xi::Tensors.Vec{1,T}, x::Vector{Tensors.Vec{2,T}}) where {T} = x[1]*(1-xi[]) + x[2]*xi[]

δ(i, j) = i == j

function polygon_area(verts)
    N = size(verts,1)
    return 0.5*abs(sum(verts[j][1]*verts[mod1(j+1,N)][2]-verts[mod1(j+1,N)][1]*verts[j][2] for j ∈ 1:N))
end

""" Compute volume of a simplex spanned by vertices `verts` """
function simplex_area(verts::Vector{Tensors.Vec{N, T}}) where {N,T}
    # Volume of reference simplex element is 1/n!
    ref_verts = reference_coordinates(RefSimplex, Val{N})
    A, b = get_affine_map(ref_verts, verts)
    F = svd(A)
    return *(F.S...)/factorial(N)
end

""" Map points points from reference simplex to simplex K """
function mapPointsFromReference(::Type{RefSimplex}, K::Vector{Tensors.Vec{N, T}}, points::Vector{Tensors.Vec{N, T}}) where {N,T}
    ref_verts = reference_coordinates(RefSimplex, Val{N})
    A, b = get_affine_map(ref_verts, K)
    return [Tensors.Vec{N}(A*x+Tensors.Vec{N}(b)) for x in points]
end

"""
get_affine_map(x, ξ)
Get (A,b) such that ξ = A * x + b is the affine
mapping from the simplex with vertices x to the simplex of vertices ξ.
"""
# function get_affine_map(x, ξ)
#     @assert length(x[1]) == length(ξ[1]) "dimension mismatch"
#     n = length(x[1])
#     T = eltype(x[1])
#
#     B = zeros(T, n*(n+1), n*(n+1))
#     L = zeros(T, n*(n+1))
#
#     for i in 1:length(x)
#         for j in 1:n
#             row = (i-1) * n + j
#             B[row, n*(j-1)+1:n*j] = x[i]
#             L[row] = ξ[i][j]
#             B[row, n * n + j] = one(T)
#         end
#     end
#     X = B\L
#     return reshape(X[1:n*n], (n, n)), X[n*n+1:end]
# end
function get_affine_map(x, ξ)
    @assert length(x[1]) == length(ξ[1]) "dimension mismatch"
    n = length(x[1])
    T = eltype(x[1])
    A = zeros(T,n,n)
    b = zeros(T,n)
    K = zeros(n+1,n+1)
    c = zeros(n+1)
    for i in 1:n
        K[:,n+1] .= 1
        for j in 1:size(x,1)
            K[j,1:n] .= x[j]
            c[j] = ξ[j][i]
        end
        a = K\c
        A[i,:] .= a[1:n]
        b[i] = a[n+1]
    end
    return A,b
end

""" Tesselate convex 2D shape with `vertices` using Delaunay """
function get_2Dtesselation(vertices)
    #Map points to min_coord <= x <= max_coord
    maxc = 0.0; minc = 0.0
    for x in vertices
        lmin,lmax = minmax(x[1],x[2])
        minc = min(lmin,minc); maxc=max(maxc,lmax)
    end
    width = VoronoiDelaunay.max_coord - VoronoiDelaunay.min_coord
    normalize(x) = (x-minc)/maxc*width+VoronoiDelaunay.min_coord

    tess = VoronoiDelaunay.DelaunayTessellation(2)
    push!(tess, VoronoiDelaunay.Point2D[VoronoiDelaunay.Point(normalize(x[1]),normalize(x[2])) for x in vertices])
    tesselation = Vector{Vector{Tensors.Vec{2,Float64}}}()
    for trig in tess
        mtrig = Vector{Tensors.Vec{2,Float64}}()
        for edge in (trig._a, trig._b, trig._c)
            x = (edge._x - VoronoiDelaunay.min_coord)*maxc*width + minc
            y = (edge._y - VoronoiDelaunay.min_coord)*maxc*width + minc
            push!(mtrig, Tensors.Vec{2}((x,y)))
        end
        push!(tesselation, mtrig)
    end
    return tesselation
end
