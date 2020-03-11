struct PoissonVirtualElement{dim} <: AbstractElement
    degree::Int
end

function PoissonVirtualElement(dim = 2, degree=1)
    PoissonVirtualElement{dim}(degree)
end

get_degree(el::PoissonVirtualElement) = el.degree

struct LocalPoissonVirtualElement
    degree::Int
    Pk_basis::Monomials
    edge_quad::QuadratureRule
end

function LocalPoissonVirtualElement(dim, degree, centroid, diameter)
    Pk_basis = Monomials(dim,degree,centroid, diameter)
    if degree == 1
        points = Vector{Float64}(); weights = Vector{Float64}()
    else
        points, weights = FastGaussQuadrature.gausslobatto(2*degree-1)
        # Shift interval from (-1,1) to (0,1)
        weights *= 0.5
        points = points .+ 1.0; points /= 2.0
    end
    quad = QuadratureRule{1,Segment,Float64}(weights, [Tensors.Tensor{1,1}([x]) for x in points])
    LocalPoissonVirtualElement(degree, Pk_basis, quad)
end

get_degree(el::LocalPoissonVirtualElement) = el.degree

function gettopology(mesh::AbstractPolytopalMesh{2},cell, element::PoissonVirtualElement)
    V = getnvertices(mesh,cell)
    deg = get_degree(element)
    Dict(0=> V, 1=> V*(deg-1), 2=> Int((deg-1)*deg/2))
end

function getnlocaldofs(element::PoissonVirtualElement{2}, cell::Cell{2,V,E,M}) where {V,E,M}
    k = get_degree(element)
    V*k+(k-1)*k/2
end

function spatial_nodal_coordinate(mesh, ci::Int,element::PoissonVirtualElement{2},i::Int) where {V,E,M}
    if i <= getnvertices(mesh,getcells(mesh)[ci]) #Vertices coordinates
        return getvertexcoords(mesh,getcells(mesh)[ci],i)
    else     #Edge points coordinates
        cell = getcells(mesh)[ci]
        nv = getnvertices(mesh,cell)
        cell_topology = gettopology(mesh,cell, element)
        dofs_per_edge = Int(cell_topology[1]/getnedges(mesh,cell))
        ei = div((i - nv)-1,dofs_per_edge)+1
        val = 1/(dofs_per_edge+1)*mod1((i-nv),dofs_per_edge)
        map_unit_to_segment(Tensors.Vec{1}((val,)), getverticescoords(mesh,EdgeIndex(ci,ei)))
    end

end
