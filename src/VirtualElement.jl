struct VirtualElement{dim}
    degree::Int
end

function VirtualElement(dim = 2, degree=1)
    VirtualElement{dim}(degree)
end

get_degree(el::VirtualElement) = el.degree

struct LocalVirtualElement
    degree::Int
    Pk_basis::Monomials
    edge_quad::QuadratureRule
end

function LocalVirtualElement(dim, degree, centroid, diameter)
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
    LocalVirtualElement(degree, Pk_basis, quad)
end

get_degree(el::LocalVirtualElement) = el.degree

function gettopology(cell::Cell{2,V,E,M}, element::VirtualElement) where {V,E,M}
    deg = get_degree(element)
    Dict(0=> V, 1=> V*(deg-1), 2=> Int((deg-1)*deg/2))
end

function spatial_nodal_coordinate(mesh, ci::Int,element::VirtualElement{2},i::Int) where {V,E,M}
    if i <= getnvertices(getcells(mesh)[ci]) #Vertices coordinates
        return getvertexcoords(mesh,ci,i)
    else     #Edge points coordinates
        cell = mesh.cells[ci]
        nv = getnvertices(cell)
        cell_topology = gettopology(cell, element)
        dofs_per_edge = Int(cell_topology[1]/getnedges(cell))
        ei = div((i - nv)-1,dofs_per_edge)+1
        val = 1/(dofs_per_edge+1)*mod1((i-nv),dofs_per_edge)
        map_unit_to_segment(Tensors.Vec{1}((val,)), getverticescoords(mesh,EdgeIndex(ci,ei)))
    end

end
