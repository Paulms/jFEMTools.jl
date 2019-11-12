struct VirtualElement
    degree::Int
    basis::Monomials
    edge_quad::QuadratureRule
end

function VirtualElement(dim = 2, degree=1)
    basis = Monomials(dim,degree)
    points, weights = FastGaussQuadrature.gausslobatto(degree + 1)
    # Shift interval from (-1,1) to (0,1)
    weights *= 0.5
    points = points .+ 1.0; points /= 2.0
    quad = QuadratureRule{1,Float64}(weights, [Tensors.Tensor{1,1}([x]) for x in points])
    VirtualElement(degree, basis, quad)
end

struct LocalVirtualElement
    basis::Monomials
end

function LocalVirtualElement(cell::Cell{2}, element)
    degree = get_degree(element)
    basis = Monomials(dim, degree)
end

get_degree(el::VirtualElement) = el.degree

function gettopology(cell::Cell{2,V,E,M}, element::VirtualElement) where {V,E,M}
    deg = get_degree(element)
    Dict(0=> V, 1=> V*(deg-1), 2=> Int((deg-1)*deg/2))
end
