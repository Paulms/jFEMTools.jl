# Patrick Keast, MODERATE-DEGREE TETRAHEDRAL QUADRATURE FORMULAS
# http://mech.fsv.cvut.cz/oofem/resources/doc/oofemrefman/gaussintegrationrule_8C_source.html

function (::Type{QuadratureRule{Tetrahedron}})(quad_type::Gauss, order::Int)
    if order == 0 || order == 1
        points = [Tensors.Vec{3}([1.0/4.0, 1.0/4.0, 1.0/1.4])]
        weigths = [1.0/6.0]
    elseif order == 2
        a = ( 5. + 3. * √(5.) ) / 20.
        b = ( 5. - √(5.) ) / 20.
        w = 1. / 24.
        points = [Tensors.Vec{3}([a,b,b]),
                  Tensors.Vec{3}([b,a,b]),
                  Tensors.Vec{3}([b,b,a]),
                  Tensors.Vec{3}([b,b,b])]
        weigths = [w,w,w,w]
    elseif order == 3
        a1 = 1. / 4.
        a2 = 1. / 2.
        b2 = 1. / 6.
        w1 = -2. / 15.
        w2 = 3. / 40.
        points = [Tensors.Vec{3}([a1,a1,a1]),
                  Tensors.Vec{3}([a2,b2,b2]),
                  Tensors.Vec{3}([b2,a2,b2]),
                  Tensors.Vec{3}([b2,b2,a2]),
                  Tensors.Vec{3}([b2,b2,b2])]
        weigths = [w1,w2,w2,w2,w2]
    elseif order == 4
        a1 = 1. / 4.;
        w1 = -74. / 5625.;

        a2 = 5. / 70.;
        b2 = 11. / 14.;
        w2 = 343. / 45000.;

        a3 = ( 1. + √(5. / 14.) ) / 4.;
        b3 = ( 1. - √(5. / 14.) ) / 4.;
        w3 = 28. / 1125.;
        points = [Tensors.Vec{3}([a1,a1,a1]),
                  Tensors.Vec{3}([b2,a2,a2]),
                  Tensors.Vec{3}([a2,b2,a2]),
                  Tensors.Vec{3}([a2,a2,b2]),
                  Tensors.Vec{3}([a2,a2,a2]),
                  Tensors.Vec{3}([a3,a3,b3]),
                  Tensors.Vec{3}([a3,b3,a3]),
                  Tensors.Vec{3}([a3,b3,b3]),
                  Tensors.Vec{3}([b3,a3,a3]),
                  Tensors.Vec{3}([b3,a3,b3]),
                  Tensors.Vec{3}([b3,b3,a3])]

        weigths = [w1,w2,w2,w2,w2,w3,w3,w3,w3,w3,w3]
    else
        throw(ArgumentError("Gauss rule of order $order for tetrahedron not available"))
    end
    return QuadratureRule{Tetrahedron,3,Float64}(weigths, points)
end
