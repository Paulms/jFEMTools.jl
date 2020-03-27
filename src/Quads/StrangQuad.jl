function (::Type{QuadratureRule{Triangle}})(quad_type::Strang, order::Int)
    if order == 0 || order == 1
        points = [Tensors.Vec{2}([1.0/3.0, 1.0/3.0])]
        weigths = [0.5]
    elseif order == 2
        points = [Tensors.Vec{2}([1.0/6.0, 1.0/6.0]),
                  Tensors.Vec{2}([1.0/6.0, 2.0/3.0]),
                  Tensors.Vec{2}([2.0/3.0, 1.0/6.0])]
        weigths = [1/6,1/6,1/6]
    elseif order == 3
        points = [Tensors.Vec{2}([0.659027622374092, 0.231933368553031]),
                  Tensors.Vec{2}([0.659027622374092, 0.109039009072877]),
                  Tensors.Vec{2}([0.231933368553031, 0.659027622374092]),
                  Tensors.Vec{2}([0.231933368553031, 0.109039009072877]),
                  Tensors.Vec{2}([0.109039009072877, 0.659027622374092]),
                  Tensors.Vec{2}([0.109039009072877, 0.231933368553031])]
        weigths = 1/12*ones(Float64,6)
    elseif order == 4
        points = [Tensors.Vec{2}([0.816847572980459, 0.091576213509771]),
                  Tensors.Vec{2}([0.091576213509771, 0.816847572980459]),
                  Tensors.Vec{2}([0.091576213509771, 0.091576213509771]),
                  Tensors.Vec{2}([0.108103018168070, 0.445948490915965]),
                  Tensors.Vec{2}([0.445948490915965, 0.108103018168070]),
                  Tensors.Vec{2}([0.445948490915965, 0.445948490915965])]
        weigths = 1/2*ones(Float64,6)
        weigths[1:3] = 0.109951743655322*weigths[1:3]
        weigths[4:6] = 0.223381589678011*weigths[4:6]
    elseif order == 5
        points = [Tensors.Vec{2}([0.33333333333333333, 0.33333333333333333]),
                  Tensors.Vec{2}([0.79742698535308720, 0.10128650732345633]),
                  Tensors.Vec{2}([0.10128650732345633, 0.79742698535308720]),
                  Tensors.Vec{2}([0.10128650732345633, 0.10128650732345633]),
                  Tensors.Vec{2}([0.05971587178976981, 0.47014206410511505]),
                  Tensors.Vec{2}([0.47014206410511505, 0.05971587178976981]),
                  Tensors.Vec{2}([0.47014206410511505, 0.47014206410511505])]
        weigths = 1/2*ones(Float64,7)
        weigths[1] = 0.22500000000000000*weigths[1]
        weigths[2:4] = 0.12593918054482717*weigths[2:4]
        weigths[5:7] = 0.13239415278850616*weigths[5:7]
    elseif order == 6
        points =[Tensors.Vec{2}([0.873821971016996, 0.063089014491502]),
                 Tensors.Vec{2}([0.063089014491502, 0.873821971016996]),
                 Tensors.Vec{2}([0.063089014491502, 0.063089014491502]),
                 Tensors.Vec{2}([0.501426509658179, 0.249286745170910]),
                 Tensors.Vec{2}([0.249286745170910, 0.501426509658179]),
                 Tensors.Vec{2}([0.249286745170910, 0.249286745170910]),
                 Tensors.Vec{2}([0.636502499121399, 0.310352451033785]),
                 Tensors.Vec{2}([0.636502499121399, 0.053145049844816]),
                 Tensors.Vec{2}([0.310352451033785, 0.636502499121399]),
                 Tensors.Vec{2}([0.310352451033785, 0.053145049844816]),
                 Tensors.Vec{2}([0.053145049844816, 0.636502499121399]),
                 Tensors.Vec{2}([0.053145049844816, 0.310352451033785])]
        weigths = 1/2*ones(Float64,12)
        weigths[1:3] = 0.050844906370207*weigths[1:3]
        weigths[4:6] = 0.116786275726379*weigths[4:6]
        weigths[7:12] = 0.082851075618374*weigths[7:12]
    else
        throw(ArgumentError("Strang rule of order $order not available"))
    end
    return QuadratureRule{Triangle,2,Float64}(weigths, points)
end
