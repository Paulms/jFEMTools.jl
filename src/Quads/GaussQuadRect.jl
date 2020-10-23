# Generate Gauss quadrature rules on cubes (0,0,...) - (1,1,...) by doing an outer product
# over all dimensions

import Base.Cartesian: @nloops, @nref, @ntuple, @nexprs

for dim in (1,2,3)
    @eval begin
        function (::Type{QuadratureRule{HyperCube{$dim}}})(quad_type::GaussLegendre, order::Int)
            p,w = FastGaussQuadrature.gausslegendre(order)
            weights = Vector{Float64}(undef, order^($dim))
            points = Vector{Vec{$dim,Float64}}(undef, order^($dim))
            count = 1
            @nloops $dim i j->(1:order) begin
                t = @ntuple $dim q-> 0.5*(p[$(Symbol("i"*"_q"))] .+ 1.0)
                points[count] = Vec{$dim,Float64}(t)
                weight = 1.0
                @nexprs $dim j->(weight *= w[i_{j}])
                weights[count] = weight
                count += 1
            end
            return QuadratureRule{HyperCube{$dim},$dim,Float64}(weights./sum(weights), points)
        end
    end
end