@testset "Test quadratures" begin

    using jFEMTools

    # Triangle Quadratures
    ###############################################

    # Test Strang Quadrature
    for i in 1:5
        quad_rule = jFEMTools.QuadratureRule{2,jFEMTools.RefSimplex}(jFEMTools.Strang(),i)
        @test sum(quad_rule.weights) ≈ 0.5
    end

end
