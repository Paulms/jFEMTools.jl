@testset "Test quadratures" begin

    using jFEMTools
    jF = jFEMTools

    # Triangle Quadratures
    ###############################################

    # Test Strang Quadrature
    for i in 1:5
        quad_rule = jF.QuadratureRule{jF.Triangle}(jF.Strang(),i)
        @test sum(quad_rule.weights) ≈ 0.5
    end

end
