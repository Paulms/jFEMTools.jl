@testset "Test interpolators" begin

    using jFEMTools
    import Tensors
    jF = jFEMTools

    # Monomials
    ###############################################
    P1_basis = jF.Monomials(2,1,Tensors.Vec{2}((0.0,0.0)),1.0)
    @test jF.value(P1_basis,1,Tensors.Vec{2}((1.0,2.0))) == 1.0
    @test jF.value(P1_basis,2,Tensors.Vec{2}((1.0,2.0))) == 1.0
    @test jF.value(P1_basis,3,Tensors.Vec{2}((1.0,2.0))) == 2.0
    @test jF.getnbasefunctions(P1_basis) == 3

end
