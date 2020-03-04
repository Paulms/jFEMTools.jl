@testset "Test BC handlers" begin

    using jFEMTools
    using SparseArrays
    using LinearAlgebra

    mesh = unitSquareMesh(TriangleCell, (2,2));
    dim = 2
    element = PoissonVirtualElement(dim,1);
    u = TrialFunction(element)
    dh = DofHandler(mesh, u);

    # ### Boundary conditions
    dbc = Dirichlet(dh,u,"boundary",x -> 0);
    @test dbc.prescribed_dofs == [1,2,3,5,6,7,8,9]
    @test dbc.values == zeros(8)
    K = sparse(I(9)); b = ones(9);
    apply!(K,b,dbc);
    @test b == [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    @test K == sparse(I(9))
end
