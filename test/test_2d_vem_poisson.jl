@testset "Test 2d vem poisson" begin

    using jFEMTools
    import Tensors

    mesh = unitSquareMesh(RectangleCell, (1,1))

    # forcing function
    rhs(x::Tensors.Vec{2}) = 2*π^2*sin(π*x[1])*sin(π*x[2])

    # Test for order = 1
    degree = 1;
    dim = 2;
    element = PoissonVirtualElement(dim,degree);
    Vh = VEMFunctionSpace(mesh,element)
    u = TrialFunction(Vh)
    dofs = DofHandler(mesh, [u]);
    operators = VEMOperators(dofs, u;load = rhs);
    # Test
    d = sqrt(2)
    Bref = 1/4*[1 1 1 1;
                -d d -d d;
                -d -d d d]
    @test Bref ≈ operators.elements[1].B

    Dref = 1/4*[4 -d -d;
                4 d -d;
                4 d  d;
                4 -d d]

    @test Dref ≈ operators.elements[1].D

    Gref = 1/2*[2 0 0;
                0 1 0;
                0 0 1]

    @test Gref ≈ operators.elements[1].G

    # Test order 2 implementation
    degree = 2
    dim = 2
    element = PoissonVirtualElement(dim,degree)
    u = TrialFunction(VEMFunctionSpace(mesh,element))
    dofs = DofHandler(mesh, [u])
    operators = VEMOperators(dofs, u;load = rhs)
    # Test
    d = sqrt(2)
    Bref = 1/12*[0 0 0 0 0 0 0  0 12;
                 -d d d -d 0 4*d 0 -4*d 0;
                 -d -d d d -4*d 0 4*d 0 0;
                 1 1 1 1 0 4 0 4 -12;
                 1 -1 1 -1 0 0 0 0 0;
                 1 1 1 1 4 0 4 0 -12]
    @test Bref ≈ operators.elements[1].B

    Dref = 1/24*[24 -6*d -6*d 3 3 3;
                24 6*d -6*d 3 -3 3;
                24 6*d 6*d 3 3 3;
                24 -6*d 6*d 3 -3 3;
                24 0 -6*d 0 0 3;
                24 6*d 0 3 0 0;
                24 0 6*d 0 0 3;
                24 -6*d 0 3 0 0;
                24 0 0 1 0 1]

    @test Dref ≈ operators.elements[1].D

    Gref = 1/24*[24 0 0 1 0 1;
                0 12 0 0 0 0;
                0 0 12 0 0 0;
                0 0 0 2 0 0;
                0 0 0 0 1 0;
                0 0 0 0 0 2]

    @test Gref ≈ operators.elements[1].G
end
