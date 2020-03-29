using Test
@time begin
    @time include("test_2d_generators.jl")
    @time include("test_mesh2.jl")
    @time include("test_boundaryhandler.jl")
    @time include("test_dofhandler.jl")
    @time include("test_interpolators.jl")
    @time include("test_quads.jl")
    @time include("test_assembler.jl")
    @time include("test_2d_vem_poisson.jl")
    @time include("test_basis.jl")
    @time include("test_fe.jl")
    @time include("test_FunctionSpace.jl")
end
