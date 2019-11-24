using Test
@time begin
    @time include("test_2d_generators.jl")
    @time include("test_boundaryhandler.jl")
    @time include("test_dofhandler.jl")
    @time include("test_interpolators.jl")
    @time include("test_quads.jl")
    @time include("test_2d_vem_poisson.jl")
end
