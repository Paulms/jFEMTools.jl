using Test
@time begin
    @time include("test_2d_generators.jl")
    @time include("test_2d_vem_poisson.jl")
end
