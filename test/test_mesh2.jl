@testset "Mesh" begin
    using jFEMTools
    import Tensors
    const jF = jFEMTools
    geometry = Dict((2,0) => jF.MeshConectivity((1,2,3,2,4,3),(1,4,7)))
    mesh2 = jF.PolytopalMesh2((4,5,2),[Tensors.Vec{2}((0.0,0.0)),Tensors.Vec{2}((1.,0.)),Tensors.Vec{2}((0.,1.)),Tensors.Vec{2}((1.,1.))],
            geometry,Dict{Int,Dict{String,Set{Int}}}())

    mesh = unitSquareMesh(TriangleCell, (1,1));

    #Tests
    jF.get_connectivity!(mesh2,2,2);
    jF.get_connectivity!(mesh2,2,1);
    jF.get_connectivity!(mesh2,1,2);
    jF.get_connectivity!(mesh2,1,0);
    jF.get_connectivity!(mesh2,0,1);

    @test all(Set(x) in [Set((1,2,3)),Set((2,3,4))] for x in jF._unpack_connectivity(mesh2.geometry[(2,0)]))
    @test all(Set(x) in [Set((1,2,3)),Set((2,4,5))] for x in jF._unpack_connectivity(mesh2.geometry[(2,1)]))
    @test all(Set(x) in [Set((2)),Set((1))] for x in jF._unpack_connectivity(mesh2.geometry[(2,2)]))
    @test all(Set(x) in [Set((1,2)),Set((2,3)),Set((3,1)),Set((3,4)),Set((4,2))] for x in jF._unpack_connectivity(mesh2.geometry[(1,0)]))
    @test all(Set(x) in [Set((1,3)),Set((1,2,5)),Set((3,2,4)),Set((4,5))] for x in jF._unpack_connectivity(mesh2.geometry[(0,1)]))

    # Operations
    @test jF.cell_volume(mesh,1) == 0.5
    @test jF.cell_volume(mesh2,1) == 0.5
    @test jF.cell_centroid(mesh,2) ≈ Tensors.Vec{2}((2/3,2/3))
    @test jF.cell_centroid(mesh2,2) ≈ Tensors.Vec{2}((2/3,2/3))
    @test jF.cell_diameter(mesh,1) ≈ sqrt(2)
    @test jF.cell_diameter(mesh2,1) ≈ sqrt(2)
end


