#@testset "Mesh" begin
    using jFEMTools
    import Tensors
    const jF = jFEMTools
    geometry = Dict((2,0) => jF.MeshConectivity((1,2,3,2,4,3),(1,4,7)))
    mesh1 = jF.PolytopalMesh2((4,5,2),[Tensors.Vec{2}((0.0,0.0)),Tensors.Vec{2}((1.,0.)),Tensors.Vec{2}((0.,1.)),Tensors.Vec{2}((1.,1.))],
            geometry,Dict{Int,Dict{String,Set{Int}}}())

    mesh2 = unitSquareMesh(TriangleCell, (1,1));

    #Tests
    jF.get_connectivity!(mesh1,2,2);
    jF.get_connectivity!(mesh1,2,1);
    jF.get_connectivity!(mesh1,1,2);
    jF.get_connectivity!(mesh1,1,0);
    jF.get_connectivity!(mesh1,0,1);

    @test all(Set(x) in [Set((1,2,3)),Set((2,3,4))] for x in jF._unpack_connectivity(mesh1.geometry[(2,0)]))
    @test all(Set(x) in [Set((1,2,3)),Set((2,4,5))] for x in jF._unpack_connectivity(mesh1.geometry[(2,1)]))
    @test all(Set(x) in [Set((2)),Set((1))] for x in jF._unpack_connectivity(mesh1.geometry[(2,2)]))
    @test all(Set(x) in [Set((1,2)),Set((2,3)),Set((3,1)),Set((3,4)),Set((4,2))] for x in jF._unpack_connectivity(mesh1.geometry[(1,0)]))
    @test all(Set(x) in [Set((1,3)),Set((1,2,5)),Set((3,2,4)),Set((4,5))] for x in jF._unpack_connectivity(mesh1.geometry[(0,1)]))
#end

# Operations
edge = EdgeIndex(2,2)
jF.get_Normal(mesh2,edge)
jF.get_Normal(mesh1,edge)

jF.get_vertices_matrix(mesh2)
jF.get_vertices_matrix(mesh1)
jF.cell_volume(mesh2,1)
jF.cell_volume2(mesh2,1)
jF.cell_volume2(mesh1,1)

using BenchmarkTools

@btime jF.cell_volume(mesh2,1)
@btime jF.cell_volume2(mesh2,1)
@btime jF.cell_volume2(mesh1,1)


function cell_centroid(mesh::PolytopalMesh2{2}, cell_idx::Int)
  cell = MeshEntity{2}(cell_idx)
  verts = getverticescoords(mesh,cell)
  vertices = [StaticArrays.SVector(x[1],x[2]) for x in verts]
  chull = PlanarConvexHulls.ConvexHull{PlanarConvexHulls.CCW}(vertices)
  return Tensors.Vec{2}(PlanarConvexHulls.centroid(chull))
end

function cell_diameter(mesh::PolytopalMesh2{dim,T}, cell_idx::Int) where {dim,T}
  cell = MeshEntity{dim}(cell_idx)
  verts = getverticescoords(mesh,cell)
  vertices = [StaticArrays.SVector(x[1],x[2]) for x in verts]
  chull = PlanarConvexHulls.ConvexHull{PlanarConvexHulls.CCW}(vertices)
  h = 0.0
  for vert in chull.vertices
      h = max(h,maximum([norm(vert-vert2) for vert2 in chull.vertices]))
  end
  h
end
