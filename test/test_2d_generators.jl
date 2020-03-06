#@testset "Test 2d mesh generators" begin
using jFEMTools
import Tensors
jF = jFEMTools
#Test generated mesh
mesh = rectangle_mesh(TriangleCell, (2,2), Tensors.Vec{2}((0.0,0.0)), Tensors.Vec{2}((1.0,1.0)))
@test getncells(mesh) == 8
@test getnvertices(mesh) == 9
for cell_idx in 1:getncells(mesh)
    @test cell_volume(mesh, cell_idx) ≈ 1/8
end
@test get_vertices_matrix(mesh) == [0.0 0.0; 0.5 0.0; 1.0 0.0; 0.0 0.5; 0.5 0.5; 1.0 0.5; 0.0 1.0; 0.5 1.0; 1.0 1.0]
@test get_cell_connectivity_list(mesh) == [(1, 2, 4), (2, 5, 4), (2, 3, 5), (3, 6, 5), (4, 5, 7), (5, 8, 7), (5, 6, 8), (6, 9, 8)]
# Check expected data for cell 1
@test mesh.cells[1].vertices == (1,2,4)
#@test mesh.cells[1].faces == (1,2,3)
#@test [face_orientation(mesh,1,i) for i in 1:3] == [true,false,true]
@test cell_diameter(mesh,1) == sqrt(2)/2
@test getnedges(mesh.cells[1]) == 3


# Quad mesh
mesh = rectangle_mesh(RectangleCell, (2,2), Tensors.Vec{2}((0.0,0.0)), Tensors.Vec{2}((1.0,1.0)))

# Hex mesh
mesh = unitSquareMesh(HexagonCell,(3,3))
@test getncells(mesh) == 18
@test getnvertices(mesh) == 38
# Test boundary sets
@test jF.getedgeset(mesh,"top") == Set((jF.EdgeIndex(15,3),jF.EdgeIndex(16,4),jF.EdgeIndex(17,4),jF.EdgeIndex(18,3)))
# Check expected data for cell 1
@test getnedges(mesh.cells[1]) == 4

# Mixed mesh
cells = [
	TriangleCell((1,2,4)),
	RectangleCell((2,3,6,5)),
	TriangleCell((5,4,2)),
    RectangleCell((4, 5, 8, 7)),
	RectangleCell((5,6,9,8))
    ];
mixmesh = PolytopalMesh(cells, mesh.vertices)

#end
