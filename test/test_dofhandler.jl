@testset "Test Dof handlers" begin

using jFEMTools
using SparseArrays

mesh = unitSquareMesh(TriangleCell, (2,2));
dim = 2
element = PoissonVirtualElement(dim,1);
u = TrialFunction(element)
dh = DofHandler(mesh, u);

@test dh.cell_dofs == [1, 2, 3, 2, 4, 3, 2, 5, 4, 5, 6, 4, 3, 4, 7, 4, 8, 7, 4, 6, 8, 6, 9, 8]
@test dh.cell_dofs_offset == [1, 4, 7, 10, 13, 16, 19, 22, 25]
@test jFEMTools.dof_range(dh, u, 1) == 1:3
@test jFEMTools.vertexdofs(dh, u) == [1,2,5,3,4,6,7,8,9]
K = jFEMTools.create_sparsity_pattern(dh)
@test size(K) == (9,9)
@test nnz(K) == 41
end
