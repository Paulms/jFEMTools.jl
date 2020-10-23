cell_to_vtkcell(::Type{LineCell}) = WriteVTK.VTKCellTypes.VTK_LINE
cell_to_vtkcell(::Type{RectangleCell}) = WriteVTK.VTKCellTypes.VTK_QUAD
cell_to_vtkcell(::Type{TriangleCell}) = WriteVTK.VTKCellTypes.VTK_TRIANGLE
cell_to_vtkcell(::Type{HexagonCell}) = WriteVTK.VTKCellTypes.VTK_POLYGON
cell_to_vtkcell(::Type{HexahedronCell}) = WriteVTK.VTKCellTypes.VTK_HEXAHEDRON
cell_to_vtkcell(::Type{TetrahedronCell}) = WriteVTK.VTKCellTypes.VTK_TETRA


"""
    vtk_grid(filename::AbstractString, grid::Grid)

Create a unstructured VTK grid from a `PolytopalMesh`. Return a `DatasetFile`
which data can be appended to, see `vtk_point_data` and `vtk_cell_data`.
"""
function vtk_grid(filename::AbstractString, mesh::AbstractPolytopalMesh{dim,T}; compress::Bool=true) where {dim,T}
    cls = WriteVTK.MeshCell[]
    for cell_idx in 1:getncells(mesh)
        celltype = cell_to_vtkcell(getCellType(mesh, cell_idx))
        push!(cls, WriteVTK.MeshCell(celltype, getverticesindices(mesh,getcell(mesh,cell_idx))))
    end
    coords = get_vertices_matrix(mesh)'
    return WriteVTK.vtk_grid(filename, coords, cls; compress=compress)
end
