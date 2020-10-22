import Tensors
using jFEMTools
import WriteVTK

nel = (2,2)
LL = Tensors.Vec{2}((0.0,0.0))
UR = Tensors.Vec{2}((1.0,3.0))
mesh = rectangle_mesh(HexagonCell,(1,3),LL,UR);
vtk_file = vtk_grid("hexagon_mesh", mesh)
outfiles = WriteVTK.vtk_save(vtk_file)


#Plot mesh
using Makie
import AbstractPlotting
include("src/plot_recipes.jl")
#popdisplay(AbstractPlotting.PlotDisplay())
#AbstractPlotting.inline!(true)
scene = Scene(resolution = (400, 200))
plot!(scene, mesh)

nel = (2,2,2)
mesh3D_1 = hyper_rectagle_mesh2(HexahedronCell,nel)
vtk_file = vtk_grid("hexahedron_mesh", mesh3D_1, compress=false)
outfiles = WriteVTK.vtk_save(vtk_file)

mesh3D_2 = hyper_rectagle_mesh2(TetrahedronCell,nel)
vtk_file = vtk_grid("tetrahedron_mesh", mesh3D_2, compress=false)
outfiles = WriteVTK.vtk_save(vtk_file)