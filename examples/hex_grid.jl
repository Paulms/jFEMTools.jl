import Tensors
using jFEMTools

nel = (3,3)
LL = Tensors.Vec{2}((0.0,0.0))
UR = Tensors.Vec{2}((2.0,2.0))
mesh = rectangle_mesh(HexagonCell,(3,3),LL,UR)

#Plot mesh
using Makie
import AbstractPlotting
include("../src/plot_recipes.jl")
#popdisplay(AbstractPlotting.PlotDisplay())
#AbstractPlotting.inline!(true)
scene = Scene(resolution = (400, 200))
plot!(scene, mesh)
