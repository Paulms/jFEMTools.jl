import Tensors
using jFEMTools

nel = (2,2)
LL = Tensors.Vec{2}((0.0,0.0))
UR = Tensors.Vec{2}((1.0,3.0))
mesh = rectangle_mesh(HexagonCell,(1,3),LL,UR);

#Plot mesh
using Makie
import AbstractPlotting
include("src/plot_recipes.jl")
#popdisplay(AbstractPlotting.PlotDisplay())
#AbstractPlotting.inline!(true)
scene = Scene(resolution = (400, 200))
plot!(scene, mesh)
