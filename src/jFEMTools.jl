module jFEMTools

using SparseArrays
using LinearAlgebra
import Base:@propagate_inbounds
import Tensors
import AbstractPlotting
import AbstractPlotting: Plot, default_theme, plot!, SceneLike, Theme, to_value, Point2f0, poly!
import Base.ht_keyindex2!
import FastGaussQuadrature
import VoronoiDelaunay

# Abstract types
abstract type AbstractQuadratureRule end

# Mesh related functions
export  rectangle_mesh, RectangleCell, TriangleCell,
        getncells, getverticesidx, getverticescoords,
        getnedges, cell_volume, cell_centroid, cell_diameter,
        mapToGlobalIdx, getvertexset, getvertexcoords,
        getnvertices, get_vertices_matrix, get_conectivity_list,
	PolytopalMesh, unitSquareMesh
# Assembler
export  start_assemble, assemble!,
        end_assemble
# Element
export VirtualElement,LocalVirtualElement

#VEM utils
export DofHandler
include("tools.jl")
include("StrangQuad.jl")
include("Interpolations.jl")
include("mesh.jl")
include("mesh_generator.jl")
include("assembler.jl")
include("plot_recipes.jl")
include("VirtualElement.jl")
include("DofHandler.jl")
end # module
