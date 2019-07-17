module jFEMTools

using SparseArrays
using LinearAlgebra
import Base:@propagate_inbounds
import Tensors
import AbstractPlotting
import AbstractPlotting: Plot, default_theme, plot!, SceneLike, Theme, to_value, Point2f0, poly!

# Mesh related functions
export  rectangle_mesh, RectangleCell, TriangleCell,
        getncells, getverticesidx, getverticescoords,
        getnedges, cell_volume, cell_centroid, cell_diameter,
        mapToGlobalIdx, getvertexset, getvertexcoords,
        getnvertices, get_vertices_matrix, get_conectivity_list,
	PolytopalMesh
# Assembler
export  start_assemble, assemble!,
        end_assemble

include("mesh.jl")
include("mesh_generator.jl")
include("assembler.jl")
include("plot_recipes.jl")
end # module
