module jFEMTools

using SparseArrays
using LinearAlgebra
import Base:@propagate_inbounds
import Tensors
import Base.ht_keyindex2!
import FastGaussQuadrature
import VoronoiDelaunay
import PlanarConvexHulls
import StaticArrays

# Abstract types
abstract type AbstractElement end

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
# TrialFunctions
export TrialFunction

# Element
export VirtualElement,LocalVirtualElement

#VEM utils
export DofHandler
export VEMOperators
export assemble_K, assemble_M, assemble_load
export Dirichlet, apply!

include("quadrature.jl")
include("tools.jl")
include("StrangQuad.jl")
include("Interpolations.jl")
include("mesh.jl")
include("mesh_generator.jl")
include("functions.jl")
include("assembler.jl")
#include("plot_recipes.jl")
include("VirtualElement.jl")
include("DofHandler.jl")
include("VEMOperators.jl")
include("boundary.jl")

# SnoopCompile output
include("precompile.jl")
_precompile_()

end # module
