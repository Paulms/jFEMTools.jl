module jFEMTools

using SparseArrays
using LinearAlgebra
import Base:@propagate_inbounds
import Tensors
import Tensors: Vec,norm,dot,gradient,det, ⊗
import Base.ht_keyindex2!
import FastGaussQuadrature
import VoronoiDelaunay
import PlanarConvexHulls
import StaticArrays

# Abstract types
abstract type AbstractElement end
abstract type AbstractPolytopalMesh{dim,T} end

function Base.show(io::IO, mesh::AbstractPolytopalMesh{dim}) where {dim}
  println(io, "$dim D PolytopalMesh")
  println(io, "type: ", typeof(mesh))
  println(io, "Number of cells: ", getncells(mesh))
  println(io, "Number of vertices: ", getnvertices(mesh))
end

# Mesh related functions
export  rectangle_mesh, RectangleCell, TriangleCell, HexagonCell,
        getncells, getverticesidx, getverticescoords,
        getnedges, cell_volume, cell_centroid, cell_diameter,
        mapToGlobalIdx, getvertexset, getvertexcoords,
        getnvertices, get_vertices_matrix, get_cell_connectivity_list,
  PolytopalMesh, unitSquareMesh, getnfacets
export FaceIndex, EdgeIndex, FacetIndex
export PolytopalMesh2, rectangle_mesh2, unitSquareMesh2
# Assembler
export  start_assemble, assemble!,
        end_assemble
# TrialFunctions
export TrialFunction

# Element
export PoissonVirtualElement,LocalPoissonVirtualElement

#VEM utils
export DofHandler
export VEMOperators
export assemble_stiffnessMat, assemble_massMat, assemble_load
export Dirichlet, apply!

#FEM
export ContinuousLagrange
export FEFunctionSpace

include("FEM/shapes.jl")
include("Quads/quadrature.jl")
include("tools.jl")
include("Quads/StrangQuad.jl")
include("Quads/GrundmannMoellerQuad.jl")
include("Interpolations.jl")
include("mesh_operations.jl")
include("mesh.jl")
include("mesh2.jl")
include("mesh_generic_funcs.jl")
include("mesh_generators.jl")
include("mesh2_generators.jl")
include("functions.jl")
include("assembler.jl")
#include("plot_recipes.jl")
include("PoissonVirtualElement.jl")
include("DofHandler.jl")
include("VEMOperators.jl")
include("boundary.jl")
include("FEM/basis.jl")
include("FEM/FiniteElement.jl")
include("FEM/LagrangeFE.jl")
include("FEM/FunctionSpaces.jl")
include("FEM/iterator.jl")

# SnoopCompile output
include("precompile.jl")
_precompile_()

end # module
