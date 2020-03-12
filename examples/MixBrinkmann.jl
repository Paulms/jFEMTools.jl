# Based on:
# E. Cáceres, G.N. Gatica, F.A. Sequeira, A mixed virtual element method for
# the Brinkman problem, Math. Models Methods Appl. Sci. 27 (2017), no. 4,
# 707–743.

using jFEMTools
import Tensors
const jF = jFEMTools;
import FastGaussQuadrature
using LinearAlgebra

order = 2;
dim = 2;

mesh = unitSquareMesh2(RectangleCell, (3,3));

# Assamble local matrices
cell = jF.getcells(mesh)[1]
ci = 1

vertices = getverticescoords(mesh, ci)
tess = jF.get_2Dtesselation(vertices) #tesselate cell
centroid = cell_centroid(mesh, ci)
diameter = cell_diameter(mesh, ci)


# Init basis
Pe_ref_basis = jF.Monomials(dim-1,order,Tensors.Vec((0.5,)),1.0)
Pkm_basis = jF.Monomials(dim,order+1,centroid,diameter)
Pk_basis = jF.Monomials(dim,order,centroid,diameter)

points, weights = FastGaussQuadrature.gausslobatto(2*order-1)
# Shift interval from (-1,1) to (0,1)
weights *= 0.5
points = points .+ 1.0; points /= 2.0

# Assemble M_{mass,1} matrix
ϕdofs = jF.getnbasefunctions(Pe_ref_basis)

Mmass1ref = zeros(Float64,ϕdofs,ϕdofs)
for i in 1:ϕdofs, j = 1:ϕdofs, k = 1:size(points,1)
  Mmass1ref[i,j] += jF.value(Pe_ref_basis,i,Tensors.Vec((points[k],)))*jF.value(Pe_ref_basis,j,Tensors.Vec((points[k],)))*weights[k]
end
Mmass1ref_inv = inv(Mmass1ref)

# Assemble M_{edge} matrix
mdofs = jF.getnbasefunctions(Pkm_basis)
Medges = Vector{Matrix{Float64}}()
for e in 1:jF.getnedges(mesh,cell)
  Medge = zeros(Float64,ϕdofs,mdofs)
  edge_coords = jF.getverticescoords(mesh, jF.EdgeIndex(ci,e))
  he = norm(edge_coords[2] - edge_coords[1]) 
  for i in 1:ϕdofs, j = 1:mdofs, k = 1:size(points,1)
    mappedPoint = jF.map_unit_to_segment(Tensors.Vec((points[k],)),edge_coords)
    Medge[i,j] += jF.value(Pe_ref_basis,i,Tensors.Vec((points[k],)))*jF.value(Pkm_basis,j,mappedPoint)*weights[k]
  end
  push!(Medges,he*Medge)
end

# Assemble M_{mass,2}^{k+1}
Mmass2m = zeros(Float64,mdofs,mdofs)
for i in 1:mdofs, j = 1:mdofs
  αβ = Pkm_basis.indices[i] .+ Pkm_basis.indices[j]
  for e in 1:jF.getnedges(mesh,cell)
    edge_coords = jF.getverticescoords(mesh, jF.EdgeIndex(ci,e))
    normal = normal = jF.get_Normal(mesh, jF.EdgeIndex(ci, e))
    for k = 1:size(points,1)
      mappedPoint = jF.map_unit_to_segment(Tensors.Vec((points[k],)),edge_coords)
      Mmass2m[i,j] += diameter/(2+sum(αβ))*
                      dot(((1/diameter*(mappedPoint-centroid)).^αβ), normal)
    end
  end
end

#Assemble M_{mass,2} matrix
mkdofs = jF.getnbasefunctions(Pk_basis)
Mmass2 = view(Mmass2m,1:mkdofs,1:mkdofs)

#Assemble M_{grad}
Mgrad = zeros(Float64,mdofs-1,mdofs-1)
for i in 2:mdofs, j = 2:mdofs
  for e in 1:jF.getnedges(mesh,cell)
    edge_coords = jF.getverticescoords(mesh, jF.EdgeIndex(ci,e))
    normal = normal = jF.get_Normal(mesh, jF.EdgeIndex(ci, e))
    for k = 1:size(points,1)
      mappedPoint = jF.map_unit_to_segment(Tensors.Vec((points[k],)),edge_coords)
      # Product of δx derivatives
      α = Pkm_basis.indices[i] .- (1,0)
      β = Pkm_basis.indices[j] .- (1,0)
      if !(any(x -> x < 0, α) || any(x -> x < 0, β))
        Mgrad[i-1,j-1] += 1/(2+sum(α.+β))*((α[1]+β[1])/diameter)*
                        dot(((1/diameter*(mappedPoint-centroid)).^(α.+β)), normal)
      end
      # Product of δy derivatives
      α = Pkm_basis.indices[i] .- (0,1)
      β = Pkm_basis.indices[j] .- (0,1)
      if !(any(x -> x < 0, α) || any(x -> x < 0, β))
        Mgrad[i-1,j-1] += 1/(2+sum(α.+β))*((α[1]+β[1])/diameter)*
                        dot(((1/diameter*(mappedPoint-centroid)).^(α.+β)), normal)
      end
    end
  end
end

# Compute basis for G_k^{⟂}
gperpdofs = Int(order*(order+1)/2)
M0 = zeros(Float64, mdofs-1,2*mdofs)

for i in 2:mdofs, j = 1:mdofs
  for e in 1:jF.getnedges(mesh,cell)
    edge_coords = jF.getverticescoords(mesh, jF.EdgeIndex(ci,e))
    normal = normal = jF.get_Normal(mesh, jF.EdgeIndex(ci, e))
    for k = 1:size(points,1)
      mappedPoint = jF.map_unit_to_segment(Tensors.Vec((points[k],)),edge_coords)
      # δx m_{i+2} m_{j}
      α = Pkm_basis.indices[i] .- (1,0)
      β = Pkm_basis.indices[j]
      if !(any(x -> x < 0, α))
        M0[i-1,j] += (α[1]+1)/(2+sum(α.+β))*
                        dot(((1/diameter*(mappedPoint-centroid)).^(α.+β)), normal)
      end
      # δy m_{i+2} m_{j}
      α = Pkm_basis.indices[i] .- (0,1)
      β = Pkm_basis.indices[j]
      if !(any(x -> x < 0, α))
        M0[i-1,j+mdofs] += (α[2]+1)/(2+sum(α.+β))*
                        dot(((1/diameter*(mappedPoint-centroid)).^(α.+β)), normal)
      end
    end
  end
end

A0 = nullspace(M0)[:,1:gperpdofs]