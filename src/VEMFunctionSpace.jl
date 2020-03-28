struct VEMFunctionSpace{dim,T<:Real,FE<:AbstractVirtualElement,mtype <: AbstractPolytopalMesh} <: AbstractVEMFunctionSpace{dim,T,FE}
  element::FE
  mesh::mtype
  components::Int
end

getelement(fs::VEMFunctionSpace) = fs.element
getmesh(fs::VEMFunctionSpace) = fs.mesh

function Base.show(io::IO, fs::VEMFunctionSpace)
  if fs.components == 1
      println(io, "Discrete scalar Virtual function space")
  else
      println(io, "Discrete vectorial Virtual function space (components = $(fs.components))")
  end
  println(io, "Virtual Element: \n", getelement(fs))
  println(io, "Domain mesh: \n", getmesh(fs))
end

function VEMFunctionSpace(mesh::AbstractPolytopalMesh{dim,T}, elem::AbstractVirtualElement; components::Int=1) where {dim,T}  
  VEMFunctionSpace{dim,T,typeof(elem),typeof(mesh)}(elem,mesh,components)
end

function getnlocaldofs(fs::VEMFunctionSpace, cell::Cell)
  getnlocaldofs(fs.element, cell)*fs.components
end

getncomponents(fs::VEMFunctionSpace) = fs.components