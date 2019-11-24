using jFEMTools
import Tensors

jF = jFEMTools

k = [Tensors.Vec{2}([1.0, 0.0]), Tensors.Vec{2}([0.0, 0.0]), Tensors.Vec{2}([0.0, 1.0])]

points = [Tensors.Vec{2}([1.0/6.0, 1.0/6.0]),
          Tensors.Vec{2}([1.0/6.0, 2.0/3.0]),
          Tensors.Vec{2}([2.0/3.0, 1.0/6.0])]

ref_verts = jF.reference_coordinates(jF.RefSimplex, Val{2})

pt = jF.mapPointsFromReference(jF.RefSimplex,k,points)

A,b = jF.get_affine_map(ref_verts, k)

x = Tensors.Vec{2}([0.0,1.0])
A*x + b
