
# Use map_func to avoid float precision differences between the same vertex
function _push_cell_vertex!(vert,cell_verts, used_vertices,vertices,nextvert, map_func)
  token = ht_keyindex2!(used_vertices, vert)
  if token > 0 # reuse dofs
      reuse_vert = used_vertices.vals[token]
      push!(cell_verts, reuse_vert)
  else # token <= 0, use new vertex
      Base._setindex!(used_vertices, nextvert, vert, -token)
      push!(cell_verts, nextvert)
      push!(vertices, map_func(vert))
      nextvert += 1
  end
  return nextvert
end

function _generate_2d_hex_centroids!(centroids::Vector{Tensors.Vec{2,T}},LL, n_centroid_rows, n_centroid_cols, hex_width, hex_heigth) where {T}
  x_coord = LL[1]; y_coord = LL[2]
  x_lim = n_centroid_cols*hex_width
  y_lim = n_centroid_rows*hex_heigth
  ss = 1
  for j in 1:(2*n_centroid_rows)
      for i in 1:(n_centroid_cols+1)
          centroid = Tensors.Vec{2,T}((x_coord,y_coord))
          push!(centroids, centroid)
          (x_coord + hex_width) > x_lim ? break : x_coord = x_coord + hex_width
      end
      (ss > 0) ? x_coord = LL[1] + hex_width/2 : x_coord = LL[1]
      ss = ss*(-1)
      (y_coord + hex_heigth/4) > y_lim ? break : y_coord = y_coord + 3/4*hex_heigth;
  end
  return centroids
end