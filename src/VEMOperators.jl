struct LocalVEMOperators{T}
  order::Int
  D::Matrix{T}
  G::Matrix{T}
  B::Matrix{T}
  H::Matrix{T}
  b::Vector{T}
end

get_order(lvem::LocalVEMOperators) = lvem.order

struct VEMOperators
    dofs::DofHandler
    elements::Vector{LocalVEMOperators}
end

struct CellCache{ct}
  cell::ct
  ci::Int
  tess
  quad
  f
end

function VEMOperators(dof::DofHandler, u::TrialFunction; load = (x->0.0)) where {dim,T}
    elements = get_local_elements(dof, getelement(getfunctionspace(u)), load)
    VEMOperators(dof,elements)
end

function get_local_elements(dof::DofHandler{2}, element::PoissonVirtualElement{2}, load)
    local_elements = Vector{LocalVEMOperators}()
    for (ci, cell) in enumerate(getcells(dof.mesh))
        vertices = getverticescoords(dof.mesh, ci)
        tess = get_2Dtesselation(vertices) #tesselate cell
        quad = QuadratureRule{Triangle}(Strang(),get_degree(element))
        cache = CellCache(cell, ci, tess, quad,load)
        centroid = cell_centroid(dof.mesh, ci)
        diameter = cell_diameter(dof.mesh, ci)
        lelement = LocalPoissonVirtualElement(2,get_degree(element),centroid,diameter)
        D = _compute_local_D(lelement, dof, cache)
        G = _compute_local_G(lelement, dof, cache)
        B = _compute_local_B(lelement, dof, cache)
        H = _compute_local_H(lelement, dof, cache)
        b = _compute_local_b(lelement, dof, cache)
        push!(local_elements, LocalVEMOperators(get_degree(element),D,G,B,H,b))
    end
    local_elements
end

# NOTE: From here the methods work only for 2D
""" Compute local D matrix """
function _compute_local_D(element::LocalPoissonVirtualElement, dof::DofHandler{2}, cache::CellCache{CT}) where {CT}
    cell = cache.cell; ci = cache.ci
    degree = get_degree(element)
    ndofs = ndofs_per_cell(dof, ci)
    nk = Int((degree+1)*(degree+2)/2)
    D = zeros(ndofs,nk)
    nv = getnvertices(dof.mesh,cell)
    vertices = getverticescoords(dof.mesh, ci)
    # compute value in vertices
    for i in 1:nv, j in 1:nk
            D[i,j] = value(element.Pk_basis, j, vertices[i])
    end
    if degree >= 2
        ne = getnedges(dof.mesh,cell)
        # Compute value for edges
        for i in (nv+1):(nv+ne), j in 1:nk
            D[i,j] = value(element.Pk_basis, j,
                map_unit_to_segment(element.edge_quad.points[2],
                        getverticescoords(dof.mesh, EdgeIndex(ci,i-nv))))
        end
        #Compute the Internal dofs
        nf = Int((degree-1)*degree/2)   #dimension of Pk_basis with order k-2
        lap_basis = Monomials(2,degree-2,element.Pk_basis.centroid, element.Pk_basis.diameter)
        for i in (nv+ne*(degree-1)+1):ndofs, j in 1:nk
             for k in cache.tess
               pt = mapPointsFromReference(Triangle,k,cache.quad.points);
               for g in 1:size(pt,1)
                 D[i,j] = D[i,j] + 2*simplex_area(k)*
                   cache.quad.weights[g]*
                   (value(element.Pk_basis, j, pt[g])*
                    value(lap_basis, i-(nv+ne*(degree-1)),
                   pt[g]))/cell_volume(dof.mesh, ci)
               end
             end
           end
         end
    return D
end

function p0m(element::LocalPoissonVirtualElement, dof::DofHandler{2},cache::CellCache{CT}, i::Int) where {CT}
  cell = cache.cell; ci = cache.ci
  degree = get_degree(element)
  nv = getnvertices(dof.mesh,cell)
  P0m = 0.0
   if degree == 1
     vertices = getverticescoords(dof.mesh, ci)
     for j in 1:nv
       P0m = P0m + value(element.Pk_basis, i, vertices[j])
     end
     P0m = P0m/nv
   else
     for k in cache.tess
       pt = mapPointsFromReference(Triangle,k,cache.quad.points);
       for g in 1:size(pt,1)
         P0m = P0m + 2*simplex_area(k)*
           cache.quad.weights[g]*value(element.Pk_basis, i, pt[g])
       end
     end
     P0m = P0m/cell_volume(dof.mesh, ci)
   end
   return P0m
end

""" Compute local G matrix """
function _compute_local_G(element::LocalPoissonVirtualElement, dof::DofHandler{2}, cache::CellCache{CT}) where {CT}
  cell = cache.cell; ci = cache.ci
  degree = get_degree(element)
  nk = Int((degree+1)*(degree+2)/2)
  ne = getnedges(dof.mesh,cell)
  G = zeros(nk, nk );
  # Gradient terms
  # First row
  for j in 1:nk
    G[1,j] = p0m(element, dof, cache, j);
  end
  if degree == 1
     # No need to integrate over the domain in case of k=1
     for i in 2:nk, j in 2:nk
         for k in 1:ne
           v = getverticescoords(dof.mesh, EdgeIndex(ci, k))
           normal = get_Normal(dof.mesh, EdgeIndex(ci, k))
           distance = norm(v[2]-v[1])

           # Use trapezium rule
           grad_m_i_v1 = gradient_value(element.Pk_basis, i, v[1])
           m_j_v1 = value(element.Pk_basis, j, v[1])

           grad_m_i_v2 = gradient_value(element.Pk_basis, i, v[2])
           m_j_v2 = value(element.Pk_basis, j, v[2])

           G[i,j] = G[i,j] + distance/2 * (
             dot(grad_m_i_v1,normal)*m_j_v1 +
             dot(grad_m_i_v2,normal)*m_j_v2 )
         end
     end
   else
     for i in 2:nk, j in 2:nk
         for k in cache.tess
           pt = mapPointsFromReference(Triangle,k,cache.quad.points);
           for g in 1:size(pt,1)
             G[i,j] = G[i,j] + 2*simplex_area(k)*
                cache.quad.weights[g]*
                dot(gradient_value(element.Pk_basis, i, pt[g]),
                    gradient_value(element.Pk_basis, j, pt[g]),)
           end
         end
       end
   end
   return G
end

function compute_d_αβ(element::LocalPoissonVirtualElement, dof::DofHandler{2},cache::CellCache{CT}, k::Int) where {CT}
  cell = cache.cell; ci = cache.ci
  degree = get_degree(element)
  nv = getnvertices(dof.mesh,cell)
  nkm2 = Int((degree-1)*degree/2)
  # Dimension of M_{k-2}: n_{k-2}
  d = zeros(nkm2,1);
  h = cell_diameter(dof.mesh, ci)
  # Build dx^2 monomial
  dxx_monomial_index = element.Pk_basis.indices[k][1]-2
  laplacian_indices = get_monomial2DIndices(degree - 2)
  if (dxx_monomial_index >= 0)
    for i in 1:nkm2
      if (laplacian_indices[i] == (dxx_monomial_index, element.Pk_basis.indices[k][2]))
        d[i] = d[i] + element.Pk_basis.indices[k][1]*(element.Pk_basis.indices[k][1]-1)/h^2;
      end
    end
  end
  # Build dy^2 monomial
  dyy_monomial_index = element.Pk_basis.indices[k][2]-2
  if (dyy_monomial_index >= 0)
    for i in 1:nkm2
      if (laplacian_indices[i] == (element.Pk_basis.indices[k][1],dyy_monomial_index))
        d[i] = d[i] + element.Pk_basis.indices[k][2]*(element.Pk_basis.indices[k][2]-1)/h^2;
      end
    end
  end
  return d
end

""" Compute the local B matrix """
function _compute_local_B(element::LocalPoissonVirtualElement, dof::DofHandler{2}, cache::CellCache{CT}) where {CT}
  cell = cache.cell; ci = cache.ci
  degree = get_degree(element)
  nk = Int((degree+1)*(degree+2)/2)
  ne = getnedges(dof.mesh,cell)
  ndofs = ndofs_per_cell(dof, ci)
  nv = getnvertices(dof.mesh,cell)
  B = zeros(nk, ndofs);

   if (degree == 1)
     # First row
     B[1,:] .= 1/nv
     for i in 2:nk, j in 1:ndofs
         for k in 1:ne
           v = getverticescoords(dof.mesh, EdgeIndex(ci, k))
           normal = get_Normal(dof.mesh, EdgeIndex(ci, k))
           distance = norm(v[2]-v[1])

           # Use trapezium rule
           grad_m_i_v1 = gradient_value(element.Pk_basis, i, v[1])
           φ_j_v1 = δ(j, getverticesindices(dof.mesh, EdgeIndex(ci,k))[1])

           grad_m_i_v2 = gradient_value(element.Pk_basis, i, v[2])
           φ_j_v2 = δ(j, getverticesindices(dof.mesh, EdgeIndex(ci,k))[2])

           B[i,j] = B[i,j] + distance/2 * (
             dot( grad_m_i_v1,normal ) * φ_j_v1 +
             dot( grad_m_i_v2,normal ) * φ_j_v2 )

         end
     end
   else
     # First row
     for j in (nv*degree+1):ndofs
       B[1,j] = p0m(element,dof,cache,j-(nv*degree))
     end
     for i in 2:nk
       d_αβ = compute_d_αβ(element,dof,cache,i)
       for j in 1:ndofs
         # First compute the components owning the Laplacian
         if j > (nv*degree)
           B[i,j] = B[i,j] -
                cell_volume(dof.mesh, ci)*d_αβ[j - (nv*degree)]
         end
         # Now the boundary part
         for k in 1:ne
           v = getverticescoords(dof.mesh, EdgeIndex(ci, k))
           normal = get_Normal(dof.mesh, EdgeIndex(ci, k))
           distance = norm(v[2]-v[1])
           for g=1:size(element.edge_quad.points,1)
             grad_m_i = gradient_value(element.Pk_basis, i,
                  map_unit_to_segment(element.edge_quad.points[g],v))
             edge_dof = -1
             if g == 1
               edge_dof = reference_edge_vertices(dof.mesh,cell)[k][1]
             elseif g == size(element.edge_quad.points,1)
               edge_dof = reference_edge_vertices(dof.mesh,cell)[k][2]
             else
               edge_dof = nv + (degree-1)*(k-1)+(g-1)
             end
             φ_j = δ(j,edge_dof)

             B[i,j] = B[i,j] +
               distance *(element.edge_quad.weights[g]*dot(grad_m_i,normal) * φ_j )
           end
        end
       end
     end
   end
   return B
end

""" Compute the local H matrix """
function _compute_local_H(element::LocalPoissonVirtualElement, dof::DofHandler{2}, cache::CellCache{CT}) where {CT}
  cell = cache.cell; ci = cache.ci
  degree = get_degree(element)
  nk = Int((degree+1)*(degree+2)/2)
  ne = getnedges(dof.mesh,cell)
  ndofs = ndofs_per_cell(dof, ci)
  nv = getnvertices(dof.mesh,cell)
  H = zeros(nk, nk)

  for i in 1:nk, j in 1:nk
      for k in cache.tess
        pt = mapPointsFromReference(Triangle,k,cache.quad.points);
        for g=1:size(pt,1)
          H[i,j] = H[i,j] + 2*simplex_area(k)*
            cache.quad.weights[g]*
            ( value(element.Pk_basis, i, pt[g])* value(element.Pk_basis, j, pt[g]))
        end
      end
  end
  return H
end

""" Compute the local load vector """
function _compute_local_b(element::LocalPoissonVirtualElement, dof::DofHandler{2}, cache::CellCache{CT}) where {CT}
  cell = cache.cell; ci = cache.ci; load_func = cache.f
  degree = get_degree(element)
  if degree == 1
    ndofs = ndofs_per_cell(dof, ci)
    b = load_func(element.Pk_basis.centroid)*ones(ndofs)*cell_volume(dof.mesh, ci)/getnedges(cell)
  else
    #TODO: Should we use k-1?
    nk = Int((degree+1)*(degree+2)/2)
    b = zeros(nk)
    for i in 1:nk
        for k in cache.tess
          pt = mapPointsFromReference(Triangle,k,cache.quad.points);
          for g=1:size(pt,1)
            b[i] = b[i] + 2*simplex_area(k)*
              cache.quad.weights[g]*
              ( value(element.Pk_basis, i, pt[g])* load_func(pt[g]))
          end
        end
    end
  end
  return b
end

function assemble_load(op::VEMOperators)
      b = zeros(ndofs(op.dofs))
      for k = 1:getncells(op.dofs.mesh)
        cell_dofs = Vector{Int}(undef, ndofs_per_cell(op.dofs, k))
        if get_order(op.elements[k]) == 1
          b_local = op.elements[k].b
        elseif get_order(op.elements[k])  < 3
          # In this case: Π^∇ = Π^0
          Π∇s = op.elements[k].G\op.elements[k].B
          b_local = Π∇s'*op.elements[k].b
        else
          throw("The L^2 operator is not yet supported for k>2.")
        end
        # Assembly
        celldofs!(cell_dofs, op.dofs, k)
        assemble!(b,cell_dofs,b_local)
      end
      return b
end

""" Compute stiffness Matrix """
function assemble_stiffnessMat(op::VEMOperators)
    K = create_sparsity_pattern(op.dofs);
    assembler = start_assemble(K)
    for k = 1:getncells(op.dofs.mesh)
      cell_dofs = Vector{Int}(undef, ndofs_per_cell(op.dofs, k))
      # Build local stiffness matrix
      Π∇s = op.elements[k].G\op.elements[k].B
      Π∇ = op.elements[k].D*Π∇s
      nk = size(op.elements[k].G,2)
      G_tilde = [ zeros(1, nk) ;
                  op.elements[k].G[2:end,:] ]
      α = 1;
      K_local = Π∇s'*G_tilde*Π∇s + α*(I-Π∇)'*(I-Π∇)

      # Assembly
      celldofs!(cell_dofs, op.dofs, k)
      assemble!(assembler, cell_dofs, K_local)
    end
    return K
end

""" Compute Mass Matrix """
function assemble_massMat(op::VEMOperators)
  M = create_sparsity_pattern(op.dofs);
  assembler = start_assemble(M)

  for k = 1:getncells(op.dofs.mesh)
    cell_dofs = Vector{Int}(undef, ndofs_per_cell(op.dofs, k))
    if get_order(op.elements[k])  < 3
      # In this case: Π^∇ = Π^0
      Π∇s = op.elements[k].G\op.elements[k].B
      Π∇ = op.elements[k].D*Π∇s

      M_local = Π∇s'*op.elements[k].H*Π∇s
    else
      throw("The L^2 operator is not yet supported for k>2.")
    end
    # Assembly
    celldofs!(cell_dofs, op.dofs, k)
    assemble!(assembler, cell_dofs, M_local)
  end
  return M
end
