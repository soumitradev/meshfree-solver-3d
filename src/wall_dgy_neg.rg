import "regent"

require "parameter"
require "data_structure"
require "file_properties"
require "qtilde_to_primitive"
require "octant_fluxes"
require "limiters"

local cstdlib = regentlib.c

task wall_dgy_neg(
  params: parameters,
  points: region(ispace(int1d), point),
  i: int
): double[5] where
reads (points.{x, y, z, tan1, tan2, nor, q, dq, qm, min_dist, yneg_nbhs, yneg_conn})
do
  var sum_delx_sqr: double = 0
  var sum_dely_sqr: double = 0
  var sum_delz_sqr: double = 0
  var sum_delx_dely: double = 0
  var sum_dely_delz: double = 0
  var sum_delz_delx: double = 0
  var sum_delx_delf: double[5]
  var sum_dely_delf: double[5]
  var sum_delz_delf: double[5]
  
  var G_i: double[5]
  var G_j: double[5]
  var G: double[5]

  var temp: double = 0
  var dist: double = 0
  var weights: double = 0
  var det: double = 0

  var tan1: double[3]
  var tan2: double[3]
  var nor: double[3]

  var phi: double[5]
  var phi_i: double[5]
  var phi_j: double[5]

  var qtilde: double[5]
  var prim: double[5]
  
  var delx: double = 0
  var dely: double = 0
  var delz: double = 0

  var dels: double = 0
  var delt: double = 0
  var deln: double = 0

  var dels_weights: double = 0
  var delt_weights: double = 0
  var deln_weights: double = 0

  var x_j: double = 0
  var y_j: double = 0
  var z_j: double = 0

  var x_i = points[i].x
  var y_i = points[i].y
  var z_i = points[i].z

  for p = 0, 3 do
    tan1[p] = points[i].tan1[p]
    tan2[p] = points[i].tan2[p]
    nor[p] = points[i].nor[p]
  end

  for p = 0, 5 do
    phi[p] = 1
  end

  if (params.limiter_flag == 1) then
    for k = 0, points[i].yneg_nbhs do
      var j = points[i].yneg_conn[k]

      x_j = points[j].x
      y_j = points[j].y
      z_j = points[j].z

      delx = x_j - x_i
      dely = y_j - y_i
      delz = z_j - z_i

      for p = 0, 5 do
        temp = delx*points[i].dq[p] + dely*points[i].dq[p + 5] + delz*points[i].dq[p + 10]
        qtilde[p] = points[i].q[p] - 0.5*temp
      end
      -- NOTE: There is probably a better way of doing this by passing this by reference.
      -- For now, we are relying on inlining this function, so it doesn't matter.
      phi_i = venkat_limiter(params, points, qtilde, i)

      for p = 0, 5 do
        temp = delx*points[j].dq[p] + dely*points[j].dq[p +5] + delz*points[j].dq[p + 10]
        qtilde[p] = points[j].q[p] - 0.5*temp
      end
      -- NOTE: There is probably a better way of doing this by passing this by reference.
      -- For now, we are relying on inlining this function, so it doesn't matter.
      phi_j = venkat_limiter(params, points, qtilde, j)
      
      for p = 0, 5 do
        if (phi_j[p] < phi[p]) then
          phi[p] = phi_j[p]
        end
        if (phi_i[p] < phi[p]) then
          phi[p] = phi_i[p]
        end
      end
    end
  end

  for k = 0, points[i].yneg_nbhs do
    var j = points[i].yneg_conn[k]

    x_j = points[j].x
    y_j = points[j].y
    z_j = points[j].z

    delx = x_j - x_i
    dely = y_j - y_i
    delz = z_j - z_i

    dels = delx*tan1[0] + dely*tan1[1] + delz*tan1[2]
    delt = delx*tan2[0] + dely*tan2[1] + delz*tan2[2]
    deln = delx*nor[0] + dely*nor[1] + delz*nor[2]

    dist = cstdlib.sqrt(dels*dels + delt*delt + deln*deln)
    weights = 1 / (cstdlib.pow(dist, params.power))

    dels_weights = dels*weights
    delt_weights = delt*weights
    deln_weights = deln*weights

    sum_delx_sqr = sum_delx_sqr + dels*dels_weights
    sum_dely_sqr = sum_dely_sqr + delt*delt_weights
    sum_delz_sqr = sum_delz_sqr + deln*deln_weights

    sum_delx_dely = sum_delx_dely + dels*delt_weights
    sum_dely_delz = sum_dely_delz + delt*deln_weights
    sum_delz_delx = sum_delz_delx + deln*dels_weights

    for p = 0, 5 do
      temp = delx*points[i].dq[p] + dely*points[i].dq[p + 5] + delz*points[i].dq[p + 10]
      qtilde[p] = points[i].q[p] - 0.5*phi[p]*temp
    end
    prim = qtilde_to_primitive(qtilde)
    G_i = flux_Gwyn(params, tan1, tan2, nor, prim)

    for p = 0, 5 do
      temp = delx*points[j].dq[p] + dely*points[j].dq[p + 5] + delz*points[j].dq[p + 10]
      qtilde[p] = points[j].q[p] - 0.5*phi[p]*temp
    end
    prim = qtilde_to_primitive(qtilde)
    G_j = flux_Gwyn(params, tan1, tan2, nor, prim)

    for p = 0, 5 do
      temp = G_j[p] - G_i[p]
      sum_delx_delf[p] = sum_delx_delf[p] + temp*dels_weights
      sum_dely_delf[p] = sum_dely_delf[p] + temp*delt_weights
      sum_delz_delf[p] = sum_delz_delf[p] + temp*deln_weights
    end
  end

  det = sum_delx_sqr*(sum_dely_sqr*sum_delz_sqr - sum_dely_delz*sum_dely_delz)
    - sum_delx_dely*(sum_delx_dely*sum_delz_sqr - sum_dely_delz*sum_delz_delx)
    + sum_delz_delx*(sum_delx_dely*sum_dely_delz - sum_dely_sqr*sum_delz_delx)

  for p = 0, 5 do
    temp = sum_delx_delf[p]*(sum_dely_sqr*sum_delz_sqr - sum_dely_delz*sum_dely_delz)
      - sum_dely_delf[p]*(sum_delx_dely*sum_delz_sqr - sum_delz_delx*sum_dely_delz)
      + sum_delz_delf[p]*(sum_delx_dely*sum_dely_delz - sum_delz_delx*sum_dely_sqr)
    G[p] = temp / det
  end

  return G
end
