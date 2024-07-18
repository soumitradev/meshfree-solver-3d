import "regent"

require "parameter"
require "data_structure"
require "file_properties"

local cstdlib = regentlib.c

task q_derivatives(
  params: parameters,
  points: region(ispace(int1d), point),
  file_props: region(ispace(int1d, 1), file_properties)
) where
reads (points.{x, y, z, q, nbhs, conn}, file_props.{max_points}),
reads writes (points.{qm, dq})
do
  if (params.second_order_flag == 0) then
    for i = 0, file_props[0].max_points do
      __demand(__vectorize)
      for p = 0, 5 do
        points[i].qm[p] = points[i].q[p]
        points[i].qm[p + 5] = points[i].q[p]
      end
      
      for j = 0, points[i].nbhs do
        var nbh = points[i].conn[j]
        for p = 0, 5 do
          if (points[nbh].q[p] > points[i].qm[p]) then
            points[i].qm[p] = points[nbh].q[p]
          end
          if (points[nbh].q[p] < points[i].qm[p + 5]) then
            points[i].qm[p + 5] = points[nbh].q[p]
          end
        end
      end

      __demand(__vectorize)
      for p = 0, 5 do
        points[i].dq[p] = 0
        points[i].dq[p + 5] = 0
        points[i].dq[p + 10] = 0
      end
    end
  elseif (params.second_order_flag == 1) then
    for i = 0, file_props[0].max_points do
      var sum_delx_sqr: double = 0
      var sum_dely_sqr: double = 0
      var sum_delz_sqr: double = 0
      var sum_delx_dely: double = 0
      var sum_dely_delz: double = 0
      var sum_delz_delx: double = 0
      var sum_delx_delq: double[5]
      var sum_dely_delq: double[5]
      var sum_delz_delq: double[5]
      var temp: double

      var x_i = points[i].x
      var y_i = points[i].y
      var z_i = points[i].z

      __demand(__vectorize)
      for p = 0, 5 do
        points[i].qm[p] = points[i].q[p]
        points[i].qm[p + 5] = points[i].q[p]
      end

      for j = 0, points[i].nbhs do
        var nbh = points[i].conn[j]
        for p = 0, 5 do
          if (points[nbh].q[p] > points[i].qm[p]) then
            points[i].qm[p] = points[nbh].q[p]
          end
          if (points[nbh].q[p] < points[i].qm[p + 5]) then
            points[i].qm[p + 5] = points[nbh].q[p]
          end
        end
        
        var x_j = points[nbh].x
        var y_j = points[nbh].y
        var z_j = points[nbh].z

        var delx = x_j - x_i
        var dely = y_j - y_i
        var delz = z_j - z_i

        var dist = cstdlib.sqrt(delx*delx + dely*dely + delz*delz)
        var weights = 1.0/(cstdlib.pow(dist, params.power))

        sum_delx_sqr = sum_delx_sqr + delx*delx*weights
        sum_dely_sqr = sum_dely_sqr + dely*dely*weights
        sum_delz_sqr = sum_delz_sqr + delz*delz*weights

        sum_delx_dely = sum_delx_dely + delx*dely*weights
        sum_dely_delz = sum_dely_delz + dely*delz*weights
        sum_delz_delx = sum_delz_delx + delz*delx*weights

        __demand(__vectorize)
        for p = 0, 5 do
          temp = points[nbh].q[p] - points[i].q[p]
          sum_delx_delq[p] = sum_delx_delq[p] + weights*delx*temp
          sum_dely_delq[p] = sum_dely_delq[p] + weights*dely*temp
          sum_delz_delq[p] = sum_delz_delq[p] + weights*delz*temp
        end
      end

      var det = sum_delx_sqr*(sum_dely_sqr*sum_delz_sqr - sum_dely_delz*sum_dely_delz)
        - sum_delx_dely*(sum_delx_dely*sum_delz_sqr - sum_dely_delz*sum_delz_delx)
        + sum_delz_delx*(sum_delx_dely*sum_dely_delz - sum_dely_sqr*sum_delz_delx)
      var one_by_det = 1.0 / det

      __demand(__vectorize)
      for p = 0, 5 do
        temp =
          sum_delx_delq[p]*(sum_dely_sqr*sum_delz_sqr - sum_dely_delz*sum_dely_delz)
          - sum_dely_delq[p]*(sum_delx_dely*sum_delz_sqr - sum_delz_delx*sum_dely_delz)
          + sum_delz_delq[p]*(sum_delx_dely*sum_dely_delz - sum_delz_delx*sum_dely_sqr)
        points[i].dq[p] = temp*one_by_det
        temp =
          sum_delx_sqr*(sum_dely_delq[p]*sum_delz_sqr - sum_dely_delz*sum_delz_delq[p])
          - sum_delx_dely*(sum_delx_delq[p]*sum_delz_sqr - sum_delz_delx*sum_delz_delq[p])
          + sum_delz_delx*(sum_delx_delq[p]*sum_dely_delz - sum_delz_delx*sum_dely_delq[p])	
        points[i].dq[p + 5] = temp*one_by_det
        temp =
          sum_delx_sqr*(sum_dely_sqr*sum_delz_delq[p] - sum_dely_delq[p]*sum_dely_delz)
          - sum_delx_dely*(sum_delx_dely*sum_delz_delq[p] - sum_delx_delq[p]*sum_dely_delz)
          + sum_delz_delx*(sum_delx_dely*sum_dely_delq[p] - sum_delx_delq[p]*sum_dely_sqr)
        points[i].dq[p + 10] = temp*one_by_det
      end
    end

    var dq_temp = region(ispace(int1d, file_props[0].max_points), double[3][5])
    for p = 0, params.inner_iterations do
      for i = 0, file_props[0].max_points do
        var sum_delx_sqr: double = 0
        var sum_dely_sqr: double = 0
        var sum_delz_sqr: double = 0
        var sum_delx_dely: double = 0
        var sum_dely_delz: double = 0
        var sum_delz_delx: double = 0
        var sum_delx_delq: double[5]
        var sum_dely_delq: double[5]
        var sum_delz_delq: double[5]
        var qtilde_i: double
        var qtilde_nbh: double
        var temp: double

        var x_i = points[i].x
        var y_i = points[i].y
        var z_i = points[i].z

        for j = 0, points[i].nbhs do
          var nbh = points[i].conn[j]
          
          var x_j = points[nbh].x
          var y_j = points[nbh].y
          var z_j = points[nbh].z

          var delx = x_j - x_i
          var dely = y_j - y_i
          var delz = z_j - z_i

          var dist = cstdlib.sqrt(delx*delx + dely*dely + delz*delz)
          var weights = 1.0/(cstdlib.pow(dist, params.power))

          sum_delx_sqr = sum_delx_sqr + delx*delx*weights
          sum_dely_sqr = sum_dely_sqr + dely*dely*weights
          sum_delz_sqr = sum_delz_sqr + delz*delz*weights

          sum_delx_dely = sum_delx_dely + delx*dely*weights
          sum_dely_delz = sum_dely_delz + dely*delz*weights
          sum_delz_delx = sum_delz_delx + delz*delx*weights

          __demand(__vectorize)
          for p = 0, 5 do
            qtilde_i = points[i].q[p] - 0.5*(delx*points[i].dq[p] + dely*points[i].dq[p + 5] + delz*points[i].dq[p + 10])
            qtilde_nbh = points[nbh].q[p] - 0.5*(delx*points[nbh].dq[p] + dely*points[nbh].dq[p + 5] + delz*points[nbh].dq[p + 10])
            temp = qtilde_nbh - qtilde_i
            sum_delx_delq[p] = sum_delx_delq[p] + weights*delx*temp
            sum_dely_delq[p] = sum_dely_delq[p] + weights*dely*temp
            sum_delz_delq[p] = sum_delz_delq[p] + weights*delz*temp
          end
        end

        var det = sum_delx_sqr*(sum_dely_sqr*sum_delz_sqr - sum_dely_delz*sum_dely_delz)
          - sum_delx_dely*(sum_delx_dely*sum_delz_sqr - sum_dely_delz*sum_delz_delx)
          + sum_delz_delx*(sum_delx_dely*sum_dely_delz - sum_dely_sqr*sum_delz_delx)
        var one_by_det = 1.0 / det

        __demand(__vectorize)
        for p = 0, 5 do
          temp =
            sum_delx_delq[p]*(sum_dely_sqr*sum_delz_sqr - sum_dely_delz*sum_dely_delz)
            - sum_dely_delq[p]*(sum_delx_dely*sum_delz_sqr - sum_delz_delx*sum_dely_delz)
            + sum_delz_delq[p]*(sum_delx_dely*sum_dely_delz - sum_delz_delx*sum_dely_sqr)
          dq_temp[i][0][p] = temp*one_by_det
          temp =
            sum_delx_sqr*(sum_dely_delq[p]*sum_delz_sqr - sum_dely_delz*sum_delz_delq[p])
            - sum_delx_dely*(sum_delx_delq[p]*sum_delz_sqr - sum_delz_delx*sum_delz_delq[p])
            + sum_delz_delx*(sum_delx_delq[p]*sum_dely_delz - sum_delz_delx*sum_dely_delq[p])	
          dq_temp[i][1][p] = temp*one_by_det
          temp =
            sum_delx_sqr*(sum_dely_sqr*sum_delz_delq[p] - sum_dely_delq[p]*sum_dely_delz)
            - sum_delx_dely*(sum_delx_dely*sum_delz_delq[p] - sum_delx_delq[p]*sum_dely_delz)
            + sum_delz_delx*(sum_delx_dely*sum_dely_delq[p] - sum_delx_delq[p]*sum_dely_sqr)
          dq_temp[i][2][p] = temp*one_by_det
        end
      end
      for i = 0, file_props[0].max_points do
        __demand(__vectorize)
        for j = 0, 5 do
          __demand(__vectorize)
          for k = 0, 3 do
            points[i].dq[j + k * 5] = dq_temp[i][k][j]
          end
        end
      end
    end
  end
end
