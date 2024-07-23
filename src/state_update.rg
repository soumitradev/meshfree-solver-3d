import "regent"

require "data_structure"
require "file_properties"
require "conserved_to_primitive"
require "primitive_to_conserved"
require "conserved_vector_ubar"
require "file_properties"

local cstdlib = regentlib.c

task state_update(
  params: parameters,
  points: region(ispace(int1d), point),
  file_props: region(ispace(int1d, 1), file_properties)
) where
reads (
  points.{status, conn, nbhs, x, y, z, flux_res, tan1, tan2, nor},
  file_props.{max_points}
),
reads writes (points.{prim}),
writes (points.{delt}, file_props.{res_new})
do
  var sum_res_sqr: double = 0
  var max_res: double = 0
  var U: double[5]
  var min_point = 0

  for i in points do
    if (points[i].status == 0) then
      U = primitive_to_conserved(points, i)
      var temp = U[0]
      for p = 0, 5 do
        U[p] = U[p] - points[i].flux_res[p]
      end
      U[3] = 0

      var res_sqr = (U[0] - temp)*(U[0] - temp)
      if (res_sqr > max_res) then
        max_res = res_sqr
      end
      sum_res_sqr = sum_res_sqr + res_sqr
      conserved_to_primitive(points, i, U)
    elseif (points[i].status == 1) then
      U = conserved_vector_ubar(params, points, i)
      for p = 0, 5 do
        U[p] = U[p] - points[i].flux_res[p]
      end
      conserved_to_primitive(points, i, U)
    elseif (points[i].status == 2) then
      U = primitive_to_conserved(points, i)
      var temp = U[0]
      for p = 0, 5 do
        U[p] = U[p] - points[i].flux_res[p]
      end

      var res_sqr = (U[0] - temp)*(U[0] - temp)
      if (res_sqr > max_res) then
        max_res = res_sqr
      end
      sum_res_sqr = sum_res_sqr + res_sqr
      conserved_to_primitive(points, i, U)
    elseif (points[i].status == 6) then
      -- Outlet points
      -- NOTE: The reason here we are computing outlet points before inlet points is
      -- because these two parts share a common variable min_point, and we want to
      -- preserve this exact order of execution that was in the serial code just to be
      -- safe.
      var min_dist: double = 100000
      for j = 0, points[i].nbhs do
        var nbh = points[i].conn[j]
        var dx = points[nbh].x - points[i].x
        var dy = points[nbh].y - points[i].y
        var dz = points[nbh].z - points[i].z
        var ds = cstdlib.sqrt(dx * dx + dy * dy + dz * dz)
        if (ds < min_dist and points[nbh].status ~= 1) then
          min_dist = ds
          min_point = nbh
        end
      end
      for p = 0, 5 do
        points[i].prim[p] = points[min_point].prim[p]
      end
    elseif (points[i].status == 5) then
      -- Inlet points
      var min_dist: double = 100000
      for j = 0, points[i].nbhs do
        var nbh = points[i].conn[j]
        var dx = points[nbh].x - points[i].x
        var dy = points[nbh].y - points[i].y
        var dz = points[nbh].z - points[i].z
        var ds = cstdlib.sqrt(dx * dx + dy * dy + dz * dz)
        if (ds < min_dist and points[nbh].status ~= 1) then
          min_dist = ds
          min_point = nbh
        end
      end
    end
  end

  file_props[0].res_new = cstdlib.sqrt(sum_res_sqr) / file_props[0].max_points
end
