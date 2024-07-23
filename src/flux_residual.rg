import "regent"

require "parameter"
require "data_structure"
require "file_properties"
require "wall_dgx_pos"
require "wall_dgx_neg"
require "wall_dgy_pos"
require "wall_dgy_neg"
require "wall_dgz_neg"
require "outer_dgx_pos"
require "outer_dgx_neg"
require "outer_dgy_pos"
require "outer_dgy_neg"
require "outer_dgz_pos"
require "interior_dgx_pos"
require "interior_dgx_neg"
require "interior_dgy_pos"
require "interior_dgy_neg"
require "interior_dgz_pos"
require "interior_dgz_neg"

local cstdlib = regentlib.c

task flux_residual(
  params: parameters,
  points: region(ispace(int1d), point)
) where
reads (points.{
  x, y, z,
  tan1, tan2, nor,
  q, dq, qm,
  min_dist, delt,
  xpos_nbhs, xpos_conn,
  xneg_nbhs, xneg_conn,
  ypos_nbhs, ypos_conn,
  yneg_nbhs, yneg_conn,
  zpos_nbhs, zpos_conn,
  zneg_nbhs, zneg_conn,
  nbhs, conn,
  status
}),
writes (points.flux_res)
do
  var Gxp: double[5]
  var Gxn: double[5]
  var Gyp: double[5]
  var Gyn: double[5]
  var Gzn: double[5]
  var Gzp: double[5]

  for i in points do
    if (points[i].status == 0) then
      Gxp = wall_dgx_pos(params, points, i)
      Gxn = wall_dgx_neg(params, points, i)
      Gyp = wall_dgy_pos(params, points, i)
      Gyn = wall_dgy_neg(params, points, i)
      Gzn = wall_dgz_neg(params, points, i)

      for p = 0, 5 do
        points[i].flux_res[p] =
          2 * points[i].delt * (Gxp[p] + Gxn[p] + Gyp[p] + Gyn[p] + Gzn[p])
      end
    elseif (points[i].status == 1) then
      Gxp = outer_dgx_pos(params, points, i)
      Gxn = outer_dgx_neg(params, points, i)
      Gyp = outer_dgy_pos(params, points, i)
      Gyn = outer_dgy_neg(params, points, i)
      Gzp = outer_dgz_pos(params, points, i)

      for p = 0, 5 do
        points[i].flux_res[p] =
          points[i].delt * (Gxp[p] + Gxn[p] + Gyp[p] + Gyn[p] + Gzp[p])
      end
    elseif (points[i].status == 2) then
      Gxp = interior_dgx_pos(params, points, i)
      Gxn = interior_dgx_neg(params, points, i)
      Gyp = interior_dgy_pos(params, points, i)
      Gyn = interior_dgy_neg(params, points, i)
      Gzn = interior_dgz_neg(params, points, i)
      Gzp = interior_dgz_pos(params, points, i)

      for p = 0, 5 do
        points[i].flux_res[p] =
          points[i].delt * (Gxp[p] + Gxn[p] + Gyp[p] + Gyn[p] + Gzp[p] + Gzn[p])
      end
    end
  end
end
