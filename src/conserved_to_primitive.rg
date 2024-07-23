import "regent"

require "data_structure"

local cstdlib = regentlib.c

task conserved_to_primitive(
  points: region(ispace(int1d), point),
  i: int,
  U: double[5]
) where reads (points.{tan1, tan2, nor}), reads writes (points.{prim}) do
  points[i].prim[0] = U[0]
  var temp = 1 / U[0]
  
  var U1_rot = U[1]
  var U2_rot = U[2]
  var U3_rot = U[3]

  U[1] =
    U1_rot * points[i].tan1[0] + U2_rot * points[i].tan2[0] + U3_rot * points[i].nor[0]
  U[2] =
    U1_rot * points[i].tan1[1] + U2_rot * points[i].tan2[1] + U3_rot * points[i].nor[1]
  U[3] =
    U1_rot * points[i].tan1[2] + U2_rot * points[i].tan2[2] + U3_rot * points[i].nor[2]

  points[i].prim[1] = temp * U[1]
  points[i].prim[2] = temp * U[2]
  points[i].prim[3] = temp * U[3]

  temp = points[i].prim[1] * points[i].prim[1]
    + points[i].prim[2] * points[i].prim[2]
    + points[i].prim[3] * points[i].prim[3]
  points[i].prim[4] = 0.4 * (U[4] - 0.5 * points[i].prim[0] * temp)
end
