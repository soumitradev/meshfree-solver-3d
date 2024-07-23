import "regent"

require "data_structure"
require "file_properties"

local cstdlib = regentlib.c

task primitive_to_conserved(
  points: region(ispace(int1d), point),
  i: int
) where reads (points.{prim, tan1, tan2, nor}) do
  var U: double[5]
  U[0] = points[i].prim[0]
  var temp1 = points[i].prim[0]*points[i].prim[1]
  var temp2 = points[i].prim[0]*points[i].prim[2]
  var temp3 = points[i].prim[0]*points[i].prim[3]

  U[1] = temp1*points[i].tan1[0] + temp2*points[i].tan1[1] + temp3*points[i].tan1[2]
  U[2] = temp1*points[i].tan2[0] + temp2*points[i].tan2[1] + temp3*points[i].tan2[2]
  U[3] = temp1*points[i].nor[0] + temp2*points[i].nor[1] + temp3*points[i].nor[2]

  temp1 = points[i].prim[1] * points[i].prim[1]
    + points[i].prim[2] * points[i].prim[2]
    + points[i].prim[3] * points[i].prim[3]

  U[4] = 2.5 * points[i].prim[4] + 0.5 * U[0] * temp1
  return U
end
