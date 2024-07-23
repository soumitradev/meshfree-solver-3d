import "regent"

require "parameter"
require "data_structure"
require "file_properties"

local cstdlib = regentlib.c

__demand(__leaf)
__demand(__inline)
task qtilde_to_primitive(q: double[5]): double[5]
  var prim: double[5]

  var beta = -q[4] * 0.5
  var temp = 0.5 / beta
  prim[1] = q[1]*temp
  prim[2] = q[2]*temp
  prim[3] = q[3]*temp

  var temp1 = q[0] + beta*(prim[1]*prim[1] + prim[2]*prim[2] + prim[3]*prim[3])
  var temp2 = temp1 - (cstdlib.log(beta)*2.5)

  prim[0] = cstdlib.exp(temp2)
  prim[4] = prim[0]*temp

  return prim
end
