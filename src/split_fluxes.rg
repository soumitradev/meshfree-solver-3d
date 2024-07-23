import "regent"

require "parameter"
require "data_structure"

local cstdlib = regentlib.c

__demand(__leaf)
__demand(__inline)
task flux_Gzn(
  params: parameters,
  tan1: double[3],
  tan2: double[3],
  nor: double[3],
  prim: double[5]
): double[5]
  var G: double[5]

  var rho = prim[0]
  var u1 = prim[1]
  var u2 = prim[2]
  var u3 = prim[3]
  var pr = prim[4]

  var ut1 = tan1[0]*u1 + tan1[1]*u2 + tan1[2]*u3
  var ut2 = tan2[0]*u1 + tan2[1]*u2 + tan2[2]*u3
  var un = nor[0]*u1 + nor[1]*u2 + nor[2]*u3

  var beta = 0.5 * rho / pr
  var S3 = un * cstdlib.sqrt(beta)
  var B3 = 0.5 * cstdlib.exp(- S3 * S3) / cstdlib.sqrt(params.pi * beta)
  var A3neg = 0.5 * (1 - cstdlib.erf(S3))

  var u_sqr = ut1*ut1 + ut2*ut2 + un*un
  var temp1 = (un*A3neg - B3)

  G[0] = rho * temp1
  G[1] = rho * ut1 * temp1
  G[2] = rho * ut2 * temp1
  G[3] = pr * A3neg + rho * un * temp1

  temp1 = 2.5 * pr + 0.5 * rho * u_sqr
  G[4] = (temp1 + pr) * un * A3neg - (temp1 + 0.5 * pr) * B3

  return G
end
