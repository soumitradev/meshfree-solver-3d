import "regent"

require "parameter"
require "data_structure"
require "file_properties"

local cstdlib = regentlib.c

__demand(__leaf)
__demand(__inline)
task flux_Gwxp(
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
  var S1 = ut1 * cstdlib.sqrt(beta)
  var S3 = un * cstdlib.sqrt(beta)
  var B1 = 0.5 * cstdlib.exp(- S1 * S1) / cstdlib.sqrt(params.pi * beta)
  var B3 = 0.5 * cstdlib.exp(- S3 * S3) / cstdlib.sqrt(params.pi * beta)
  var A1pos = 0.5 * (1 + cstdlib.erf(S1))
  var A3neg = 0.5 * (1 - cstdlib.erf(S3))

  var u_sqr = ut1*ut1 + ut2*ut2 + un*un
  var temp1 = (ut1*A1pos + B1)

  G[0] = rho * temp1 * A3neg
  G[1] = (pr * A1pos + rho * ut1 * temp1) * A3neg
  G[2] = rho * ut2 * temp1 * A3neg
  G[3] = rho * temp1 * (un * A3neg - B3)

  var temp2 = 2.5 * pr + 0.5 * rho * u_sqr
  var temp3 = (temp2 + pr) * ut1 * A1pos + (temp2 + 0.5 * pr) * B1
  G[4] = temp3 * A3neg - 0.5 * rho * un * B3 * temp1

  return G
end

__demand(__leaf)
__demand(__inline)
task flux_Gwxn(
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
  var S1 = ut1 * cstdlib.sqrt(beta)
  var S3 = un * cstdlib.sqrt(beta)
  var B1 = 0.5 * cstdlib.exp(- S1 * S1) / cstdlib.sqrt(params.pi * beta)
  var B3 = 0.5 * cstdlib.exp(- S3 * S3) / cstdlib.sqrt(params.pi * beta)
  var A1neg = 0.5 * (1 - cstdlib.erf(S1))
  var A3neg = 0.5 * (1 - cstdlib.erf(S3))

  var u_sqr = ut1*ut1 + ut2*ut2 + un*un
  var temp1 = (ut1*A1neg - B1)

  G[0] = rho * temp1 * A3neg
  G[1] = (pr * A1neg + rho * ut1 * temp1) * A3neg
  G[2] = rho * ut2 * temp1 * A3neg
  G[3] = rho * temp1 * (un * A3neg - B3)

  var temp2 = 2.5 * pr + 0.5 * rho * u_sqr
  var temp3 = (temp2 + pr) * ut1 * A1neg - (temp2 + 0.5 * pr) * B1
  G[4] = temp3 * A3neg - 0.5 * rho * un * B3 * temp1

  return G
end

__demand(__leaf)
__demand(__inline)
task flux_Goxp(
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
  var S1 = ut1 * cstdlib.sqrt(beta)
  var S3 = un * cstdlib.sqrt(beta)
  var B1 = 0.5 * cstdlib.exp(- S1 * S1) / cstdlib.sqrt(params.pi * beta)
  var B3 = 0.5 * cstdlib.exp(- S3 * S3) / cstdlib.sqrt(params.pi * beta)
  var A1pos = 0.5 * (1 + cstdlib.erf(S1))
  var A3pos = 0.5 * (1 + cstdlib.erf(S3))

  var u_sqr = ut1*ut1 + ut2*ut2 + un*un
  var temp1 = (ut1*A1pos + B1)

  G[0] = rho * temp1 * A3pos
  G[1] = (pr * A1pos + rho * ut1 * temp1) * A3pos
  G[2] = rho * ut2 * temp1 * A3pos
  G[3] = rho * temp1 * (un * A3pos + B3)

  var temp2 = 2.5 * pr + 0.5 * rho * u_sqr
  var temp3 = (temp2 + pr) * ut1 * A1pos + (temp2 + 0.5 * pr) * B1
  G[4] = temp3 * A3pos + 0.5 * rho * un * B3 * temp1

  return G
end

__demand(__leaf)
__demand(__inline)
task flux_Goxn(
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
  var S1 = ut1 * cstdlib.sqrt(beta)
  var S3 = un * cstdlib.sqrt(beta)
  var B1 = 0.5 * cstdlib.exp(- S1 * S1) / cstdlib.sqrt(params.pi * beta)
  var B3 = 0.5 * cstdlib.exp(- S3 * S3) / cstdlib.sqrt(params.pi * beta)
  var A1neg = 0.5 * (1 - cstdlib.erf(S1))
  var A3pos = 0.5 * (1 + cstdlib.erf(S3))

  var u_sqr = ut1*ut1 + ut2*ut2 + un*un
  var temp1 = (ut1*A1neg - B1)

  G[0] = rho * temp1 * A3pos
  G[1] = (pr * A1neg + rho * ut1 * temp1) * A3pos
  G[2] = rho * ut2 * temp1 * A3pos
  G[3] = rho * temp1 * (un * A3pos + B3)

  var temp2 = 2.5 * pr + 0.5 * rho * u_sqr
  var temp3 = (temp2 + pr) * ut1 * A1neg - (temp2 + 0.5 * pr) * B1
  G[4] = temp3 * A3pos + 0.5 * rho * un * B3 * temp1

  return G
end

__demand(__leaf)
__demand(__inline)
task flux_Gwyp(
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
  var S2 = ut2 * cstdlib.sqrt(beta)
  var S3 = un * cstdlib.sqrt(beta)
  var B2 = 0.5 * cstdlib.exp(- S2 * S2) / cstdlib.sqrt(params.pi * beta)
  var B3 = 0.5 * cstdlib.exp(- S3 * S3) / cstdlib.sqrt(params.pi * beta)
  var A2pos = 0.5 * (1 + cstdlib.erf(S2))
  var A3neg = 0.5 * (1 - cstdlib.erf(S3))

  var u_sqr = ut1*ut1 + ut2*ut2 + un*un
  var temp1 = (ut2*A2pos + B2)

  G[0] = rho * temp1 * A3neg
  G[1] = rho * ut1 * temp1 * A3neg
  G[2] = (pr * A2pos + rho * ut2 * temp1) * A3neg
  G[3] = rho * temp1 * (un * A3neg - B3)

  var temp2 = 2.5 * pr + 0.5 * rho * u_sqr
  var temp3 = (temp2 + pr) * ut2 * A2pos + (temp2 + 0.5 * pr) * B2
  G[4] = temp3 * A3neg - 0.5 * rho * un * B3 * temp1

  return G
end

__demand(__leaf)
__demand(__inline)
task flux_Gwyn(
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
  var S2 = ut2 * cstdlib.sqrt(beta)
  var S3 = un * cstdlib.sqrt(beta)
  var B2 = 0.5 * cstdlib.exp(- S2 * S2) / cstdlib.sqrt(params.pi * beta)
  var B3 = 0.5 * cstdlib.exp(- S3 * S3) / cstdlib.sqrt(params.pi * beta)
  var A2neg = 0.5 * (1 - cstdlib.erf(S2))
  var A3neg = 0.5 * (1 - cstdlib.erf(S3))

  var u_sqr = ut1*ut1 + ut2*ut2 + un*un
  var temp1 = (ut2*A2neg - B2)

  G[0] = rho * temp1 * A3neg
  G[1] = rho * ut1 * temp1 * A3neg
  G[2] = (pr * A2neg + rho * ut2 * temp1) * A3neg
  G[3] = rho * temp1 * (un * A3neg - B3)

  var temp2 = 2.5 * pr + 0.5 * rho * u_sqr
  var temp3 = (temp2 + pr) * ut2 * A2neg - (temp2 + 0.5 * pr) * B2
  G[4] = temp3 * A3neg - 0.5 * rho * un * B3 * temp1

  return G
end

__demand(__leaf)
__demand(__inline)
task flux_Goyp(
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
  var S2 = ut2 * cstdlib.sqrt(beta)
  var S3 = un * cstdlib.sqrt(beta)
  var B2 = 0.5 * cstdlib.exp(- S2 * S2) / cstdlib.sqrt(params.pi * beta)
  var B3 = 0.5 * cstdlib.exp(- S3 * S3) / cstdlib.sqrt(params.pi * beta)
  var A2pos = 0.5 * (1 + cstdlib.erf(S2))
  var A3pos = 0.5 * (1 + cstdlib.erf(S3))

  var u_sqr = ut1*ut1 + ut2*ut2 + un*un
  var temp1 = (ut2*A2pos + B2)

  G[0] = rho * temp1 * A3pos
  G[1] = rho * ut1 * temp1 * A3pos
  G[2] = (pr * A2pos + rho * ut2 * temp1) * A3pos
  G[3] = rho * temp1 * (un * A3pos + B3)

  var temp2 = 2.5 * pr + 0.5 * rho * u_sqr
  var temp3 = (temp2 + pr) * ut2 * A2pos + (temp2 + 0.5 * pr) * B2
  G[4] = temp3 * A3pos + 0.5 * rho * un * B3 * temp1

  return G
end

__demand(__leaf)
__demand(__inline)
task flux_Goyn(
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
  var S2 = ut2 * cstdlib.sqrt(beta)
  var S3 = un * cstdlib.sqrt(beta)
  var B2 = 0.5 * cstdlib.exp(- S2 * S2) / cstdlib.sqrt(params.pi * beta)
  var B3 = 0.5 * cstdlib.exp(- S3 * S3) / cstdlib.sqrt(params.pi * beta)
  var A2neg = 0.5 * (1 - cstdlib.erf(S2))
  var A3pos = 0.5 * (1 + cstdlib.erf(S3))

  var u_sqr = ut1*ut1 + ut2*ut2 + un*un
  var temp1 = (ut2*A2neg - B2)

  G[0] = rho * temp1 * A3pos
  G[1] = rho * ut1 * temp1 * A3pos
  G[2] = (pr * A2neg + rho * ut2 * temp1) * A3pos
  G[3] = rho * temp1 * (un * A3pos + B3)

  var temp2 = 2.5 * pr + 0.5 * rho * u_sqr
  var temp3 = (temp2 + pr) * ut2 * A2neg - (temp2 + 0.5 * pr) * B2
  G[4] = temp3 * A3pos + 0.5 * rho * un * B3 * temp1

  return G
end
