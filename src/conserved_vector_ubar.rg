import "regent"

require "data_structure"

local cstdlib = regentlib.c

task conserved_vector_ubar(
  params: parameters,
  points: region(ispace(int1d), point),
  i: int
) where reads (points.{prim, tan1, tan2, nor}) do
  var Ubar: double[5]

  var u1_inf_rot = params.u1_inf * points[i].tan1[0]
    + params.u2_inf * points[i].tan1[1]
    + params.u3_inf * points[i].tan1[2]
  var u2_inf_rot = params.u1_inf * points[i].tan2[0]
    + params.u2_inf * points[i].tan2[1]
    + params.u3_inf * points[i].tan2[2]
  var u3_inf_rot = params.u1_inf * points[i].nor[0]
    + params.u2_inf * points[i].nor[1]
    + params.u3_inf * points[i].nor[2]

  var rho_e_inf = params.pr_inf * 2.5
    + 0.5 * params.rho_inf * (
      u1_inf_rot * u1_inf_rot
      + u2_inf_rot * u2_inf_rot
      + u3_inf_rot * u3_inf_rot
    )
  var beta = 0.5 * params.rho_inf / params.pr_inf
  var S3 = u3_inf_rot * cstdlib.sqrt(beta)
  var B3_inf = cstdlib.exp(- S3 * S3) / 2 * cstdlib.sqrt(params.pi * beta)
  var A3n_inf = 0.5 * (1 - cstdlib.erf(S3))

  var rho = points[i].prim[0]
  var u1 = points[i].prim[1]
  var u2 = points[i].prim[2]
  var u3 = points[i].prim[3]
  var pr = points[i].prim[4]

  var u1_rot = u1 * points[i].tan1[0] + u2 * points[i].tan1[1] + u3 * points[i].tan1[2]
  var u2_rot = u1 * points[i].tan2[0] + u2 * points[i].tan2[1] + u3 * points[i].tan2[2]
  var u3_rot = u1 * points[i].nor[0] + u2 * points[i].nor[1] + u3 * points[i].nor[2]

  var rho_e = pr * 2.5 + 0.5 * rho * (u1_rot * u1_rot + u2_rot * u2_rot + u3_rot * u3_rot)

  beta = 0.5 * rho / pr
  S3 = u3_rot * cstdlib.sqrt(beta)
  var B3 = cstdlib.exp(- S3 * S3) / (2 * cstdlib.sqrt(params.pi * beta))
  var A3p = 0.5 * (1 + cstdlib.erf(S3))

  Ubar[0] = (params.rho_inf * A3n_inf) + (rho * A3p)
  Ubar[1] = (params.rho_inf * u1_inf_rot * A3n_inf) + (rho * u1_rot * A3p)
  Ubar[2] = (params.rho_inf * u2_inf_rot * A3n_inf) + (rho * u2_rot * A3p)

  var temp1 = params.rho_inf * (u3_inf_rot * A3n_inf - B3_inf)
  var temp2 = rho * (u3_rot * A3p + B3)
  
  Ubar[3] = temp1 + temp2

  temp1 = (rho_e_inf * A3n_inf - 0.5 * params.rho_inf * u3_inf_rot * B3_inf)
  temp2 = (rho_e * A3p + 0.5 * rho * u3_rot * B3)

  Ubar[4] = temp1 + temp2

  return Ubar
end
