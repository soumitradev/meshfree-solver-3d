import "regent"

require "data_structure"
require "file_properties"

local cstdlib = regentlib.c

task q_variables(
  points: region(ispace(int1d), point),
  file_props: region(ispace(int1d, 1), file_properties)
) where reads (points.{prim}, file_props.{max_points}), writes (points.{q}) do
  var rho: double = 0
  var u1: double = 0
  var u2: double = 0
  var u3: double = 0
  var pr: double = 0
  var beta: double = 0
  var two_times_beta: double = 0

  for i = 0, file_props[0].max_points do
    rho = points[i].prim[0]
    u1 = points[i].prim[1]
    u2 = points[i].prim[2]
    u3 = points[i].prim[3]
    pr = points[i].prim[4]
    beta = 0.5 * rho / pr
    two_times_beta = 2.0*beta

    points[i].q[0] =
      cstdlib.log(rho) + (cstdlib.log(beta)*2.5) - beta*(u1*u1 + u2*u2 + u3*u3)
    points[i].q[1] = two_times_beta*u1
    points[i].q[2] = two_times_beta*u2
    points[i].q[3] = two_times_beta*u3
    points[i].q[4] = -two_times_beta
  end
end
