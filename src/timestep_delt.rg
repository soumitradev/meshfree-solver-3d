import "regent"

require "parameter"
require "data_structure"
require "file_properties"

local cstdlib = regentlib.c

task timestep_delt(
  params: parameters,
  points: region(ispace(int1d), point)
) where
reads (points.{prim, conn, nbhs, x, y, z}),
writes (points.{delt})
do
  var rho: double = 0
  var u1: double = 0
  var u2: double = 0
  var u3: double = 0
  var pr: double = 0

  var mod_u: double = 0
  var delt: double = 0
  var min_delt: double = 0

  var delx: double = 0
  var dely: double = 0
  var delz: double = 0
  var dist: double = 0

  for i in points do
    min_delt = 1
    for k = 0, points[i].nbhs do
      var j = points[i].conn[k]

      rho = points[j].prim[0]
      u1 = points[j].prim[1]
      u2 = points[j].prim[2]
      u3 = points[j].prim[3]
      pr = points[j].prim[4]
      
      delx = points[j].x - points[i].x
      dely = points[j].y - points[i].y
      delz = points[j].z - points[i].z

      dist = cstdlib.sqrt(delx*delx + dely*dely + delz*delz)
  
      mod_u = cstdlib.sqrt(u1*u1 + u2*u2 + u3*u3)
      delt = (params.CFL * dist) / (mod_u + (3 * cstdlib.sqrt(pr / rho)))

      if (delt < min_delt) then
        min_delt = delt
      end
    end
    points[i].delt = min_delt
  end
end
