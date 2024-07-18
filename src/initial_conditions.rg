import "regent"

require "file_properties"
require "data_structure"
require "parameter"

local cstdlib = regentlib.c

__demand(__inline)
task initial_conditions(
  params: parameters,
  points: region(ispace(int1d), point),
  file_props: region(ispace(int1d, 1), file_properties)
) where reads (file_props.{max_points}), reads writes (points.{prim}) do
  if (params.initial_conditions_flag == 0) then
    for i = 0, file_props[0].max_points do
      points[i].prim[0] = params.rho_inf
      points[i].prim[1] = params.u1_inf
      points[i].prim[2] = params.u2_inf
      points[i].prim[3] = params.u3_inf
      points[i].prim[4] = params.pr_inf
    end
  elseif (params.initial_conditions_flag == 1) then
    var restart = cstdlib.fopen("/data/restart_file", "r")
    var temp: regentlib.array(double, 5)
    for i = 0, file_props[0].max_points do
      cstdlib.fscanf(
        restart,
        "%lf %lf %lf %lf %lf",
        &temp[0],
        &temp[1],
        &temp[2],
        &temp[3],
        &temp[4]
      )
      points[i].prim = temp
    end
    cstdlib.fclose(restart)
  end
end
