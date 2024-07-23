import "regent"

require "parameter"
require "data_structure"
require "file_properties"
require "q_variables"
require "q_derivatives"
require "timestep_delt"
require "flux_residual"

local cstdlib = regentlib.c

task fpi_solver(
  params: parameters,
  points: region(ispace(int1d), point),
  file_props: region(ispace(int1d, 1), file_properties),
  t: int
) where reads writes (points, file_props) do
  q_variables(points, file_props)
  q_derivatives(params, points, file_props)
  timestep_delt(params, points)
  flux_residual(params, points)
  -- state_update()

  if (t <= 2) then
    file_props[0].res_old = file_props[0].res_new
    file_props[0].residue = 0
  else
    file_props[0].residue = cstdlib.log10(file_props[0].res_new / file_props[0].res_old)
  end
end
