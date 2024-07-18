import "regent"

require "parameter"
require "data_structure"
require "file_properties"

local cstdlib = regentlib.c

task fpi_solver(
  params: parameters,
  points: region(ispace(int1d), point),
  file_props: region(ispace(int1d, 1), file_properties),
  t: int
) where reads writes (points, file_props) do
  -- q_variables()
  -- q_derivatives()
  -- timestep_delt()
  -- flux_residual()
  -- state_update()

  if (t <= 2) then
    file_props[0].res_old = file_props[0].res_new
    file_props[0].residue = 0
  else
    file_props[0].residue = cstdlib.log10(file_props[0].res_new / file_props[0].res_old)
  end
end
