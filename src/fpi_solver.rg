import "regent"

require "parameter"
require "data_structure"
require "file_properties"
require "q_variables"
require "q_derivatives"
require "timestep_delt"
require "flux_residual"
require "state_update"

local cstdlib = regentlib.c
local format = require("std/format")

task fpi_solver(
  params: parameters,
  points: region(ispace(int1d), point),
  file_props: region(ispace(int1d, 1), file_properties),
  t: int
) where reads writes (points, file_props) do
  q_variables(points, file_props)
  format.println("Finished q_variables for {}", t)  
  q_derivatives(params, points, file_props)
  format.println("Finished q_derivatives for {}", t)  
  timestep_delt(params, points)
  flux_residual(params, points)
  state_update(params, points, file_props)

  if (t <= 2) then
    file_props[0].res_old = file_props[0].res_new
    file_props[0].residue = 0
  else
    file_props[0].residue = cstdlib.log10(file_props[0].res_new / file_props[0].res_old)
  end
end
