import "regent"

require "file_properties"
require "parameter"
require "data_structure"
require "point_preprocessor"
require "initial_conditions"
require "q_lskum"
require "print_output"

local format = require("std/format")

task main()
  var params = init_params()
  -- Define file properties as a region as of now so we can pass it down to tasks
  var file_props = region(ispace(int1d, 1), file_properties)
  file_props[0].max_points = 580485
  var points = region(ispace(int1d, file_props[0].max_points), point)

  format.println("Reading grid from file")
  read_input_point_data(points, file_props)
  format.println("Finished reading grid from file")
  
  format.println("Setting initial conditions")
  initial_conditions(params, points, file_props)
  format.println("Initial conditions set")
  
  format.println("Starting solver")
  q_lskum(params, points, file_props)
  format.println("Finished running solver for {} iterations", params.max_iters)

  format.println("Writing output file and restart file")
  print_output(points)
  format.println("Finished writing output and restart files")
end

regentlib.saveobj(main, '/build/meshfree_solver_test.out', 'executable', nil, {'-lm'})
