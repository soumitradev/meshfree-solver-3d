import "regent"

require "file_properties"
require "parameter"
require "data_structure"
require "point_preprocessor"
require "initial_conditions"

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

end

regentlib.saveobj(main, '/build/meshfree_solver_test.out', 'executable', nil, {'-lm'})
