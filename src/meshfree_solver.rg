import "regent"

require "file_properties"
require "parameter"
require "data_structure"
require "point_preprocessor"

local format = require("std/format")

task main()
  var params = init_params()
  -- Define file properties as a region as of now so we can pass it down to tasks
  var file_props = region(ispace(int1d, 1), file_properties)
  file_props[0].max_points = 580485
  var points = region(ispace(int1d, file_props[0].max_points), point)

  format.println("Reading grid from file")
  read_input_point_data(params, points, file_props)
  format.println("Finished reading grid from file")

end

regentlib.saveobj(main, '/build/meshfree_solver_test.out', 'executable', nil, {'-lm'})
