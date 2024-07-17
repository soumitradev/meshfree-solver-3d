import "regent"

require "data_structure"
require "file_properties"

local cstdlib = regentlib.c

task read_input_point_data(
  params: parameters,
  points: region(ispace(int1d), point),
  file_props: region(ispace(int1d, 1), file_properties)
) where reads writes (points, file_props) do
  var file = cstdlib.fopen("/data/3d_input_data", "r")

  var p: point
  var temp: int
  var counter: int
  
  cstdlib.fscanf(file, "%d", &temp)

  for i = 0, file_props[0].max_points do
    cstdlib.fscanf(
      file,
      "%d %lf %lf %lf %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d",
      &counter, &p.x, &p.y, &p.z, &p.status, &p.min_dist,
      &p.tan1[0], &p.tan1[1], &p.tan1[2], &p.tan2[0], &p.tan2[1], &p.tan2[2],
      &p.nor[0], &p.nor[1], &p.nor[2], &p.nbhs
    )
    for j = 0, p.nbhs do
      cstdlib.fscanf(file, "%d", &p.conn[j])
    end
    points[i] = p
  end

  file_props[0].wall_points = 0
  file_props[0].outer_points = 0
  file_props[0].interior_points = 0
  file_props[0].supersonic_outlet_points = 0
  file_props[0].supersonic_inlet_points = 0

  for i = 0, file_props[0].max_points do
    if (points[i].status == 0) then
      file_props[0].interior_points = file_props[0].interior_points + 1
    elseif (points[i].status == 1) then
      file_props[0].wall_points = file_props[0].wall_points + 1
    elseif (points[i].status == 2) then
      file_props[0].outer_points = file_props[0].outer_points + 1
    elseif (points[i].status == 5) then
      file_props[0].supersonic_inlet_points = file_props[0].supersonic_inlet_points + 1
    elseif (points[i].status == 6) then
      file_props[0].supersonic_outlet_points = file_props[0].supersonic_outlet_points + 1
    end
  end
end
