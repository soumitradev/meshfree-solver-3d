import "regent"

require "data_structure"

local cstdlib = regentlib.c

task print_output(points: region(ispace(int1d), point)) where
reads (points.{x, y, z, prim})
do
  var solution_file = cstdlib.fopen("/out/solution.dat", "w")
  var restart_file = cstdlib.fopen("/out/restart.dat", "w")

  for i in points do
    cstdlib.fprintf(
      solution_file,
      "%lf %lf %lf %lf %lf %lf %lf %lf\n",
      points[i].x,
      points[i].y,
      points[i].z,
      points[i].prim[0],
      points[i].prim[1],
      points[i].prim[2],
      points[i].prim[3],
      points[i].prim[4]
    )
    cstdlib.fprintf(
      restart_file,
      "%lf %lf %lf %lf %lf\n",
      points[i].prim[0],
      points[i].prim[1],
      points[i].prim[2],
      points[i].prim[3],
      points[i].prim[4]
    )
  end

  cstdlib.fclose(solution_file)
  cstdlib.fclose(restart_file)
end
