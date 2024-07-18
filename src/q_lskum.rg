import "regent"

require "parameter"
require "data_structure"
require "file_properties"
require "generate_split_stencils"
require "fpi_solver"

local format = require("std/format")
local ctime = terralib.includec("time.h")
local cstdlib = regentlib.c

task q_lskum(
  params: parameters,
  points: region(ispace(int1d), point),
  file_props: region(ispace(int1d, 1), file_properties)
) where reads writes (points, file_props) do
  var file = cstdlib.fopen("/out/residue.dat", "w")
  var start = ctime.clock()

  generate_split_stencils(points, file_props)

  for t = 0, params.max_iters do
    fpi_solver(params, points, file_props, t)
    format.println("{} {} {}", t + 1, file_props[0].res_new, file_props[0].residue)
    cstdlib.fprintf(
      file,
      "%d %lf %lf\n",
      t + 1,
      file_props[0].res_new,
      file_props[0].residue
    )
  end

  var finish = ctime.clock()

  -- WARNING: The constant we divide by is a platform-dependent constant
  -- The ctime variable here unfortunately doesn't contain the value of this macro, and
  -- so we have resorted to hardcoding it. The most reliable way is to print the value
  -- of this macro through a compiled C program via the same installation of Clang that
  -- regent uses for its LLVM runtime on the same Apptainer container, but it should
  -- usually take the same value as it does on the host using GCC.
  format.println("Total time: {}ms", (finish - start) / 1000)
  cstdlib.fclose(file)
end
