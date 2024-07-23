import "regent"

require "parameter"
require "data_structure"
require "file_properties"

local cstdlib = regentlib.c

__demand(__leaf)
__demand(__inline)
task venkat_limiter(
  params: parameters,
  points: region(ispace(int1d), point),
  qtilde: double[5],
  i: int
): double[5] where reads (points.{q, qm, min_dist}) do
  var phi: double[5]

  var q: double = 0
  var del_neg: double = 0
  var del_pos: double = 0
  var epsi: double = 0
  var num: double = 0
  var den: double = 0
  var temp: double = 0

  for p = 0, 5 do
    q = points[i].q[p]
    del_neg = qtilde[p] - q

    if (cstdlib.abs(del_neg) <= 10e-6) then
      phi[p] = 1
    elseif (cstdlib.abs(del_neg) > 10e-6) then
      if (del_neg > 0) then
        del_pos = points[i].qm[p] - q
      elseif (del_neg < 0) then
        del_pos = points[i].qm[p + 5] - q
      end

      epsi = cstdlib.pow(params.VL_CONST * points[i].min_dist, 3)
      num = (del_neg*del_pos) + (epsi*epsi)
      num = num*del_neg + 2*del_neg*del_neg*del_pos
      den = (del_pos*del_pos + 2*del_neg*del_neg + del_neg*del_pos + epsi*epsi)*del_neg
      temp = num / den

      if (temp < 1) then
        phi[p] = temp
      else
        phi[p] = 1
      end
    end
  end

  return phi
end
