import "regent"

local sin = regentlib.sin(double)
local atan = regentlib.atan(double)
local cos = regentlib.cos(double)

struct parameters {
  max_iters: int,

  Mach: double,
  aoa: double,

  gamma: double,
  pi: double,
  theta: double,

  rho_inf: double,
  u1_inf: double,
  u2_inf: double,
  u3_inf: double,
  pr_inf: double,

  power: double,

  second_order_flag: int,
  limiter_flag: int,

  VL_CONST: double,
  CFL: double,
  
  inner_iterations: int,
  initial_conditions_flag: int
}

task init_params(): parameters
  var params: parameters
  params.max_iters =  50000

  params.Mach = 0.7
  params.aoa = 0.0

  params.gamma = 1.4
  params.pi = 4.0*atan(1.0)
  params.theta = params.aoa*params.pi/180.0

  params.rho_inf = 1.0
  params.u1_inf = params.Mach*cos(params.theta)
  params.u2_inf = 0
  params.u3_inf = params.Mach*sin(params.theta)
  params.pr_inf = 1.0/1.4

  params.power = 2.0

  params.second_order_flag = 1
  params.limiter_flag = 1

  params.VL_CONST = 100.0
  params.CFL = 0.4

  params.inner_iterations = 0
  params.initial_conditions_flag = 0

  return params
end
