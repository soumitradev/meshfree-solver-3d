import "regent"

fspace file_properties {
  max_points: int,
  wall_points: int,
  outer_points: int,
  interior_points: int,
  supersonic_outlet_points: int,
  supersonic_inlet_points: int,
  shape_points: int,
  number_of_shapes: int,

  res_old: double,
  res_new: double,
  residue: double,
  max_res: double,

  max_res_point: int,
  cfv: double,
  Cl: double,
  Cd: double,
  Cm: double,
  total_entropy: double,
  enstrophy: double
}
