import "regent"

-- Don't worry about Struct-of-Arrays vs Array-of-Struct for now, by default Regent follows SOA.
fspace point {
  x: double,
  y: double,
  z: double,

  tan1: regentlib.array(double, 3),
  tan2: regentlib.array(double, 3),
  nor: regentlib.array(double, 3),

  status: int,
  nbhs: int,
  
  conn: regentlib.array(int, 40),

  min_dist: double,

  xpos_nbhs: int,
  xpos_conn: regentlib.array(int, 25),

  ypos_nbhs: int,
  ypos_conn: regentlib.array(int, 25),
  
  zpos_nbhs: int,
  zpos_conn: regentlib.array(int, 25),

  xneg_nbhs: int,
  xneg_conn: regentlib.array(int, 25),

  yneg_nbhs: int,
  yneg_conn: regentlib.array(int, 25),
  
  zneg_nbhs: int,
  zneg_conn: regentlib.array(int, 25),

  prim: regentlib.array(double, 5),
  flux_res: regentlib.array(double, 5),
  
  q: regentlib.array(double, 5),
  qm: regentlib.array(double, 10),
  dq: regentlib.array(double, 15),

  delt: double,
}
