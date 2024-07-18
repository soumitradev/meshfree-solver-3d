import "regent"

require "data_structure"
require "file_properties"

__demand(__inline)
task get_wall_point_neighbours(points: region(ispace(int1d), point), i: int) where
reads (points.{x, y, z, nbhs, conn, tan1, tan2, nor}),
writes (points.{xpos_conn, xneg_conn, ypos_conn, yneg_conn, zneg_conn}),
reads writes (points.{xpos_nbhs, xneg_nbhs, ypos_nbhs, yneg_nbhs, zneg_nbhs}) do
  var xi = points[i].x
  var yi = points[i].y
  var zi = points[i].z

  points[i].xpos_nbhs = 0
  points[i].xneg_nbhs = 0
  points[i].ypos_nbhs = 0
  points[i].yneg_nbhs = 0
  points[i].zneg_nbhs = 0

  for j = 0, points[i].nbhs do
    var k = points[i].conn[j]
    var xk = points[k].x
    var yk = points[k].y
    var zk = points[k].z

    var delx = xk - xi
    var dely = yk - yi
    var delz = zk - zi

    var dels = delx*points[i].tan1[0] + dely*points[i].tan1[1] +delz*points[i].tan1[2]
    var delt = delx*points[i].tan2[0] + dely*points[i].tan2[1] +delz*points[i].tan2[2]
    var deln = delx*points[i].nor[0] + dely*points[i].nor[1] +delz*points[i].nor[2]

    if (dels <= 0) then
      points[i].xpos_conn[points[i].xpos_nbhs] = k
      points[i].xpos_nbhs = points[i].xpos_nbhs + 1
    else
      points[i].xneg_conn[points[i].xneg_nbhs] = k
      points[i].xneg_nbhs = points[i].xneg_nbhs + 1
    end

    if (delt <= 0) then
      points[i].ypos_conn[points[i].ypos_nbhs] = k
      points[i].ypos_nbhs = points[i].ypos_nbhs + 1
    else
      points[i].yneg_conn[points[i].yneg_nbhs] = k
      points[i].yneg_nbhs = points[i].yneg_nbhs + 1
    end

    points[i].zneg_conn[points[i].zneg_nbhs] = k
    points[i].zneg_nbhs = points[i].zneg_nbhs + 1
  end
end

__demand(__inline)
task get_interior_point_neighbours(points: region(ispace(int1d), point), i: int) where
reads (points.{x, y, z, nbhs, conn, tan1, tan2, nor}),
writes (points.{xpos_conn, xneg_conn, ypos_conn, yneg_conn, zpos_conn, zneg_conn}),
reads writes (points.{xpos_nbhs, xneg_nbhs, ypos_nbhs, yneg_nbhs, zpos_nbhs, zneg_nbhs})
do
  var xi = points[i].x
  var yi = points[i].y
  var zi = points[i].z

  points[i].xpos_nbhs = 0
  points[i].xneg_nbhs = 0
  points[i].ypos_nbhs = 0
  points[i].yneg_nbhs = 0
  points[i].zpos_nbhs = 0
  points[i].zneg_nbhs = 0

  for j = 0, points[i].nbhs do
    var k = points[i].conn[j]
    var xk = points[k].x
    var yk = points[k].y
    var zk = points[k].z

    var delx = xk - xi
    var dely = yk - yi
    var delz = zk - zi

    var dels = delx*points[i].tan1[0] + dely*points[i].tan1[1] +delz*points[i].tan1[2]
    var delt = delx*points[i].tan2[0] + dely*points[i].tan2[1] +delz*points[i].tan2[2]
    var deln = delx*points[i].nor[0] + dely*points[i].nor[1] +delz*points[i].nor[2]

    if (dels <= 0) then
      points[i].xpos_conn[points[i].xpos_nbhs] = k
      points[i].xpos_nbhs = points[i].xpos_nbhs + 1
    else
      points[i].xneg_conn[points[i].xneg_nbhs] = k
      points[i].xneg_nbhs = points[i].xneg_nbhs + 1
    end

    if (delt <= 0) then
      points[i].ypos_conn[points[i].ypos_nbhs] = k
      points[i].ypos_nbhs = points[i].ypos_nbhs + 1
    else
      points[i].yneg_conn[points[i].yneg_nbhs] = k
      points[i].yneg_nbhs = points[i].yneg_nbhs + 1
    end

    if (deln <= 0) then
      points[i].zpos_conn[points[i].zpos_nbhs] = k
      points[i].zpos_nbhs = points[i].zpos_nbhs + 1
    else
      points[i].zneg_conn[points[i].zneg_nbhs] = k
      points[i].zneg_nbhs = points[i].zneg_nbhs + 1
    end
  end
end

__demand(__inline)
task get_outer_point_neighbours(points: region(ispace(int1d), point), i: int) where
reads (points.{x, y, z, nbhs, conn, tan1, tan2, nor}),
writes (points.{xpos_conn, xneg_conn, ypos_conn, yneg_conn, zpos_conn}),
reads writes (points.{xpos_nbhs, xneg_nbhs, ypos_nbhs, yneg_nbhs, zpos_nbhs}) do
  var xi = points[i].x
  var yi = points[i].y
  var zi = points[i].z

  points[i].xpos_nbhs = 0
  points[i].xneg_nbhs = 0
  points[i].ypos_nbhs = 0
  points[i].yneg_nbhs = 0
  points[i].zpos_nbhs = 0

  for j = 0, points[i].nbhs do
    var k = points[i].conn[j]
    var xk = points[k].x
    var yk = points[k].y
    var zk = points[k].z

    var delx = xk - xi
    var dely = yk - yi
    var delz = zk - zi

    var dels = delx*points[i].tan1[0] + dely*points[i].tan1[1] +delz*points[i].tan1[2]
    var delt = delx*points[i].tan2[0] + dely*points[i].tan2[1] +delz*points[i].tan2[2]
    var deln = delx*points[i].nor[0] + dely*points[i].nor[1] +delz*points[i].nor[2]

    if (dels <= 0) then
      points[i].xpos_conn[points[i].xpos_nbhs] = k
      points[i].xpos_nbhs = points[i].xpos_nbhs + 1
    else
      points[i].xneg_conn[points[i].xneg_nbhs] = k
      points[i].xneg_nbhs = points[i].xneg_nbhs + 1
    end

    if (delt <= 0) then
      points[i].ypos_conn[points[i].ypos_nbhs] = k
      points[i].ypos_nbhs = points[i].ypos_nbhs + 1
    else
      points[i].yneg_conn[points[i].yneg_nbhs] = k
      points[i].yneg_nbhs = points[i].yneg_nbhs + 1
    end

    points[i].zpos_conn[points[i].zpos_nbhs] = k
    points[i].zpos_nbhs = points[i].zpos_nbhs + 1
  end
end

task generate_split_stencils(
  points: region(ispace(int1d), point),
  file_props: region(ispace(int1d, 1), file_properties)
) where
reads (points.{x, y, z, nbhs, conn, status, tan1, tan2, nor}, file_props.{max_points}),
writes (points.{xpos_conn, xneg_conn, ypos_conn, yneg_conn, zpos_conn, zneg_conn}),
reads writes (points.{xpos_nbhs, xneg_nbhs, ypos_nbhs, yneg_nbhs, zpos_nbhs, zneg_nbhs})
do
  for i = 0, file_props[0].max_points do
    if (points[i].status == 0) then
      get_wall_point_neighbours(points, i)
    elseif (points[i].status == 1) then
      get_interior_point_neighbours(points, i)
    elseif (points[i].status == 2) then
      get_outer_point_neighbours(points, i)
    end
  end
end
