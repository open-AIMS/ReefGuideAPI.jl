using GLMakie, GeoMakie
using LinearAlgebra
import ArchGDAL as AG
import GeoInterface as GI
import GeoDataFrames as GDF

# Testing
test_points = [
    [225, 170],
    [250, 125],
    [225, 70],
    [175, 70],
    [155, 85],
    [120, 160]
]
poly_points = [
    [200, 200],
    [300, 150],
    [250, 75],
    [200, 50],
    [150, 70],
    [100, 150]
]
hex = AG.createpolygon()

ring = AG.createlinearring([poly_points]...)
AG.addgeom!(hex, ring)

box = AG.createpolygon()
ring = AG.createlinearring([(0.0,0.0), (25.0,0.0), (25.0, 7.5), (0.0, 7.5), (0.0, 0.0)])
AG.addgeom!(box, ring)

rot_geom = Vector{Any}(missing, size(test_points, 1))
mv_boxes = Vector{Any}(missing, size(test_points, 1))
for (j, point) in enumerate(test_points)
    pixel = AG.createpoint()
    pixel = AG.addpoint!(pixel, point[1], point[2])
    mv_box = move_geom(box, (point[1], point[2]))
    mv_boxes[j] = mv_box
    search_box_line = find_horiz(mv_box)

    edge_line = identify_closest_edge(pixel, hex)

    edge_bearing = angle_cust([(0.0,5.0), (0.0,0.0)], from_zero(edge_line))
    rot_angle = angle_cust(from_zero(search_box_line), from_zero(edge_line))
    if edge_bearing > 90
        rot_angle = -rot_angle
    end

    rot_geom[j] = rotate_geom(GO.simplify(mv_box; number = 5), rot_angle)
end

poly(hex)
poly!(mv_boxes)
poly!(rot_geom)
