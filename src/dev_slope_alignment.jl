using GLMakie, GeoMakie
using LinearAlgebra
import ArchGDAL as AG
import GeoInterface as GI
import GeoDataFrames as GDF

include("geom_handlers/site_assessment.jl")

reg = "Townsville-Whitsunday"
MPA_OUTPUT_DIR = "c:/Users/bgrier/Documents/Projects/GBR-reef-guidance-assessment/outputs/MPA/"
# rst_stack = (
#     Depth = joinpath(MPA_OUTPUT_DIR, "$(reg)_bathy.tif"),
#     # Slope = joinpath(MPA_OUTPUT_DIR, "$(reg)_slope.tif"),
#     # Benthic = joinpath(MPA_OUTPUT_DIR, "$(reg)_benthic.tif"),
#     # Geomorphic = joinpath(MPA_OUTPUT_DIR, "$(reg)_geomorphic.tif"),
#     # WavesHs = joinpath(MPA_OUTPUT_DIR, "$(reg)_waves_Hs.tif"),
#     WavesTp = joinpath(MPA_OUTPUT_DIR, "$(reg)_waves_Tp.tif")
#     # Turbidity = joinpath(MPA_OUTPUT_DIR, "$(reg)_turbid.tif"),
#     # Rugosity = joinpath(MPA_OUTPUT_DIR, "$(reg)_rugosity.tif"),
#     # ValidSlopes = joinpath(MPA_OUTPUT_DIR, "$(reg)_valid_slopes.tif"),
#     # ValidFlats = joinpath(MPA_OUTPUT_DIR, "$(reg)_valid_flats.tif")
# )
scan_locs = joinpath(MPA_OUTPUT_DIR, "$(reg)_suitable_slopes_rugosity.tif")
scan_locs = Raster(scan_locs)
scan_locs = Rasters.resample(scan_locs; crs=EPSG(7856))
threshold = 95
scan_locs = trim(scan_locs .> threshold)

# rst_stack = RasterStack(rst_stack; lazy=true)
# rst_stack = Rasters.reproject(rst_stack; crs=EPSG(7856))

# Define the polygon shape to search for (and auto-rotate)
xs = (1, 450)
ys = (1, 10)
search_plot = create_poly(create_bbox(xs, ys), EPSG(7856))
poly(search_plot)
geom=search_plot
# geom = GeometryOps.reproject(geom, EPSG(7856), EPSG(7844))

gdf = GDF.read("../../canonical-reefs/output/rrap_canonical_2024-07-24-T12-38-38.gpkg")
gdf.geometry = AG.reproject(gdf.geometry, EPSG(7844), EPSG(7856); order=:trad)
gdf.geometry = GO.simplify.(gdf.geometry; number=20)
valid_pixels = Rasters.trim(mask(scan_locs; with=gdf.geometry))
indices = findall(valid_pixels)

function find_horiz(geom)
    coords = collect(GI.coordinates(geom)...)
    first_coord = first(coords)
    second_coord = coords[
        (getindex.(coords, 2) .∈ first_coord[2]) .&&
        (getindex.(coords, 1) .∉ first_coord[1])
    ]

    return [tuple(first_coord...), tuple(first(second_coord)...)]
end

function from_zero(v)
    max_coord = maximum(v)
    if first(v) == max_coord
        v = reverse(v)
    end

    new_coords = [Vector{Union{Missing, Float64}}(missing, size(max_coord, 1)), Vector{Union{Missing, Float64}}(missing, size(max_coord, 1))]
    for (j, coords) in enumerate(new_coords)
        for val in eachindex(coords)
            coords[val] = max_coord[val] - v[j][val]
        end
    end

    return new_coords
end

function angle_cust(a, b)
    return acosd(clamp(a⋅b/(norm(a)*norm(b)), -1, 1))
end

function identify_closest_edge(pixel, reef)
   reef_points = GI.coordinates(reef...)
   reef_points = vcat(reef_points...)
   reef_points = Vector{Union{Missing, AG.IGeometry{AG.wkbPoint}}}(missing, size(point_coords, 1))
   for (z, point) in enumerate(point_coords)
       reef_points[z] = AG.createpoint(point)
   end

    distances = GO.distance.([pixel], reef_points)
    # Find the two closest vertices to a pixel
    point_a = tuple(point_coords[distances .== sort(distances)[1]][1]...)
    point_b = tuple(point_coords[distances .== sort(distances)[2]][1]...)
    edge_line = [point_a, point_b]

    return edge_line
end



function identify_potential_sites(rst_stack, scan_locs, indices, threshold, geom)
    res = abs(step(dims(scan_locs, X)))
    #geom_buff = GO.simplify(geom_buff; number = 4) # simplify buffer to 4 vertex polygon

    # # TODO: Dynamically build this ruleset
    # ruleset = Dict(
    #     :Depth => (data) -> within_thresholds(data, -9.0, -2.0),
    #     :WavesTp => (data) -> within_thresholds(data, 0.0, 5.9)
    # )

    # Search each location to assess
    best_score = zeros(length(indices))
    best_poly = Vector(undef, length(indices))
    best_degree = zeros(Int64, length(indices))
    for (i, index) in enumerate(indices)
        # Move geom to new centroid
        lon = dims(valid_pixels, X)[index[1]]
        lat = dims(valid_pixels, Y)[index[2]]

        mv_geom = move_geom(geom, (lon, lat))
        search_box_line = find_horiz(mv_geom)
        geom_buff = GO.buffer(mv_geom, res)

        pixel = AG.createpoint()
        pixel = AG.addpoint!(pixel, lon, lat)
        reef = gdf[GO.within.([pixel], gdf.geometry),:].geometry
        edge_line = identify_closest_edge(pixel, reef)

        # Calculate the angle between the two lines
        edge_bearing = angle_cust([(0.0,5.0), (0.0,0.0)], from_zero(edge_line))
        rot_angle = angle_cust(from_zero(search_box_line), from_zero(edge_line))
        if edge_bearing > 90
            rot_angle = -rot_angle
        end

        b_score, b_deg, b_poly = assess_reef_site(rst_stack, mv_geom, ruleset; degree_step=15.0, start_rot=rot_angle)

        best_score[i] = b_score
        best_degree[i] = b_deg
        best_poly[i] = b_poly
    end

    return DataFrame(score=best_score, orientation=best_degree, poly=best_poly)
end

function assess_reef_site(rst, geom, ruleset; degree_step=15.0, start_rot)
    n_rotations = length(start_rot-degree_step:degree_step:start_rot+degree_step)
    score = zeros(n_rotations)
    best_poly = Vector(undef, n_rotations)

    for (j, r) in enumerate(start_rot-degree_step:degree_step:start_rot+degree_step)
        rot_geom = rotate_geom(geom, r)
        c_rst = crop(rst; to=rot_geom)
        if !all(size(c_rst) .> (0, 0))
            @warn "No data found!"
            continue
        end

        window = trues(size(c_rst))
        for (n, crit_rule) in ruleset
            window .= window .& crit_rule(c_rst[n])
            if count(window) < ceil(length(window) / 3)
                # Stop checking other rules if below hard threshold
                break
            end
        end

        score[j] = mean(window)
        best_poly[j] = rot_geom
    end

    return score[argmax(score)], argmax(score)-1, best_poly[argmax(score)]
end
