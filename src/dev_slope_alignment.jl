using GLMakie, GeoMakie
using LinearAlgebra
using Rasters
import ArchGDAL as AG
import GeoInterface as GI
import GeoDataFrames as GDF
import GeoFormatTypes as GFT

using
    Glob,
    TOML

using DataFrames
using OrderedCollections
using Memoization
using SparseArrays

import GeoDataFrames as GDF
using
    GeoParquet,
    Rasters

using
    FLoops,
    HTTP,
    Oxygen

include("assessment/criteria.jl")
include("geom_handlers/site_assessment.jl")
include("assessment/query_thresholds.jl")

reg = "Townsville-Whitsunday"
MPA_OUTPUT_DIR = "c:/Users/bgrier/Documents/Projects/GBR-reef-guidance-assessment/outputs/MPA/"
rst_stack = (
     Depth = joinpath(MPA_OUTPUT_DIR, "$(reg)_bathy.tif"),
    #  Slope = joinpath(MPA_OUTPUT_DIR, "$(reg)_slope.tif"),
    #  Benthic = joinpath(MPA_OUTPUT_DIR, "$(reg)_benthic.tif"),
    #  Geomorphic = joinpath(MPA_OUTPUT_DIR, "$(reg)_geomorphic.tif"),
    #  WavesHs = joinpath(MPA_OUTPUT_DIR, "$(reg)_waves_Hs.tif"),
     WavesTp = joinpath(MPA_OUTPUT_DIR, "$(reg)_waves_Tp.tif")
    # Turbidity = joinpath(MPA_OUTPUT_DIR, "$(reg)_turbid.tif"),
    # Rugosity = joinpath(MPA_OUTPUT_DIR, "$(reg)_rugosity.tif"),
#     # ValidSlopes = joinpath(MPA_OUTPUT_DIR, "$(reg)_valid_slopes.tif"),
#     # ValidFlats = joinpath(MPA_OUTPUT_DIR, "$(reg)_valid_flats.tif")
)
scan_locs = joinpath(MPA_OUTPUT_DIR, "$(reg)_suitable_slopes_rugosity.tif")
scan_locs = Raster(scan_locs)
threshold = 95
scan_locs = trim(scan_locs .> threshold)

rst_stack = RasterStack(rst_stack; lazy=true)

gdf = GDF.read("../canonical-reefs/output/rrap_canonical_2024-07-24-T12-38-38.gpkg")

gdf.geometry = GO.simplify.(gdf.geometry; number=20)
valid_pixels = Rasters.trim(mask(scan_locs; with=gdf.geometry))
indices = findall(valid_pixels)

x_dist = 450
y_dist = 10

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

function meters_to_degrees(x, lat)
    return x / (111.1*1000 * cosd(lat))
end

function identify_closest_edge(pixel, reef)
    if typeof(reef) == Vector{GI.Wrappers.WrapperGeometry{false, false}}
        reef_points = GI.coordinates(reef...)
    else
        reef_points = GI.coordinates(reef)
    end

    reef_points = vcat(reef_points...)
    point_coords = Vector{Union{Missing, AG.IGeometry{AG.wkbPoint}}}(missing, size(reef_points, 1))
    for (z, point) in enumerate(reef_points)
        point_coords[z] = AG.createpoint(point)
    end

    distances = GO.distance.([pixel], point_coords)
    # Find the two closest vertices to a pixel
    point_a = tuple(reef_points[distances .== sort(distances)[1]][1]...)
    point_b = tuple(reef_points[distances .== sort(distances)[2]][1]...)
    edge_line = [point_a, point_b]

    return edge_line
end

function initial_search_box(lon, lat, x_dist, y_dist, res, crs, gdf)
    lon_dist = meters_to_degrees(x_dist, lat)
    xs = (lon - lon_dist/2, lon + lon_dist/2)
    lat_dist = meters_to_degrees(y_dist, lat)
    ys = (lat - lat_dist/2, lat + lat_dist/2)

    search_plot = create_poly(create_bbox(xs, ys), crs)
    search_box_line = find_horiz(search_plot)
    geom_buff = GO.buffer(search_plot, res)

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

    return geom_buff, rot_angle
end

"""
    identify_potential_sites(
        rst_stack::RasterStack,
        search_pixels::Raster,
        indices::Vector{CartesianIndex{2}},
        gdf::DataFrame,
        x_dist::Union{Int64, Float64},
        y_dist::Union{Int64, Float64}
    )::DataFrame

Identify the most suitable site polygons for each pixel in the `search_pixels` raster where
`indices` denotes which pixels to check for suitability. `x_dist` and `y_dist` are x and y
lengths of the search polygon. A buffer of `rst_stack` resolution is applied to the search box.

# Arguments
- `rst_stack` : RasterStack containing environmental variables for assessment.
- `search_pixels` : Raster that matches indices for lon/lat information.
- `indices` : Vector of CartesianIndices noting pixels to assess sites.
- `gdf` : GeoDataFrame containing the reef outlines used to align the search box edge.
- `x_dist` : Length of horizontal side of search box.
- `y_dist` : Length of vertical side of search box.
- `crs` : CRS of the input Rasters. Using GeoFormatTypes.EPSG().

# Returns
DataFrame containing highest score, rotation and polygon for each assessment at pixels in indices.
"""
function identify_potential_sites(
    rst_stack::RasterStack,
    indices_pixels::Raster,
    indices::Vector{CartesianIndex{2}},
    gdf::DataFrame,
    x_dist::Union{Int64, Float64},
    y_dist::Union{Int64, Float64},
    crs::GeoFormatTypes.CoordinateReferenceSystemFormat
)::DataFrame
    res = abs(step(dims(rst_stack, X)))
    #geom_buff = GO.simplify(geom_buff; number = 4) # simplify buffer to 4 vertex polygon

    # # TODO: Dynamically build this ruleset
    ruleset = Dict(
         :Depth => (data) -> within_thresholds(data, -9.0, -2.0),
         :WavesTp => (data) -> within_thresholds(data, 0.0, 5.9)
    )

    # Search each location to assess
    best_score = zeros(length(indices))
    best_poly = Vector(undef, length(indices))
    best_rotation = zeros(Int64, length(indices))
    @floop for (i, index) in enumerate(indices)
        # Move geom to new centroid
        lon = dims(indices_pixels, X)[index[1]]
        lat = dims(indices_pixels, Y)[index[2]]
        geom_buff, rot_angle = initial_search_box(lon, lat, x_dist, y_dist, res, crs, gdf)

        b_score, b_rot, b_poly = assess_reef_site(
            rst_stack,
            geom_buff,
            ruleset;
            degree_step=15.0,
            start_rot=rot_angle,
            n_per_side=2
        )

        best_score[i] = b_score
        best_rotation[i] = b_rot
        best_poly[i] = b_poly
    end

    return DataFrame(score=best_score, orientation=best_rotation, poly=best_poly)
end

function assess_reef_site(rst, geom, ruleset; degree_step=15.0, start_rot=0.0, n_per_side=1)
    rotations = start_rot-(degree_step*n_per_side):degree_step:start_rot+(degree_step*n_per_side)
    n_rotations = length(rotations)
    score = zeros(n_rotations)
    best_poly = Vector(undef, n_rotations)

    for (j, r) in enumerate(rotations)
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

    return score[argmax(score)], argmax(score), best_poly[argmax(score)]
end

TSV_95_slopes = identify_potential_sites(rst_stack, valid_pixels, indices, gdf, 450, 10, EPSG(7844))
reduced_df = TSV_95_slopes
for (ind, row) in enumerate(eachrow(TSV_95_slopes))
    intersecting_polys = TSV_95_slopes[(GO.intersects.([row.poly], TSV_95_slopes[:, :poly])), :]
    reduced_df[ind, :] = TSV_95_slopes[argmax(intersecting_polys.score), :]
end
test = unique(TSV_95_slopes)
