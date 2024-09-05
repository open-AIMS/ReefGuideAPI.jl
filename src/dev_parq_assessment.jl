"""Geometry-based assessment methods."""

using Base.Threads
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
using LinearAlgebra

include("geom_handlers/geom_ops.jl")
include("assessment/criteria.jl")
include("geom_handlers/site_assessment.jl")
include("assessment/query_thresholds.jl")

function identify_potential_sites_edges(
    parq_df::DataFrame,
    indices_pixels::Raster,
    indices::Vector{CartesianIndex{2}},
    gdf::DataFrame,
    x_dist::Union{Int64, Float64},
    y_dist::Union{Int64, Float64},
    t_crs::GeoFormatTypes.CoordinateReferenceSystemFormat;
    geometry_col::Symbol=:geometry,
    lines_col::Symbol=:lines,
    degree_step::Float64=15.0,
    n_rot_p_side::Int64=2
)::DataFrame
    res = abs(step(dims(indices_pixels, X)))

    # Search each location to assess
    best_score = zeros(length(indices))
    best_poly = Vector(undef, length(indices))
    best_rotation = zeros(Int64, length(indices))
    @floop for (i, index) in enumerate(indices)
        lon = dims(indices_pixels, X)[index[1]]
        lat = dims(indices_pixels, Y)[index[2]]
        geom_buff, rot_angle = initial_search_box(lon, lat, x_dist, y_dist, res, t_crs, gdf, geometry_col, lines_col)

        bounds = [
            lon - meters_to_degrees(x_dist/2, lat),
            lon + meters_to_degrees(x_dist/2, lat),
            lat - meters_to_degrees(x_dist/2, lat),
            lat + meters_to_degrees(x_dist/2, lat)
        ]

        rel_pix = parq_df[
            (parq_df.lon .> bounds[1]) .&
            (parq_df.lon .< bounds[2]) .&
            (parq_df.lat .> bounds[3]) .&
            (parq_df.lat .< bounds[4])
        , :]

        b_score, b_rot, b_poly = assess_reef_site_df(
            rel_pix,
            geom_buff;
            degree_step=degree_step,
            start_rot=rot_angle,
            n_per_side=n_rot_p_side
        )

        best_score[i] = b_score
        best_rotation[i] = b_rot
        best_poly[i] = b_poly
    end

    return DataFrame(score=best_score, orientation=best_rotation, poly=best_poly)
end

function assess_reef_site_df(rel_pix, geom; degree_step=15.0, start_rot=0.0, n_per_side=2)
    rotations = start_rot-(degree_step*n_per_side):degree_step:start_rot+(degree_step*n_per_side)
    n_rotations = length(rotations)
    score = zeros(n_rotations)
    best_poly = Vector(undef, n_rotations)

    for (j, r) in enumerate(rotations)
        rot_geom = rotate_geom(geom, r)
        score[j] = size(rel_pix[GO.intersects.([geom], rel_pix.geometry), :], 1) / 100
        best_poly[j] = rot_geom
    end

    return score[argmax(score)], argmax(score), best_poly[argmax(score)]
end

gdf = GDF.read("../../canonical-reefs/output/rrap_canonical_2024-07-24-T12-38-38.gpkg")
gdf.buffer = GO.simplify(gdf.geometry; number=30)
for row in eachrow(gdf)
    lat = GO.centroid(row.buffer)[2]
    row.buffer = GO.simplify(GO.buffer(row.buffer, meters_to_degrees(40, lat)); number=30)
end
gdf.lines = polygon_to_lines.(gdf.buffer)

valid_pix = GeoParquet.read("../../GBR-reef-guidance-assessment/outputs/MPA/Mackay-Capricorn_valid_slopes_lookup.parq")
rename!(valid_pix, [:PortDistSlopes => :error_col,:PortDistFlats => :PortDistSlopes])

valid_pix = valid_pix[
    (ruleset[:Depth](valid_pix.Depth) .&
    ruleset[:WavesTp](valid_pix.WavesTp) .&
    ruleset[:Slope](valid_pix.Slope) .&
    ruleset[:WavesHs](valid_pix.WavesHs) .&
    ruleset[:Turbidity](valid_pix.Turbidity) .&
    ruleset[:PortDistSlopes](valid_pix.PortDistSlopes)), :
]
valid_pix.all_crit = (
    ruleset[:Depth](valid_pix.Depth) .&
    ruleset[:WavesTp](valid_pix.WavesTp) .&
    ruleset[:Slope](valid_pix.Slope) .&
    ruleset[:WavesHs](valid_pix.WavesHs) .&
    ruleset[:Turbidity](valid_pix.Turbidity) .&
    ruleset[:PortDistSlopes](valid_pix.PortDistSlopes)
)
valid_pix.lon = first.(GI.coordinates.(valid_pix.geometry))
valid_pix.lat = last.(GI.coordinates.(valid_pix.geometry))

reg = "Mackay-Capricorn"
MPA_OUTPUT_DIR = "c:/Users/bgrier/Documents/Projects/GBR-reef-guidance-assessment/outputs/MPA/"
scan_locs = joinpath(MPA_OUTPUT_DIR, "$(reg)_suitable_slopes.tif")
scan_locs = Raster(scan_locs)
threshold = 95
scan_locs = trim(scan_locs .> threshold)

indices = findall(scan_locs)

x_dist = 450
y_dist = 10

identify_potential_sites_edges(valid_pix, scan_locs, indices, gdf, 450, 10, EPSG(7844))

# RasterStack Comparisons
rst_stack = (
    Depth = joinpath(MPA_OUTPUT_DIR, "$(reg)_bathy.tif"),
    Benthic = joinpath(MPA_OUTPUT_DIR, "$(reg)_benthic.tif"),
    Geomorphic = joinpath(MPA_OUTPUT_DIR, "$(reg)_geomorphic.tif"),
    Slope = joinpath(MPA_OUTPUT_DIR, "$(reg)_slope.tif"),
    Turbidity = joinpath(MPA_OUTPUT_DIR, "$(reg)_turbid.tif"),
    WavesTp = joinpath(MPA_OUTPUT_DIR, "$(reg)_waves_Tp.tif"),
    WavesHs = joinpath(MPA_OUTPUT_DIR, "$(reg)_waves_Hs.tif"),
    #Rugosity = joinpath(MPA_OUTPUT_DIR, "$(reg)_rugosity.tif"),
    PortDistSlopes = joinpath(MPA_OUTPUT_DIR, "$(reg)_port_distance_flats.tif")
    #ValidSlopes = joinpath(MPA_OUTPUT_DIR, "$(reg)_valid_slopes.tif")
    #ValidFlats = joinpath(MPA_OUTPUT_DIR, "$(reg)_valid_flats.tif")
)
rst_stack = RasterStack(rst_stack; lazy=true)

ruleset = Dict(
         :Depth => (data) -> within_thresholds(data, -9.0, -2.0),
         :Slope => (data) -> within_thresholds(data, 0.0, 40.0),
         :WavesTp => (data) -> within_thresholds(data, 0.0, 5.9),
         :WavesHs => (data) -> within_thresholds(data, 0.0, 1.0),
         :Turbidity => (data) -> within_thresholds(data, 0, 58),
         :PortDistSlopes => (data) -> within_thresholds(data, 0.0, 199.9)
)
