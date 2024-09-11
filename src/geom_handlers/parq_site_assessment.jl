"""Geometry-based assessment methods."""

include("geom_ops.jl")
include("common_assessment.jl")

function assess_reef_site(
    rel_pix::DataFrame,
    geom::GI.Wrappers.Polygon,
    max_count::Float64;
    degree_step::Float64=15.0,
    start_rot::Float64=0.0,
    n_per_side::Int64=2,
    surr_threshold::Float64=0.33
)::Tuple{Float64,Int64,GI.Wrappers.Polygon,Int64}
    rotations = (start_rot-(degree_step*n_per_side)):degree_step:(start_rot+(degree_step*n_per_side))
    n_rotations = length(rotations)
    score = zeros(n_rotations)
    best_poly = Vector(undef, n_rotations)
    qc_flag = zeros(Int64, n_rotations)

    for (j, r) in enumerate(rotations)
        rot_geom = rotate_geom(geom, r)
        score[j] = size(rel_pix[GO.intersects.([rot_geom], rel_pix.geometry), :], 1) / max_count
        best_poly[j] = rot_geom

        if score[j] < surr_threshold
            qc_flag[j] = 1
            break
        end
    end

    return score[argmax(score)], argmax(score), best_poly[argmax(score)], maximum(qc_flag)
end

"""
    identify_potential_sites_edges(
        parq_df::DataFrame,
        indices_pixels::Raster,
        indices::Vector{CartesianIndex{2}},
        gdf::DataFrame,
        x_dist::Union{Int64,Float64},
        y_dist::Union{Int64,Float64},
        target_crs::GeoFormatTypes.CoordinateReferenceSystemFormat,
        reef_lines::Vector{Vector{GeometryBasics.Line{2,Float64}}},
        reg::String;
        degree_step::Float64=15.0,
        n_rot_p_side::Int64=2,
        surr_threshold::Float64=0.33
    )::DataFrame

Identify the most suitable site polygons for each pixel in the `search_pixels` raster where
`indices` denotes which pixels to check for suitability. `x_dist` and `y_dist` are x and y
lengths of the search polygon in meters. A buffer of `rst_stack` resolution is applied to the search box.
And angle from a pixel to a reef edge is identified and used for searching with custom rotation
parameters. Method is currently opperating for CRS in degrees units.

# Arguments
- `parq_df` : DataFrame containing environmental variables for assessment.
- `indices_pixels` : Raster that matches indices for lon/lat information.
- `indices` : Vector of CartesianIndices noting pixels to assess sites.
- `gdf` : GeoDataFrame containing the reef outlines used to align the search box edge.
- `x_dist` : Length of horizontal side of search box (in meters).
- `y_dist` : Length of vertical side of search box (in meters).
- `t_crs` : CRS of the input Rasters. Using GeoFormatTypes.EPSG().
- `reg` : Management region name in GBRMPA format - e.g. "Mackay/Capricorn Management Area"
- `geometry_col` : Column name containing target geometries for edge detection in gdf. Should be the same geometry column used to trim and mask the valid search pixels.
- `lines_col` : Column name containing perimeter lines for edge detection in gdf.
- `degree_step` : Degree to perform rotations around identified edge angle.
- `n_rot_p_side` : Number of rotations to perform clockwise and anticlockwise around the identified edge angle. Default 2 rotations.
- `surr_threshold` : Theshold used to skip searching where the proportion of suitable pixels is too low.

# Returns
DataFrame containing highest score, rotation and polygon for each assessment at pixels in indices.
"""
function identify_potential_sites_edges(
    parq_df::DataFrame,
    indices_pixels::Raster,
    indices::Vector{CartesianIndex{2}},
    gdf::DataFrame,
    x_dist::Union{Int64,Float64},
    y_dist::Union{Int64,Float64},
    target_crs::GeoFormatTypes.CoordinateReferenceSystemFormat,
    reef_lines::Vector{Vector{GeometryBasics.Line{2,Float64}}},
    reg::String;
    degree_step::Float64=15.0,
    n_rot_p_side::Int64=2,
    surr_threshold::Float64=0.33
)::DataFrame
    reef_lines = reef_lines[gdf.management_area .== reg]
    gdf = gdf[gdf.management_area .== reg, :]
    res = abs(step(dims(indices_pixels, X)))
    max_count = (
        (x_dist / degrees_to_meters(res, mean(indices_pixels.dims[2]))) *
        ((y_dist + 2 * degrees_to_meters(res, mean(indices_pixels.dims[2]))) /
         degrees_to_meters(res, mean(indices_pixels.dims[2])))
    )

    # Search each location to assess
    best_score = zeros(length(indices))
    best_poly = Vector(undef, length(indices))
    best_rotation = zeros(Int64, length(indices))
    quality_flag = zeros(Int64, length(indices))
    for (i, index) in enumerate(indices)
        lon = dims(indices_pixels, X)[index[1]]
        lat = dims(indices_pixels, Y)[index[2]]
        geom_buff = initial_search_box((lon, lat), x_dist, y_dist, target_crs, res)

        pixel = GO.Point(lon, lat)
        rot_angle = initial_search_rotation(pixel, geom_buff, gdf, reef_lines)

        bounds = [
            lon - meters_to_degrees(x_dist / 2, lat),
            lon + meters_to_degrees(x_dist / 2, lat),
            lat - meters_to_degrees(x_dist / 2, lat),
            lat + meters_to_degrees(x_dist / 2, lat)
        ]

        rel_pix = parq_df[
            (parq_df.lon .> bounds[1]) .& (parq_df.lon .< bounds[2]) .& (parq_df.lat .> bounds[3]).&(parq_df.lat .< bounds[4]), :]

        b_score, b_rot, b_poly, qc_flag = assess_reef_site(
            rel_pix,
            geom_buff,
            max_count;
            degree_step=degree_step,
            start_rot=rot_angle,
            n_per_side=n_rot_p_side,
            surr_threshold=surr_threshold
        )

        best_score[i] = b_score
        best_rotation[i] = b_rot
        best_poly[i] = b_poly
        quality_flag[i] = qc_flag
    end

    return DataFrame(score=best_score, orientation=best_rotation, qc_flag=quality_flag, poly=best_poly)
end
