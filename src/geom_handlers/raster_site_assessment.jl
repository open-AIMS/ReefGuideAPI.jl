"""Geometry-based assessment methods."""

using FLoops

include("geom_ops.jl")
include("common_assessment.jl")

function identify_potential_sites(rst_stack, scan_locs, threshold, geom)
    res = abs(step(dims(rst_stack, X)))
    geom_buff = GO.buffer(geom, res)

    # TODO: Dynamically build this ruleset
    ruleset = Dict(
        :Depth => (data) -> within_thresholds(data, -9.0, -2.0),
        :WavesTp => (data) -> within_thresholds(data, 0.0, 5.9)
    )

    # Search each location to assess
    best_score = zeros(length(scan_locs))
    best_poly = Vector(undef, length(scan_locs))
    best_degree = zeros(Int64, length(scan_locs))
    for (i, (lon_idx, lat_idx)) in enumerate(scan_locs)
        # Move geom to new centroid
        lon = dims(rst_stack, X)[lon_idx]
        lat = dims(rst_stack, Y)[lat_idx]
        mv_geom = move_geom(geom_buff, (lon, lat))

        b_score, b_deg, b_poly = assess_reef_site(rst_stack, mv_geom, ruleset; degree_step=15.0)

        best_score[i] = b_score
        best_degree[i] = b_deg
        best_poly[i] = b_poly
    end

    return DataFrame(score=best_score, orientation=best_degree, poly=best_poly)
end

function assess_reef_site(rst, geom, ruleset; degree_step=15.0)
    n_rotations = length(0.0:degree_step:359.0)
    score = zeros(n_rotations)
    best_poly = Vector(undef, n_rotations)

    for (j, r) in enumerate(0.0:degree_step:359.0)
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

    return score[argmax(score)], argmax(score) - 1, best_poly[argmax(score)]
end

"""
    assess_reef_site(
        rst::Union{Raster,RasterStack},
        geom::GI.Wrappers.Polygon,
        ruleset::Dict{Symbol,Function};
        degree_step::Float64=15.0,
        start_rot::Float64=0.0,
        n_per_side::Int64=1
    )::Tuple{Float64,Int64,GI.Wrappers.Polygon}

Assess given reef site.
"""
function assess_reef_site(
    rst::Union{Raster,RasterStack},
    geom::GI.Wrappers.Polygon,
    ruleset::Dict{Symbol,Function};
    degree_step::Float64=15.0,
    start_rot::Float64=0.0,
    n_per_side::Int64=1
)::Tuple{Float64,Int64,GI.Wrappers.Polygon}
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

    return score[argmax(score)], argmax(score)-(n_per_side+1), best_poly[argmax(score)]
end

"""
    identify_potential_sites_edges(
        rst_stack::RasterStack,
        search_pixels::DataFrame,
        gdf::DataFrame,
        x_dist::Union{Int64,Float64},
        y_dist::Union{Int64,Float64},
        target_crs::GeoFormatTypes.CoordinateReferenceSystemFormat,
        region::String,
        reef_lines::Vector{Vector{GeometryBasics.Line{2,Float64}}};
        degree_step::Float64=15.0,
        n_rot_per_side::Int64=2
    )::DataFrame

Identify the most suitable site polygons for each pixel in the `search_pixels` DataFrame.
`x_dist` and `y_dist` are x and y lengths of the search polygon. A buffer of `rst_stack`
resolution is applied to the search box. And angle from a pixel to a reef edge is identified
and used for searching with custom rotation parameters.

# Arguments
- `rst_stack` : RasterStack containing environmental variables for assessment.
- `search_pixels` : DataFrame containing lon and lat values for each pixel intended for
- `gdf` : GeoDataFrame containing the reef outlines used to align the search box edge.
- `x_dist` : Length of horizontal side of search box.
- `y_dist` : Length of vertical side of search box.
- `target_crs` : CRS of the input Rasters. Using GeoFormatTypes.EPSG().
- `region` : Management region name in GBRMPA format - e.g. "Mackay/Capricorn Management Area"
- `reef_lines` : Vector containing reef outline segments for each reef in `gdf.geometry` (Must be separate object to `gdf` rather than a column).
- `degree_step` : Degree to perform rotations around identified edge angle.
- `n_rot_per_side` : Number of rotations to perform clockwise and anticlockwise around the identified edge angle. Default 2 rotations.

# Returns
DataFrame containing highest score, rotation and polygon for each assessment at pixels in indices.
"""
function identify_potential_sites_edges(
    rst_stack::RasterStack,
    search_pixels::DataFrame,
    gdf::DataFrame,
    x_dist::Union{Int64,Float64},
    y_dist::Union{Int64,Float64},
    target_crs::GeoFormatTypes.CoordinateReferenceSystemFormat,
    region::String,
    reef_lines::Vector{Vector{GeometryBasics.Line{2,Float64}}};
    degree_step::Float64=15.0,
    n_rot_per_side::Int64=2
)::DataFrame
    reef_lines = reef_lines[gdf.management_area .== region]
    gdf = gdf[gdf.management_area .== region, :]
    res = abs(step(dims(rst_stack, X)))

    # # TODO: Dynamically build this ruleset
    ruleset = Dict(
        :Depth => (data) -> within_thresholds(data, -9.0, -2.0),
        :WavesTp => (data) -> within_thresholds(data, 0.0, 5.9)
    )

    # Search each location to assess
    best_score = zeros(length(search_pixels.lon))
    best_poly = Vector(undef, length(search_pixels.lon))
    best_rotation = zeros(Int64, length(search_pixels.lon))
    @floop for (i, index) in enumerate(eachrow(search_pixels))
        lon = index.lon
        lat = index.lat
        geom_buff = initial_search_box((lon, lat), x_dist, y_dist, target_crs, res)

        pixel = GO.Point(lon, lat)
        rot_angle = initial_search_rotation(pixel, geom_buff, gdf, reef_lines)

        b_score, b_rot, b_poly = assess_reef_site(
            rst_stack,
            geom_buff,
            ruleset;
            degree_step=degree_step,
            start_rot=rot_angle,
            n_per_side=n_rot_per_side
        )

        best_score[i] = b_score
        best_rotation[i] = b_rot
        best_poly[i] = b_poly
    end

    return DataFrame(score=best_score, orientation=best_rotation, poly=best_poly)
end