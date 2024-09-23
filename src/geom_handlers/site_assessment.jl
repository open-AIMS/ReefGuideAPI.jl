"""Geometry-based assessment methods."""

using LinearAlgebra
using FLoops

include("geom_ops.jl")


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

# Additional functions for reef-edge alignment processing.

"""
    meters_to_degrees(x, lat)

Convert meters to degrees.
"""
function meters_to_degrees(x, lat)
    return x / (111.1 * 1000 * cosd(lat))
end

"""
    from_zero(v::Vector{Tuple{Float64,Float64}})::Vector{Tuple{Float64, Float64}}

Translates Vector of points `v` to begin from (0, 0), retaining direction and length.
"""
function from_zero(v::Vector{Tuple{Float64,Float64}})::Vector{Tuple{Float64,Float64}}
    max_coord = maximum(v)
    if first(v) == max_coord
        v = reverse(v)
    end

    new_coords = [
        Vector{Union{Missing,Float64}}(missing, size(max_coord, 1)),
        Vector{Union{Missing,Float64}}(missing, size(max_coord, 1))
    ]
    for (j, coords) in enumerate(new_coords)
        for val in eachindex(coords)
            coords[val] = max_coord[val] - v[j][val]
        end
    end

    return Tuple.(new_coords)
end

"""
    line_angle(a::T, b::T)::Float64 where {T <: Vector{Tuple{Float64,Float64}}}

Calculate the angle between two lines.

# Arguments
- `a` : Line between point coordinates.
- `b` : Line between point coordinates.

# Returns
Angle between the two lines.

# Examples
```julia
line_angle([(0.0,5.0), (0.0,0.0)], from_zero(edge_line))
line_angle([(0.0,5.0), (0.0,0.0)], [(1.0, 4.0), (7.0, 8.0)])
```
"""
function line_angle(a::T, b::T)::Float64 where {T<:Vector{Tuple{Float64,Float64}}}
    return acosd(clamp(a ⋅ b / (norm(a) * norm(b)), -1, 1))
end

"""
    identify_closest_edge(
        pixel::GeometryBasics.Point{2, Float64},
        reef_lines::Vector{GeometryBasics.Line{2, Float64}}
    )::Vector{Tuple{Float64, Float64}}

Find the nearest line in `reef_lines` to a point `pixel`.

# Arguments
- `pixel` : Target point geometry.
- `reef_lines` : Vector containing lines for comparison.
"""
function identify_closest_edge(
    pixel::GeometryBasics.Point{2,Float64},
    reef_lines::Vector{GeometryBasics.Line{2,Float64}}
)::Vector{Tuple{Float64,Float64}}
    nearest_edge = reef_lines[argmin(GO.distance.([pixel], reef_lines))]

    return [tuple(x...) for x in nearest_edge]
end

"""
    initial_search_box(
        (lon, lat)::Tuple{Float64,Float64},
        x_dist::Union{Int64, Float64},
        y_dist::Union{Int64, Float64},
        target_crs::GeoFormatTypes.CoordinateReferenceSystemFormat,
        res::Float64
    )::GI.Wrappers.Polygon

Create an initial search box that is centered around the point `(lon, lat)` in `target_crs`,
and is buffered by `res` distance.

# Arguments
- `(lon, lat)` : Longitude and latitude coordinates of the center target pixel.
- `x_dist` : x (longitude) dimension length of initial search box.
- `y_dist` : y (latitude) dimension length of initial search box.
- `target_crs` : Target CRS of box to match input data types.
- `res` : Buffer distance (resolution of input raster search data).

# Returns
Initial search box geometry.
"""
function initial_search_box(
    (lon, lat)::Tuple{Float64,Float64},
    x_dist::Union{Int64,Float64},
    y_dist::Union{Int64,Float64},
    target_crs::GeoFormatTypes.CoordinateReferenceSystemFormat,
    res::Float64
)::GI.Wrappers.Polygon
    lon_dist = meters_to_degrees(x_dist, lat)
    xs = (lon - lon_dist / 2, lon + lon_dist / 2)
    lat_dist = meters_to_degrees(y_dist, lat)
    ys = (lat - lat_dist / 2, lat + lat_dist / 2)

    search_plot = create_poly(create_bbox(xs, ys), target_crs)
    geom_buff = GO.buffer(search_plot, res)

    return geom_buff
end

"""
    initial_search_rotation(
        pixel::GeometryBasics.Point{2, Float64},
        geom_buff::GI.Wrappers.Polygon,
        gdf::DataFrame,
        reef_lines::Vector{Vector{GeometryBasics.Line{2, Float64}}}
    )::Float64

Identifies the closest edge to the target `pixel`/'`geom_buff` and returns the initial rotation
angle required to match the edge line.

# Arguments
- `pixel` : Target point at the center of the search polygon.
- `geom_buff` : Initial search box with zero rotation.
- `gdf` : GeoDataFrame containing a geometry column used for pixel masking.
- `reef_lines` : Line segments for the outlines of each reef in `gdf`.

# Returns
Rotation angle of initial search polygon
"""
function initial_search_rotation(
    pixel::GeometryBasics.Point{2,Float64},
    geom_buff::GI.Wrappers.Polygon,
    gdf::DataFrame,
    reef_lines::Vector{Vector{GeometryBasics.Line{2,Float64}}}
)::Float64
    reef_lines = reef_lines[GO.within.([pixel], gdf[:, first(GI.geometrycolumns(gdf))])]
    reef_lines = vcat(reef_lines...)
    edge_line = identify_closest_edge(pixel, reef_lines)

    # Calculate the angle between the two lines
    edge_bearing = line_angle([(0.0, 5.0), (0.0, 0.0)], from_zero(edge_line))
    rot_angle = line_angle(from_zero(find_horizontal(geom_buff)), from_zero(edge_line))
    if edge_bearing > 90
        rot_angle = -rot_angle
    end

    return rot_angle
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

    return score[argmax(score)], argmax(score), best_poly[argmax(score)]
end

"""
    identify_potential_sites_edges(
        rst_stack::RasterStack,
        indices_pixels::Raster,
        indices::Vector{CartesianIndex{2}},
        gdf::DataFrame,
        x_dist::Union{Int64, Float64},
        y_dist::Union{Int64, Float64},
        target_crs::GeoFormatTypes.CoordinateReferenceSystemFormat,
        reef_lines::Vector{Vector{GeometryBasics.Line{2, Float64}}};
        degree_step::Float64=15.0,
        n_rot_per_side::Int64=2
    )::DataFrame

Identify the most suitable site polygons for each pixel in the `indices_pixels` raster where
`indices` denotes which pixels to check for suitability. `x_dist` and `y_dist` are x and y
lengths of the search polygon. A buffer of `rst_stack` resolution is applied to the search box.
And angle from a pixel to a reef edge is identified and used for searching with custom rotation
parameters.

# Arguments
- `rst_stack` : RasterStack containing environmental variables for assessment.
- `indices_pixels` : Raster that matches indices for lon/lat information.
- `indices` : Vector of CartesianIndices noting pixels to assess sites.
- `gdf` : GeoDataFrame containing the reef outlines used to align the search box edge.
- `x_dist` : Length of horizontal side of search box.
- `y_dist` : Length of vertical side of search box.
- `target_crs` : CRS of the input Rasters. Using GeoFormatTypes.EPSG().
- `reef_lines` : Vector containing reef outline segments for each reef in `gdf`.
- `degree_step` : Degree to perform rotations around identified edge angle.
- `n_rot_per_side` : Number of rotations to perform clockwise and anticlockwise around the identified edge angle. Default 2 rotations.

# Returns
DataFrame containing highest score, rotation and polygon for each assessment at pixels in indices.
"""
function identify_potential_sites_edges(
    rst_stack::RasterStack,
    indices_pixels::Raster,
    indices::Vector{CartesianIndex{2}},
    gdf::DataFrame,
    x_dist::Union{Int64,Float64},
    y_dist::Union{Int64,Float64},
    target_crs::GeoFormatTypes.CoordinateReferenceSystemFormat,
    reef_lines::Vector{Vector{GeometryBasics.Line{2,Float64}}};
    degree_step::Float64=15.0,
    n_rot_per_side::Int64=2
)::DataFrame
    res = abs(step(dims(rst_stack, X)))

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
        lon = dims(indices_pixels, X)[index[1]]
        lat = dims(indices_pixels, Y)[index[2]]
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

"""
    filter_intersecting_sites(res_df::DataFrame)::DataFrame

Identify and keep the highest scoring site polygon where site polygons are overlapping.

# Arguments
- `res_df` : Results DataFrame containing potential site polygons (output from
`identify_potential_sites()` or `identify_potential_sites_edges()`).
"""
function filter_intersecting_sites(res_df::DataFrame)::DataFrame
    res_df.row_ID = 1:size(res_df, 1)
    ignore_list = []

    for row in eachrow(res_df)
        if row.row_ID ∈ ignore_list
            continue
        end

        if any(GO.intersects.([row.poly], res_df[:, :poly]))
            intersecting_polys = res_df[(GO.intersects.([row.poly], res_df[:, :poly])), :]
            if maximum(intersecting_polys.score) <= row.score
                for x_row in eachrow(intersecting_polys[intersecting_polys.row_ID .!= row.row_ID, :])
                    push!(ignore_list, x_row.row_ID)
                end
            else
                push!(ignore_list, row.row_ID)
            end
        end
    end

    return res_df[res_df.row_ID .∉ [unique(ignore_list)], Not(:row_ID)]
end