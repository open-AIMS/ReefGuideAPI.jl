"""Geometry-based assessment methods."""

using LinearAlgebra

include("geom_ops.jl")

# Additional functions for reef-edge alignment processing.
function meters_to_degrees(x, lat)
    return x / (111.1*1000 * cosd(lat))
end

function degrees_to_meters(x, lat)
    return x * (111.1*1000 * cosd(lat))
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
    return acosd(clamp(a⋅b / (norm(a) * norm(b)), -1, 1))
end

"""
    filter_far_polygons(gdf, pixel, lat; geometry_col=:geometry)

Filter out reefs that are > 10km from the target pixel (currently hardcoded threshold).
"""
function filter_far_polygons(gdf, pixel, lat; geometry_col=:geometry)
    return gdf[(
        GO.distance.(GO.centroid.(gdf[:, geometry_col]), [pixel]) .< meters_to_degrees(10000, lat)
    ), :]
end

"""
    initial_search_box(
        (lon::Float64, lat::Float64),
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
- Initial search box geometry.
"""
function initial_search_box(
    (lon, lat),
    x_dist::Union{Int64, Float64},
    y_dist::Union{Int64, Float64},
    target_crs::GeoFormatTypes.CoordinateReferenceSystemFormat,
    res::Float64
)::GI.Wrappers.Polygon
    lon_dist = meters_to_degrees(x_dist, lat)
    xs = (lon - lon_dist/2, lon + lon_dist/2)
    lat_dist = meters_to_degrees(y_dist, lat)
    ys = (lat - lat_dist/2, lat + lat_dist/2)

    search_plot = create_poly(create_bbox(xs, ys), target_crs)
    geom_buff = GO.buffer(search_plot, res)

    return geom_buff
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
    pixel::GeometryBasics.Point{2, Float64},
    reef_lines::Vector{GeometryBasics.Line{2, Float64}}
)::Vector{Tuple{Float64, Float64}}
    nearest_edge = reef_lines[argmin(GO.distance.([pixel], reef_lines))]

    return [tuple(x...) for x in nearest_edge]
end

"""
    initial_search_rotation(
        pixel::GeometryBasics.Point{2, Float64},
        geom_buff::GI.Wrappers.Polygon,
        gdf::DataFrame,
        reef_outlines::Vector{Vector{GeometryBasics.Line{2, Float64}}}
    )::Float64

Identifies the closest edge to the target `pixel`/'`geom_buff` and returns the initial rotation
angle required to match the edge line.

# Arguments
- `pixel` : Target point at the center of the search polygon.
- `geom_buff` : Initial search box with zero rotation.
- `gdf` : GeoDataFrame containing a geometry column used for pixel masking.
- `reef_outlines` : Line segments for the outlines of each reef in `gdf`.
"""
function initial_search_rotation(
    pixel::GeometryBasics.Point{2, Float64},
    geom_buff::GI.Wrappers.Polygon,
    gdf::DataFrame,
    reef_outlines::Vector{Vector{GeometryBasics.Line{2, Float64}}}
)::Float64
    reef_lines = reef_outlines[GO.within.([pixel], gdf[:, first(GI.geometrycolumns(gdf))])]
    reef_lines = vcat(reef_lines...)

    # If a pixel is outside of a polygon, use the closest polygon instead.
    if isempty(reef_lines)
        reef_distances = GO.distance.([pixel], gdf[:, first(GI.geometrycolumns(gdf))])
        reef_lines = reef_outlines[argmin(reef_distances)]
        reef_lines = vcat(reef_lines...)
    end

    edge_line = identify_closest_edge(pixel, reef_lines)

    # Calculate the angle between the two lines
    edge_bearing = line_angle([(0.0,5.0), (0.0,0.0)], from_zero(edge_line))
    rot_angle = line_angle(from_zero(find_horizontal(geom_buff)), from_zero(edge_line))
    if edge_bearing > 90
        rot_angle = -rot_angle
    end

    return rot_angle
end

function assess_reef_site(
    rel_pix::DataFrame,
    geom::GI.Wrappers.Polygon,
    max_count::Float64;
    degree_step::Float64=15.0,
    start_rot::Float64=0.0,
    n_per_side::Int64=2,
    surr_threshold::Float64=0.33
)::Tuple{Float64, Int64, GI.Wrappers.Polygon, Int64}
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
        x_dist::Union{Int64, Float64},
        y_dist::Union{Int64, Float64},
        t_crs::GeoFormatTypes.CoordinateReferenceSystemFormat,
        reg::String;
        geometry_col::Symbol=:geometry,
        lines_col::Symbol=:lines,
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
    x_dist::Union{Int64, Float64},
    y_dist::Union{Int64, Float64},
    target_crs::GeoFormatTypes.CoordinateReferenceSystemFormat,
    reef_lines::Vector{Vector{GeometryBasics.Line{2, Float64}}},
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
        ((y_dist + 2*degrees_to_meters(res, mean(indices_pixels.dims[2]))) /
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

"""
    filter_intersecting_sites(res_df::DataFrame)::DataFrame

Filter out sites where the qc_flag indicates a suitabiltiy < `surr_threshold` in searching.
Identify and keep the highest scoring site polygon where site polygons are overlapping.

# Arguments
- `res_df` : Results DataFrame containing potential site polygons (output from
`identify_potential_sites()` or `identify_potential_sites_edges()`).
"""
function filter_intersecting_sites(res_df::DataFrame)::DataFrame
    res_df.row_ID = 1:size(res_df,1)
    ignore_list = []

    for (ind, row) in enumerate(eachrow(res_df))
        if row.row_ID ∈ ignore_list
            continue
        end

        if row.qc_flag == 1
            push!(ignore_list, row.row_ID)
            continue
        end

        if any(GO.intersects.([row.poly], res_df[:,:poly]))
            intersecting_polys = res_df[(GO.intersects.([row.poly], res_df[:,:poly])),:]
            if maximum(intersecting_polys.score) <= row.score
                for x_row in eachrow(intersecting_polys[intersecting_polys.row_ID .!= row.row_ID,:])
                    push!(ignore_list, x_row.row_ID)
                end
            else
                push!(ignore_list, row.row_ID)
            end
        end
    end

    return res_df[res_df.row_ID .∉ [unique(ignore_list)], Not(:row_ID)]
end
