""" Functions common to both site_assessment methods."""

using LinearAlgebra
using GeoDataFrames

"""
    meters_to_degrees(x, lat)

Convert meters to degrees at target latitude.
"""
function meters_to_degrees(x, lat)
    return x / (111.1 * 1000 * cosd(lat))
end

"""
    degrees_to_meters(x, lat)

Convert degrees to meters at target latitude.
"""
function degrees_to_meters(x, lat)
    return x * (111.1 * 1000 * cosd(lat))
end

"""
    from_zero(v::Vector{Tuple{Float64,Float64}})::Vector{Tuple{Float64, Float64}}

Translates Vector of points `v` to begin from (0, 0), retaining direction and length.

# Argument
- `v` : Vector of point coordinates (`Tuple{Float64, Float64}`).
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
    filter_far_reefs(gdf::DataFrame, pixel::GIWrap.Point, lat::Float64, dist::Union{Int64,Float64})::BitVector

Filter out reefs that are > `dist` (meters) from the target pixel (currently `dist` is hardcoded in `initial_search_rotation()`).
"""
function filter_far_reefs(
    geoms,
    pixel::GeometryBasics.Point,
    lat::Float64,
    dist::Union{Int64,Float64}
)::BitVector
    return (GO.distance.(GO.centroid.(geoms), [pixel]) .< meters_to_degrees(dist, lat))
end

"""
    initial_search_box(
        (lon::Float64, lat::Float64),
        x_dist::Union{Int64, Float64},
        y_dist::Union{Int64, Float64},
        target_crs::GeoFormatTypes.EPSG
    )::GI.Wrappers.Polygon

Create an initial search box that is centered around the point `(lon, lat)` in `target_crs`.

# Arguments
- `(lon, lat)` : Longitude and latitude coordinates of the center target pixel.
- `x_dist` : x (longitude) dimension length of initial search box.
- `y_dist` : y (latitude) dimension length of initial search box.
- `target_crs` : Target CRS of box to match input data types.

# Returns
Initial search box geometry.
"""
function initial_search_box(
    (lon, lat),
    x_dist::Union{Int64,Float64},
    y_dist::Union{Int64,Float64},
    target_crs::GeoFormatTypes.EPSG
)::GI.Wrappers.Polygon
    lon_dist = meters_to_degrees(x_dist, lat)
    xs = (lon - lon_dist / 2, lon + lon_dist / 2)
    lat_dist = meters_to_degrees(y_dist, lat)
    ys = (lat - lat_dist / 2, lat + lat_dist / 2)

    search_plot = create_poly(create_bbox(xs, ys), target_crs)

    return search_plot
end

"""
    closest_reef_edge(
        pixel::GeometryBasics.Point{2, Float64},
        reef_lines::Vector{T}
    )::Vector{Tuple{Float64, Float64}} where {T}

Find the nearest line in `reef_lines` to a point `pixel`.

# Arguments
- `pixel` : Target point geometry.
- `reef_lines` : Vector containing lines for comparison.

# Returns
Coordinates of the reef edge line that is closest to the target `pixel`. Returned in Tuples.
"""
function closest_reef_edge(
    pixel::GeometryBasics.Point{2,Float64},
    reef_lines::Vector{T}
)::Vector{Tuple{Float64,Float64}} where {T}
    nearest_edge = reef_lines[argmin(GO.distance.([pixel], reef_lines))]
    return Tuple.(coordinates(nearest_edge))
end

"""
    initial_search_rotation(
        pixel::GeometryBasics.Point{2, Float64},
        geom_buff::GI.Wrappers.Polygon,
        gdf::DataFrame,
        reef_outlines::Vector{Vector{T}};
        search_buffer::Union{Int64,Float64}=20000.0
    )::Float64 where {T}

Identifies the closest edge to the target `pixel`/`geom_buff`. The angle required to rotate
`geom_buff` by to match this reef edge is calculated.

# Extended help
The returned angle is the angle relative to the default `geom_buff` horizontal orientation.
Therefore, if returned angle = 45 degrees, `rot_geom(geom_buff, 45)` will rotate `geom_buff`
by 45 degrees to match the identified reef edge.

# Arguments
- `pixel` : Target point at the center of the search polygon.
- `geom_buff` : Initial search box with zero rotation.
- `reef_geoms` : Geometries representing reefs
- `reef_outlines` : Line segments for the outlines of each reef in `gdf`.

# Returns
Rotation angle required to match reef edge when used in `rotate_geom(geom_buff, rot_angle)`.
"""
function initial_search_rotation(
    pixel::GeometryBasics.Point{2,Float64},
    geom_buff::GI.Wrappers.Polygon,
    reef_geoms::Vector{ArchGDAL.IGeometry},
    reef_outlines::Vector{Vector{T}}
)::Float64 where {T}
    # Get closest reef outline
    closest_reef_outline_idx = try
        in_reef = findall(GO.within.([pixel], reef_geoms))
        first(in_reef)
    catch err
        if !(err isa BoundsError)
            rethrow(err)
        end

        # Pixel is outside a reef outline so get the closest one instead
        # Note: distance will be 0 if point is inside the geometry
        argmin(GO.distance.([pixel], reef_geoms))
    end

    reef_lines = reef_outlines[closest_reef_outline_idx]
    edge_line = closest_reef_edge(pixel, reef_lines)

    # Calculate the angle between the two lines
    edge_bearing = line_angle([(0.0, 5.0), (0.0, 0.0)], from_zero(edge_line))
    rot_angle = line_angle(from_zero(find_horizontal(geom_buff)), from_zero(edge_line))
    if edge_bearing > 90
        rot_angle = -rot_angle
    end

    return rot_angle
end

"""
    filter_sites(res_df::DataFrame)::DataFrame

Filter out sites where the qc_flag indicates a suitability < `surr_threshold` in searching.
Where site polygons are overlapping, keep only the highest scoring site polygon.

# Arguments
- `res_df` : Results DataFrame containing potential site polygons
             (output from `identify_potential_sites()` or `identify_edge_aligned_sites()`).

# Returns
DataFrame containing only the highest scoring sites where site polygons intersect, and
containing only sites with scores greater than the `surr_threshold` specified in
`identify_edge_aligned_sites()` (default=0.6).
"""
function filter_sites(res_df::DataFrame)::DataFrame
    if nrow(res_df) == 0
        return res_df
    end

    res_df.row_ID = 1:size(res_df, 1)
    ignore_list = Int64[]

    for row in eachrow(res_df)
        if row.row_ID ∈ ignore_list
            continue
        end

        not_ignored = res_df.row_ID .∉ Ref(ignore_list)
        poly = row.geometry
        poly_interx = GO.intersects.(Ref(poly), res_df[not_ignored, :geometry])

        if count(poly_interx) > 1
            intersecting_polys = res_df[not_ignored, :][poly_interx, :]

            # Find the ID of the polygon with the best score.
            # Add polygon IDs with non-maximal score to the ignore list
            best_poly_idx = argmax(intersecting_polys.score)
            best_score = intersecting_polys[best_poly_idx, :score]
            if best_score < 0.7
                # Ignore everything as they are below useful threshold
                append!(ignore_list, intersecting_polys[:, :row_ID])
                continue
            end

            best_ID = intersecting_polys[best_poly_idx, :row_ID]
            not_best_ID = intersecting_polys.row_ID .!= best_ID
            append!(ignore_list, intersecting_polys[not_best_ID, :row_ID])
        elseif count(poly_interx) == 1 && (row.score < 0.7)
            # Remove current poly if below useful threshold
            append!(ignore_list, Ref(row.row_ID))
        end
    end

    return res_df[res_df.row_ID .∉ [unique(ignore_list)], :]
end

"""
    output_geojson(destination_path::String, df::DataFrame)::Nothing

Writes out GeoJSON file to a target directory. Output file will be located at location:
`destination_path`.

# Arguments
- `destination_path` : File path to write geojson file to.
- `df` : DataFrame intended for writing to geojson file.
"""
function output_geojson(
    destination_path::String, df::DataFrame; output_crs=EPSG(4326)
)::Nothing
    out_df = copy(df)
    out_df.geometry = GO.reproject(
        out_df.geometry, GI.crs(first(out_df.geometry)), output_crs
    )

    GDF.write(destination_path, out_df; crs=GI.crs(first(out_df.geometry)))

    return nothing
end


"""
    buffer_simplify(
        gdf::DataFrame;
        number_verts::Int64=30,
        buffer_dist_m::Int64=40
    )::Vector{GeoInterface.Wrappers.WrapperGeometry}

Simplify and buffer the polygons in a GeoDataFrame to account for uncertainty and inaccuracies
in the reef outlines.

# Arguments
- `gdf` : GeoDataFrame containing the reef polygons in `gdf.geometry`.
- `number_verts` : Number of vertices to simplify the reefs to. Default is 30 vertices.
- `buffer_dist_m` : Buffering distance in meters to account for innacuracies in reef outlines. Default distance is 20m.

# Returns
Vector containing buffered and simplified reef polygons
"""
function buffer_simplify(
    gdf::DataFrame;
    number_verts::Int64=30,
    buffer_dist_m::Int64=20
)::Vector{GIWrap.WrapperGeometry}
    reef_buffer = GO.simplify(gdf.geometry; number=number_verts)
    for row in eachrow(reef_buffer)
        lat = GO.centroid(row)[2]
        row = GO.simplify(
            GO.buffer(row, meters_to_degrees(buffer_dist_m, lat));
            number=number_verts
        )
    end

    return reef_buffer
end
