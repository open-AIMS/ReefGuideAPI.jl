"""Geometry-based assessment methods."""

# Tabular data assessment methods

"""
    assess_reef_site(
        rel_pix::DataFrame,
        geom::GI.Wrappers.Polygon,
        max_count::Float64,
        target_crs::GeoFormatTypes.CoordinateReferenceSystemFormat;
        degree_step::Float64=15.0,
        start_rot::Float64=0.0,
        n_per_side::Int64=2,
        surr_threshold::Float64=0.33
    )::Tuple{Float64,Int64,GI.Wrappers.Polygon,Int64}

Assesses the rotations of a search box `geom` for their suitability score (calculated as the
proportion of pixels that meet all specified criteria thresholds). Search box rotation steps
are returned so that the `start_rot` angle is 0, rotations anti-clockwise are negative and
rotations clockwise are
positive.

# Arguments
- `rel_pix` : DataFrame containing the point data for pixels that are within maxmimum user search box dimensions from a pixel.
- `geom` : Starting search box for assessment.
- `max_count` : The maximum number of pixels that can intersect the search box (used to standardise scores between 0 and 1).
- `target_crs` : Coordinate Reference System used for analysis vector and raster data.
- `degree_step` : Step to vary the search box rotations.
- `start_rot` : Starting angle rotation that aligns the box with the closest reef edge.
- `n_per_side` : Number of rotations to perform around the starting search box angle.
- `surr_threshold` : Suitability threshold, below which sites are excluded from result sets.

# Returns
- Highest score
- Highest scoring rotation step
- Highest scoring polygon
- Quality control flag for site, indicating if `surr_threshold` was met in the highest scoring rotation.
"""
function assess_reef_site(
    rel_pix::DataFrame,
    geom::GI.Wrappers.Polygon,
    max_count::Float64,
    target_crs::GeoFormatTypes.CoordinateReferenceSystemFormat;
    degree_step::Float64=15.0,
    start_rot::Float64=0.0,
    n_per_side::Int64=2,
    surr_threshold::Float64=0.33
)::Tuple{Float64,Int64,GI.Wrappers.Polygon,Int64}
    rotations =
        (start_rot - (degree_step * n_per_side)):degree_step:(start_rot + (degree_step * n_per_side))
    n_rotations = length(rotations)
    score = zeros(n_rotations)
    best_poly = Vector(undef, n_rotations)
    qc_flag = zeros(Int64, n_rotations)

    for (j, r) in enumerate(rotations)
        rot_geom = rotate_geom(geom, r, target_crs)
        score[j] =
            size(rel_pix[GO.intersects.([rot_geom], rel_pix.geometry), :], 1) / max_count
        best_poly[j] = rot_geom

        if score[j] < surr_threshold
            qc_flag[j] = 1
            break
        end
    end

    return score[argmax(score)],
    argmax(score) - (n_per_side + 1), best_poly[argmax(score)],
    maximum(qc_flag)
end

"""
    identify_edge_aligned_sites(
        df::DataFrame,
        search_pixels::DataFrame,
        res::Float64,
        gdf::DataFrame,
        x_dist::Union{Int64,Float64},
        y_dist::Union{Int64,Float64},
        target_crs::GeoFormatTypes.CoordinateReferenceSystemFormat,
        reef_lines::Vector{Vector{GeometryBasics.Line{2,Float64}}},
        region::String;
        degree_step::Float64=15.0,
        n_rot_p_side::Int64=2,
        surr_threshold::Float64=0.33
    )::DataFrame

Identify the most suitable site polygons for each pixel in the `search_pixels` DataFrame.
`x_dist` and `y_dist` are x and y lengths of the search polygon in meters. A buffer of the
raster files' resolution is applied to the search box. And angle from a pixel to a reef edge
is identified and used for searching with custom rotation parameters.
Method is currently opperating for CRS in degrees units.

# Arguments
- `df` : DataFrame containing environmental variables for assessment.
- `search_pixels` : DataFrame containing lon and lat columns for each pixel that is intended for analysis.
- `res` : Resolution of the original raster pixels. Can by found via `abs(step(dims(raster, X)))`.
- `gdf` : GeoDataFrame containing the reef outlines used to align the search box edge.
- `x_dist` : Length of horizontal side of search box (in meters).
- `y_dist` : Length of vertical side of search box (in meters).
- `target_crs` : CRS of the input Rasters. Using GeoFormatTypes.EPSG().
- `reef_lines` : Vector containing reef outline segments created from polygons in `gdf.geometry` (Must be separate object to `gdf` rather than column).
- `region` : Region name, e.g. "Cairns-Cooktown" or "FarNorthern".
- `degree_step` : Degree to perform rotations around identified edge angle.
- `n_rot_p_side` : Number of rotations to perform clockwise and anticlockwise around the identified edge angle. Default 2 rotations.
- `surr_threshold` : Theshold used to skip searching where the proportion of suitable pixels is too low.

# Returns
DataFrame containing highest score, rotation and polygon for each assessment at pixels in indices.
"""
function identify_edge_aligned_sites(
    df::DataFrame,
    search_pixels::DataFrame,
    res::Float64,
    gdf::DataFrame,
    x_dist::Union{Int64,Float64},
    y_dist::Union{Int64,Float64},
    target_crs::GeoFormatTypes.CoordinateReferenceSystemFormat,
    reef_lines::Vector{Vector{GeometryBasics.Line{2,Float64}}},
    region::String;
    degree_step::Float64=15.0,
    n_rot_p_side::Int64=2,
    surr_threshold::Float64=0.33
)::DataFrame
    region_long = REGIONAL_DATA["region_long_names"][region]
    reef_lines = reef_lines[gdf.management_area .== region_long]
    gdf = gdf[gdf.management_area .== region_long, :]
    max_count = (
        (x_dist / degrees_to_meters(res, mean(search_pixels.lon))) *
        (
            (y_dist + 2 * degrees_to_meters(res, mean(search_pixels.lat))) /
            degrees_to_meters(res, mean(search_pixels.lat))
        )
    )

    # Search each location to assess
    best_score = zeros(length(search_pixels.lon))
    best_poly = Vector(undef, length(search_pixels.lon))
    best_rotation = zeros(Int64, length(search_pixels.lon))
    quality_flag = zeros(Int64, length(search_pixels.lon))
    for (i, index) in enumerate(eachrow(search_pixels))
        lon = index.lon
        lat = index.lat
        geom_buff = initial_search_box((lon, lat), x_dist, y_dist, target_crs, res)

        pixel = GO.Point(lon, lat)
        rot_angle = initial_search_rotation(pixel, geom_buff, gdf, reef_lines)

        bounds = [
            lon - meters_to_degrees(x_dist / 2, lon),
            lon + meters_to_degrees(x_dist / 2, lon),
            lat - meters_to_degrees(y_dist / 2, lat),
            lat + meters_to_degrees(y_dist / 2, lat)
        ]

        rel_pix = df[
            (df.lon .> bounds[1]) .& (df.lon .< bounds[2]) .& (df.lat .> bounds[3]) .& (df.lat .< bounds[4]),
            :]

        b_score, b_rot, b_poly, qc_flag = assess_reef_site(
            rel_pix,
            geom_buff,
            max_count,
            target_crs;
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

    return DataFrame(;
        score=best_score, orientation=best_rotation, qc_flag=quality_flag, poly=best_poly
    )
end

# Raster based assessment methods

"""
    assess_reef_site(
        rst::Union{Raster,RasterStack},
        geom::GI.Wrappers.Polygon,
        ruleset::Dict{Symbol,Function};
        degree_step::Float64=15.0,
        start_rot::Float64=0.0,
        n_per_side::Int64=1
    )::Tuple{Float64,Int64,GI.Wrappers.Polygon}

Assess given reef site for it's suitability score at different specified rotations around the
initial reef-edge rotation.

# Arguments
- `rst` : Raster or RasterStack object used to assess site suitability.
- `geom` : Initial site polygon with no rotation applied.
- `ruleset` : Criteria ruleset to apply to `rst` pixels when assessing which pixels are suitable.
- `degree_step` : Degree value to vary each rotation by. Default = 15 degrees.
- `start_rot` : Initial rotation used to align the site polygon with the nearest reef edge. Default = 0 degrees.
- `n_per_side` : Number of times to rotate polygon on each side (clockwise and anticlockwise). Default = 2 rotations on each side.

# Returns
- Highest score identified with rotating polygons.
- The index of the highest scoring rotation.
- The polygon with the highest score out of the assessed rotated polygons.
"""
function assess_reef_site(
    rst::Union{Raster,RasterStack},
    geom::GI.Wrappers.Polygon,
    ruleset::Dict{Symbol,Function};
    degree_step::Float64=15.0,
    start_rot::Float64=0.0,
    n_per_side::Int64=2
)::Tuple{Float64,Int64,GI.Wrappers.Polygon}
    rotations =
        (start_rot - (degree_step * n_per_side)):degree_step:(start_rot + (degree_step * n_per_side))
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

    return score[argmax(score)], argmax(score) - (n_per_side + 1), best_poly[argmax(score)]
end

"""
    identify_edge_aligned_sites(
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
function identify_edge_aligned_sites(
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
    for (i, index) in enumerate(eachrow(search_pixels))
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

    return DataFrame(; score=best_score, orientation=best_rotation, poly=best_poly)
end
