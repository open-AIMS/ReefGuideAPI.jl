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
        suit_threshold::Float64=0.3
    )::Tuple{Float64,Int64,GI.Wrappers.Polygon,Int64}

Assesses the rotations of a search box `geom` for their suitability score (calculated as the
proportion of pixels that meet all specified criteria thresholds). Search box rotation steps
are returned so that the `start_rot` angle is 0, rotations anti-clockwise are negative and
rotations clockwise are
positive.

# Extended help
The scores produced are a proportion of the polygon that are covered by valid pixel points,
relative to the maximum number of points (`max_count`) (0-1). `max_count` is approximate
(determined by user `x_dist` and `y_dist` box dimensions) and doesn't account for buffering
and rotation of the search box. In rare cases scores could be > 1, however returned values
are capped at max 1.

# Arguments
- `rel_pix` : The point data for relevant pixels that are within the search area around a pixel.
- `geom` : Starting search box for assessment.
- `max_count` : The maximum number of pixels that can intersect the search box (used to standardise scores between 0 and 1).
- `target_crs` : Coordinate Reference System used for analysis vector and raster data.
- `degree_step` : Step to vary the search box rotations.
- `start_rot` : Starting angle rotation that aligns the box with the closest reef edge.
- `n_per_side` : Number of rotations to perform around the starting search box angle.
- `suit_threshold` : Suitability threshold, below which sites are excluded from result sets.

# Returns
- Highest score
- Highest scoring rotation step
- Highest scoring polygon
- Quality control flag for site, indicating if `suit_threshold` was met in the highest scoring rotation.
"""
function assess_reef_site(
    rel_pix::DataFrame,
    geom::GI.Wrappers.Polygon,
    max_count::Float64,
    target_crs::GeoFormatTypes.CoordinateReferenceSystemFormat;
    degree_step::Float64=15.0,
    start_rot::Float64=0.0,
    n_per_side::Int64=2,
    suit_threshold::Float64=0.3
)::Tuple{Float64,Int64,GI.Wrappers.Polygon,Int64}
    rot_start = (start_rot - (degree_step * n_per_side))
    rot_end = (start_rot + (degree_step * n_per_side))
    rotations = rot_start:degree_step:rot_end
    n_rotations = length(rotations)
    score = zeros(n_rotations)
    best_poly = Vector(undef, n_rotations)
    qc_flag = zeros(Int64, n_rotations)

    for (j, r) in enumerate(rotations)
        rot_geom = rotate_geom(geom, r, target_crs)
        score[j] =
            length(
                rel_pix[
                    GO.intersects.([rot_geom], rel_pix.geometry), :lon_idx
                ]
            ) / max_count
        best_poly[j] = rot_geom

        if score[j] < suit_threshold
            # Early exit as there's no point in searching further.
            # Changing the rotation is unlikely to improve the score.
            qc_flag[j] = 1
            break
        end
    end

    return min(score[argmax(score)], 1),
    argmax(score) - (n_per_side + 1),
    best_poly[argmax(score)],
    maximum(qc_flag)
end
function assess_reef_site(
    rel_pix::DataFrame,
    rotated::Vector{GI.Wrappers.Polygon},
    max_count::Float64;
    n_per_side::Int64=2,
    suit_threshold::Float64=0.3
)::Tuple{Float64,Int64,GI.Wrappers.Polygon,Int64}
    # Implementation with pre-rotations
    n_rotations = length(rotated)
    score = @MVector zeros(n_rotations)
    best_poly = Vector(undef, n_rotations)
    qc_flag = @MVector zeros(Int64, n_rotations)

    @floop for (j, r) in enumerate(rotated)
        score[j] =
            length(
                rel_pix[
                    GO.intersects.([r], rel_pix.geometry), :lon_idx
                ]
            ) / max_count
        best_poly[j] = r

        if score[j] < suit_threshold
            # Early exit as there's no point in searching further.
            # Changing the rotation is unlikely to improve the score.
            qc_flag[j] = 1
            break
        end
    end

    return min(score[argmax(score)], 1),
    argmax(score) - (n_per_side + 1),
    best_poly[argmax(score)],
    maximum(qc_flag)
end

"""
    identify_edge_aligned_sites(
        env_lookup::DataFrame,
        search_pixels::DataFrame,
        res::Float64,
        gdf::DataFrame,
        x_dist::Union{Int64,Float64},
        y_dist::Union{Int64,Float64},
        target_crs::GeoFormatTypes.CoordinateReferenceSystemFormat,
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
- `env_lookup` : DataFrame containing environmental variables for assessment.
- `search_pixels` : DataFrame containing lon and lat columns for each pixel that is intended for analysis.
- `res` : Resolution of the original raster pixels. Can by found via `abs(step(dims(raster, X)))`.
- `reef_outlines` : GeoDataFrame of reef outlines used to align the search box edge.
- `x_dist` : Length of horizontal side of search box (in meters).
- `y_dist` : Length of vertical side of search box (in meters).
- `target_crs` : CRS of the input Rasters. Using GeoFormatTypes.EPSG().
- `region` : Region name, e.g. "Cairns-Cooktown" or "FarNorthern".
- `degree_step` : Degree to perform rotations around identified edge angle.
- `n_rot_p_side` : Number of rotations to perform clockwise and anticlockwise around the identified edge angle. Default 2 rotations.
- `surr_threshold` : Theshold used to skip searching where the proportion of suitable pixels is too low.

# Returns
DataFrame containing highest score, rotation and polygon for each assessment at pixels in indices.
"""
function identify_edge_aligned_sites(
    env_lookup::DataFrame,
    search_pixels::DataFrame,
    res::Float64,
    reef_outlines::DataFrame,
    x_dist::Union{Int64,Float64},
    y_dist::Union{Int64,Float64},
    target_crs::GeoFormatTypes.CoordinateReferenceSystemFormat,
    region::String;
    degree_step::Float64=15.0,
    n_rot_p_side::Int64=2,
    suit_threshold::Float64=0.33
)::DataFrame
    region_long = REGIONAL_DATA["region_long_names"][region]
    target_reefs = reef_outlines[reef_outlines.management_area .== region_long, :]
    max_count = (
        (x_dist * y_dist) / (degrees_to_meters(res, mean(search_pixels.lats))^2)
    )

    reef_lines = polygon_to_lines.(buffer_simplify(target_reefs))

    # Search each location to assess
    n_pixels = nrow(search_pixels)
    best_score = zeros(n_pixels)
    best_poly = Vector(undef, n_pixels)
    best_rotation = zeros(Int64, n_pixels)
    quality_flag = zeros(Int64, n_pixels)
    bounds = zeros(4)
    target_geoms = target_reefs[!, first(GI.geometrycolumns(target_reefs))]
    loc_constraint = BitVector(falses(nrow(env_lookup)))

    search_box = initial_search_box(
        (search_pixels.lons[1], search_pixels.lats[1]),
        x_dist,
        y_dist,
        target_crs,
        res
    )
    start_deg_angle = initial_search_rotation(
        GO.Point(search_pixels.lons[1], search_pixels.lats[1]),
        search_box,
        target_geoms,
        reef_lines
    )
    search_box = rotate_polygon(search_box, start_deg_angle)

    # Precalculate rotations
    rot_start = (start_deg_angle - (degree_step * n_rot_p_side))
    rot_end = (start_deg_angle + (degree_step * n_rot_p_side))
    rotations = rot_start:degree_step:rot_end
    rotated_geoms = Vector{GI.Wrappers.Polygon}(undef, length(rotations))
    for (j, r) in enumerate(rotations)
        rotated_geoms[j] = rotate_polygon(search_box, r)
    end

    for (i, pix) in enumerate(eachrow(search_pixels))
        lon = pix.lons
        lat = pix.lats

        rotated_geoms .= move_geom.(rotated_geoms, [(lon, lat)])

        max_offset = (
            abs(meters_to_degrees(maximum([x_dist, y_dist]) / 2, lat)) +
            (2 * res)
        )
        bounds .= Float64[
            lon - max_offset,
            lon + max_offset,
            lat - max_offset,
            lat + max_offset
        ]

        # Identify relevant pixels to assess
        loc_constraint .= (
            (env_lookup.lons .>= bounds[1]) .&  # left
            (env_lookup.lons .<= bounds[2]) .&  # right
            (env_lookup.lats .>= bounds[3]) .&  # bottom
            (env_lookup.lats .<= bounds[4])     # top
        )
        rel_pix = env_lookup[loc_constraint, :]

        if nrow(rel_pix) == 0
            best_score[i] = 0.0
            best_rotation[i] = 0
            best_poly[i] = geom_buff
            quality_flag[i] = 1
            continue
        end

        best_score[i], best_rotation[i], best_poly[i], quality_flag[i] = assess_reef_site(
            rel_pix,
            rotated_geoms,
            max_count;
            n_per_side=n_rot_p_side,
            suit_threshold=suit_threshold
        )
    end

    return DataFrame(;
        score=best_score,
        orientation=best_rotation,
        qc_flag=quality_flag,
        geometry=best_poly
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
- `rst` : Raster of suitability scores.
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
    geom::GI.Wrappers.Polygon;
    degree_step::Float64=15.0,
    start_rot::Float64=0.0,
    n_per_side::Int64=2
)::Tuple{Float64,Int64,GI.Wrappers.Polygon}
    rotations =
        (start_rot - (degree_step * n_per_side)):degree_step:(start_rot + (degree_step * n_per_side))
    n_rotations = length(rotations)
    score = zeros(n_rotations)
    best_poly = Vector(undef, n_rotations)

    target_crs = GI.crs(rst)
    for (j, r) in enumerate(rotations)
        rot_geom = rotate_geom(geom, r, target_crs)
        c_rst = crop(rst; to=rot_geom)
        if !all(size(c_rst) .> (0, 0))
            @warn "No data found!"
            continue
        end

        score[j] = mean(c_rst)
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
    region::String;
    degree_step::Float64=15.0,
    n_rot_per_side::Int64=2
)::DataFrame
    gdf = gdf[gdf.management_area .== region, :]
    reef_lines = polygon_to_lines.(buffer_simplify(gdf))

    res = abs(step(dims(rst_stack, X)))

    # Search each location to assess
    best_score = zeros(length(search_pixels.lons))
    best_poly = Vector(undef, length(search_pixels.lons))
    best_rotation = zeros(Int64, length(search_pixels.lons))
    for (i, index) in enumerate(eachrow(search_pixels))
        lon = index.lons
        lat = index.lats
        geom_buff = initial_search_box((lon, lat), x_dist, y_dist, target_crs, res)

        pixel = GO.Point(lon, lat)
        rot_angle = initial_search_rotation(pixel, geom_buff, gdf, reef_lines)

        b_score, b_rot, b_poly = assess_reef_site(
            rst_stack,
            geom_buff;
            degree_step=degree_step,
            start_rot=rot_angle,
            n_per_side=n_rot_per_side
        )

        best_score[i] = b_score
        best_rotation[i] = b_rot
        best_poly[i] = b_poly
    end

    return DataFrame(; score=best_score, orientation=best_rotation, geometry=best_poly)
end
