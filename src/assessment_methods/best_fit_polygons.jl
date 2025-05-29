"""Geometry-based assessment methods."""

using NearestNeighbors

# Tabular data assessment methods

"""
    assess_reef_site(
        rel_pix::DataFrame,
        geom::GI.Wrappers.Polygon,
        max_count::Float64,
        target_crs::GeoFormatTypes.EPSG;
        degree_step::Float64=15.0,
        start_rot::Float64=0.0,
        n_per_side::Int64=2,
        suit_threshold::Float64=0.7
    )::Tuple{Float64,Int64,GI.Wrappers.Polygon,Int64}

Assesses the rotations of a search box `geom` for their suitability score (calculated as the
proportion of pixels that meet all specified criteria thresholds). Search box rotation steps
are returned so that the `start_rot` angle is 0, rotations anti-clockwise are negative and
rotations clockwise are positive.

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
    target_crs::GeoFormatTypes.EPSG;
    degree_step::Float64=15.0,
    start_rot::Float64=0.0,
    n_per_side::Int64=1,
    suit_threshold::Float64=0.7
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
    max_count::Float64,
    n_per_side::Int64,
    suit_threshold::Float64
)::Tuple{Float64,Int64,GI.Wrappers.Polygon,Int64}
    # Implementation with pre-rotations
    n_rotations = length(rotated)
    score = @MVector zeros(n_rotations)
    qc_flag = @MVector zeros(Int64, n_rotations)

    for (j, r) in enumerate(rotated)
        score[j] = count(GO.intersects.([r], rel_pix.geometry)) / max_count

        if score[j] < suit_threshold
            # Early exit as there's no point in searching further.
            # Changing the rotation is unlikely to improve the score.
            qc_flag[j] = 1
            break
        end
    end

    return min(score[argmax(score)], 1),
    argmax(score) - (n_per_side + 1),
    rotated[argmax(score)],
    maximum(qc_flag)
end

"""
    assess_reef_site(
        rel_pix::DataFrame,
        rotated::Vector{GI.Wrappers.Polygon},
        suitability_threshold::Float64
    )::Tuple{Float64,Int64,GI.Wrappers.Polygon,Int64}

Assesses the rotations of a search box `geom` for their suitability score (calculated as the
proportion of pixels that meet all specified criteria thresholds).

# Extended help
The scores produced are a proportion of the polygon that are covered by valid pixel points,
relative to the maximum number of points (`max_count`) (0-1). `max_count` is approximate
(determined by user `x_dist` and `y_dist` box dimensions) and doesn't account for buffering
and rotation of the search box. In rare cases scores could be > 1, however returned values
are capped at max 1.

# Arguments
- `rel_pix` : The point data for relevant pixels that are within the search area around a pixel.
- `rotated` : Pre-rotated geometries.
- `suitability_threshold` : Suitability threshold, below which sites are excluded from result sets.

# Returns
- Highest score
- Highest scoring rotation step
- Highest scoring polygon
- Quality control flag for site, indicating if `suitability_threshold` was breached.
"""
function assess_reef_site(
    rel_pix::DataFrame,
    rotated::Vector{GI.Wrappers.Polygon},
    suitability_threshold::Float64
)::Tuple{Float64,Int64,GI.Wrappers.Polygon,Int64}
    # Implementation with pre-rotations
    n_rotations = length(rotated)
    score = @MVector zeros(n_rotations)
    qc_flag = @MVector zeros(Int64, n_rotations)

    for (j, r) in enumerate(rotated)
        score[j] = count(GO.contains.(Ref(r), rel_pix.geometry))

        if score[j] < suitability_threshold
            # Early exit as there's no point in searching further.
            # Changing the rotation is unlikely to improve the score.
            qc_flag[j] = 1
            break
        end

        if score[j] >= (nrow(rel_pix) * 0.98)
            # Found near the best possible score so might as well exit
            break
        end
    end

    best_idx = argmax(score)
    return (
        score[best_idx],
        best_idx,
        rotated[best_idx],
        qc_flag[best_idx]
    )
end

"""
    identify_edge_aligned_sites(
        env_lookup::DataFrame,
        search_pixels::DataFrame,
        res::Float64,
        gdf::DataFrame,
        x_dist::Union{Int64,Float64},
        y_dist::Union{Int64,Float64},
        target_crs::GeoFormatTypes.EPSG,
        region::String;
        degree_step::Float64=15.0,
        n_rot_p_side::Int64=2,
        suit_threshold::Float64=0.7
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
- `suit_threshold` : Theshold used to skip searching where the proportion of suitable pixels is too low.

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
    target_crs::GeoFormatTypes.EPSG,
    region::String;
    degree_step::Float64=15.0,
    n_rot_p_side::Int64=1,
    suit_threshold::Float64=0.7
)::DataFrame
    region_long = REGIONAL_DATA["region_long_names"][region]
    target_reefs = reef_outlines[reef_outlines.management_area .== region_long, :]
    max_count = (
        (x_dist * y_dist) / ceil(degrees_to_meters(res, mean(search_pixels.lats)))^2
    )

    reef_lines = polygon_to_lines.(buffer_simplify(target_reefs))
    search_box = initial_search_box(
        (search_pixels.lons[1], search_pixels.lats[1]),
        x_dist,
        y_dist,
        target_crs
    )
    start_deg_angle = initial_search_rotation(
        GO.Point(search_pixels.lons[1], search_pixels.lats[1]),
        search_box,
        target_reefs[!, first(GI.geometrycolumns(target_reefs))],
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

    # Create KD-tree to identify pixels for assessment that are close to each other.
    # We can then apply a heuristic to avoid near-identical assessments.
    sp_inds = Tuple.(search_pixels.indices)
    inds = Matrix{Float32}([first.(sp_inds) last.(sp_inds)]')
    kdtree = KDTree(inds; leafsize=25)
    ignore_idx = Int64[]
    for (i, coords) in enumerate(eachcol(inds))
        if i in ignore_idx
            continue
        end

        # If there are a group of pixels close to each other, only assess the one closest to
        # the center of the group.
        # Select pixels within ~22 meters (1 pixel = 10m² ∴ 3.3 ≈ 33m horizontal distance)
        idx, dists = knn(kdtree, coords, 30)  # retrieve 30 closest locations
        sel = (dists == 0) .| (dists .< 3.3)  # select current pixel and those ~33m away
        xs = mean(inds[1, idx[sel]])
        ys = mean(inds[2, idx[sel]])

        # Find index of pixel closest to the center of the group
        closest_idx = argmin(
            map(p2 -> sum((p2 .- (xs, ys)) .^ 2), eachcol(inds[:, idx[sel]]))
        )

        to_keep = idx[sel][closest_idx]
        to_ignore = idx[sel][idx[sel] .!= to_keep]

        append!(ignore_idx, to_ignore)
    end

    # Search each location to assess
    ignore_idx = unique(ignore_idx)
    assessment_locs = search_pixels[Not(ignore_idx), :]
    n_pixels = nrow(assessment_locs)
    @debug "KD-tree filtering : removed $(length(ignore_idx)) near-identical locations, assessing $(n_pixels) locations"

    best_score = zeros(n_pixels)
    best_poly = Vector(undef, n_pixels)
    best_rotation = zeros(Int64, n_pixels)
    quality_flag = zeros(Int64, n_pixels)
    FLoops.assistant(false)
    @floop for (i, pix) in enumerate(eachrow(assessment_locs))
        lon = pix.lons
        lat = pix.lats

        rotated_copy::Vector{GI.Wrappers.Polygon} = move_geom.(rotated_geoms, [(lon, lat)])

        max_offset = (
            abs(meters_to_degrees(maximum([x_dist, y_dist]) / 2, lat)) +
            (2 * res)
        )
        bounds = Float64[
            lon - max_offset,
            lon + max_offset,
            lat - max_offset,
            lat + max_offset
        ]

        # Identify relevant pixels to assess
        loc_constraint = (
            (env_lookup.lons .>= bounds[1]) .&  # left
            (env_lookup.lons .<= bounds[2]) .&  # right
            (env_lookup.lats .>= bounds[3]) .&  # bottom
            (env_lookup.lats .<= bounds[4])     # top
        )
        rel_pix = env_lookup[loc_constraint, :]

        # Skip if no relevant pixels or if the number of suitable pixels are below required
        # threshold
        if nrow(rel_pix) == 0 || ((nrow(rel_pix) / max_count) < suit_threshold)
            best_score[i] = 0.0
            best_rotation[i] = 0
            best_poly[i] = rotated_copy[1]
            quality_flag[i] = 1
            continue
        end

        best_score[i], best_rotation[i], best_poly[i], quality_flag[i] = assess_reef_site(
            rel_pix,
            rotated_copy,
            max_count,
            n_rot_p_side,
            suit_threshold
        )
    end

    return DataFrame(;
        score=best_score,
        orientation=best_rotation,
        qc_flag=quality_flag,
        geometry=best_poly
    )
end

"""
    find_optimal_site_alignment(
        env_lookup::DataFrame,
        search_pixels::DataFrame,
        res::Float64,
        x_dist::Union{Int64,Float64},
        y_dist::Union{Int64,Float64},
        target_crs::GeoFormatTypes.EPSG;
        degree_step::Float64=20.0,
        suit_threshold::Float64=0.7
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
- `x_dist` : Length of horizontal side of search box (in meters).
- `y_dist` : Length of vertical side of search box (in meters).
- `target_crs` : CRS of the input Rasters. Using GeoFormatTypes.EPSG().
- `degree_step` : Degree to perform rotations around identified edge angle.
- `suit_threshold` : Theshold used to skip searching where the proportion of suitable pixels is too low.

# Returns
DataFrame containing highest score, rotation and polygon for each assessment at pixels in indices.
"""
function find_optimal_site_alignment(
    env_lookup::DataFrame,
    search_pixels::DataFrame,
    res::Float64,
    x_dist::Union{Int64,Float64},
    y_dist::Union{Int64,Float64},
    target_crs::GeoFormatTypes.EPSG;
    degree_step::Float64=10.0,
    suit_threshold::Float64=0.7
)::DataFrame
    max_count = (
        (x_dist * y_dist) / ceil(degrees_to_meters(res, mean(search_pixels.lats)))^2
    )

    search_box = initial_search_box(
        (search_pixels.lons[1], search_pixels.lats[1]),
        x_dist,
        y_dist,
        target_crs
    )

    # Precalculate rotations
    rotations = 0.0:degree_step:179.0
    rotated_geoms = Vector{GI.Wrappers.Polygon}(undef, length(rotations))
    for (j, r) in enumerate(rotations)
        rotated_geoms[j] = rotate_polygon(search_box, r)
    end

    # Create KD-tree to identify pixels for assessment that are close to each other.
    # We can then apply a heuristic to avoid near-identical assessments.
    sp_inds = Tuple.(search_pixels.indices)
    inds = Matrix{Float32}([first.(sp_inds) last.(sp_inds)]')
    kdtree = KDTree(inds; leafsize=25)
    ignore_idx = Int64[]
    for (i, coords) in enumerate(eachcol(inds))
        if i in ignore_idx
            continue
        end

        # If there are a group of pixels close to each other, only assess the one closest to
        # the center of the group.
        # Select pixels within ~22 meters (1 pixel = 10m² ∴ 3.3 ≈ 33m horizontal distance)
        idx, dists = knn(kdtree, coords, 30)  # retrieve 30 closest locations
        sel = (dists == 0) .| (dists .< 3.3)  # select current pixel and those ~33m away
        xs = mean(inds[1, idx[sel]])
        ys = mean(inds[2, idx[sel]])

        # Find index of pixel closest to the center of the group
        closest_idx = argmin(
            map(p2 -> sum((p2 .- (xs, ys)) .^ 2), eachcol(inds[:, idx[sel]]))
        )

        to_keep = idx[sel][closest_idx]
        to_ignore = idx[sel][idx[sel] .!= to_keep]

        append!(ignore_idx, to_ignore)
    end

    # Search each location to assess
    ignore_idx = unique(ignore_idx)
    assessment_locs = search_pixels[Not(ignore_idx), :]
    n_pixels = nrow(assessment_locs)
    @debug "$(now()) : KD-tree filtering - removed $(length(ignore_idx)) near-identical locations, now assessing $(n_pixels) locations"

    best_score = zeros(n_pixels)
    best_poly = Vector(undef, n_pixels)
    best_rotation = zeros(Int64, n_pixels)
    quality_flag = zeros(Int64, n_pixels)
    max_x_y_dist = maximum([x_dist, y_dist])
    FLoops.assistant(false)
    @floop for (i, pix) in enumerate(eachrow(assessment_locs))
        lon = pix.lons
        lat = pix.lats

        rotated_copy::Vector{GI.Wrappers.Polygon} = move_geom.(
            rotated_geoms,
            Ref((lon, lat))
        )

        max_offset = (
            abs(meters_to_degrees(max_x_y_dist * 0.5, lat)) +
            (2.0 * res)
        )

        bounds = Float64[
            lon - max_offset,
            lon + max_offset,
            lat - max_offset,
            lat + max_offset
        ]

        # Identify relevant pixels to assess
        loc_constraint = (
            (env_lookup.lons .>= bounds[1]) .&  # left
            (env_lookup.lons .<= bounds[2]) .&  # right
            (env_lookup.lats .>= bounds[3]) .&  # bottom
            (env_lookup.lats .<= bounds[4])     # top
        )
        rel_pix = env_lookup[loc_constraint, :]
        n_matches = nrow(rel_pix)

        # Skip if no relevant pixels or if the number of suitable pixels are below required
        # threshold
        if n_matches == 0 || ((n_matches / max_count) < suit_threshold)
            best_score[i] = 0.0
            best_rotation[i] = 0
            best_poly[i] = rotated_copy[1]
            quality_flag[i] = 1
            continue
        end

        best_score[i], best_rotation[i], best_poly[i], quality_flag[i] = assess_reef_site(
            rel_pix,
            rotated_copy,
            10.0
        )
    end

    if maximum(best_score) > 0.0
        best_score = min.(best_score ./ max_count, 1.0)
    end

    return DataFrame(;
        score=best_score,
        orientation=best_rotation,
        qc_flag=quality_flag,
        geometry=best_poly
    )
end

"""
    assess_reef_site(
        rst::Union{Raster,RasterStack},
        geom::GI.Wrappers.Polygon,
        ruleset::Dict{Symbol,Function},
        degree_step::Float64,
        start_rot::Float64
    )::Tuple{Float64,Int64,GI.Wrappers.Polygon}

Assess given reef site for it's suitability score at different specified rotations.

# Arguments
- `rst` : Raster of suitability scores.
- `geom` : Initial site polygon with no rotation applied.
- `ruleset` : Criteria ruleset to apply to `rst` pixels when assessing which pixels are suitable.
- `degree_step` : Degree value to vary each rotation by. Default = 20 degrees.
- `start_rot` : Initial rotation used to align the site polygon with the nearest reef edge. Default = 0 degrees.

# Returns
- Highest score identified with rotating polygons.
- The index of the highest scoring rotation.
- The polygon with the highest score out of the assessed rotated polygons.
"""
function assess_reef_site(
    rst::Union{Raster,RasterStack},
    geom::GI.Wrappers.Polygon,
    degree_step::Float64,
    start_rot::Float64
)::Tuple{Float64,Int64,GI.Wrappers.Polygon}
    rotations = start_rot:degree_step:360.0
    n_rotations = length(rotations)
    score = zeros(n_rotations)
    best_poly = Vector(undef, n_rotations)

    target_crs = convert(EPSG, GI.crs(rst))
    for (j, r) in enumerate(rotations)
        rot_geom = rotate_geom(geom, r, target_crs)
        c_rst = crop(rst; to=rot_geom)
        if !all(size(c_rst) .> (0, 0))
            @warn "No data found!"
            continue
        end

        score[j] = mean(c_rst)
        best_poly[j] = rot_geom

        if score[j] > 0.95
            # Found a good enough rotation
            break
        end
    end

    return score[argmax(score)], argmax(score) - (n_per_side + 1), best_poly[argmax(score)]
end

"""
    identify_edge_aligned_sites(
        rst_stack::Raster,
        search_pixels::DataFrame,
        gdf::DataFrame,
        x_dist::Union{Int64,Float64},
        y_dist::Union{Int64,Float64},
        target_crs::GeoFormatTypes.EPSG,
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
- `rst_stack` : Raster containing environmental variables for assessment.
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
    rst_stack::Raster,
    search_pixels::DataFrame,
    gdf::DataFrame,
    x_dist::Union{Int64,Float64},
    y_dist::Union{Int64,Float64},
    target_crs::GeoFormatTypes.EPSG,
    region::String;
    degree_step::Float64=15.0,
    n_rot_per_side::Int64=1
)::DataFrame
    gdf = gdf[gdf.management_area_short .== region, :]

    reef_lines = polygon_to_lines.(buffer_simplify(gdf))

    res = abs(step(dims(rst_stack, X)))

    # Search each location to assess
    best_score = zeros(length(search_pixels.lons))
    best_poly = Vector(undef, length(search_pixels.lons))
    best_rotation = zeros(Int64, length(search_pixels.lons))
    for (i, index) in enumerate(eachrow(search_pixels))
        lon = index.lons
        lat = index.lats
        geom_buff = initial_search_box((lon, lat), x_dist, y_dist, target_crs)

        pixel = GO.Point(lon, lat)
        rot_angle = initial_search_rotation(pixel, geom_buff, gdf.geometry, reef_lines)

        b_score, b_rot, b_poly = assess_reef_site(
            rst_stack,
            geom_buff,
            degree_step,
            rot_angle
        )

        best_score[i] = b_score
        best_rotation[i] = b_rot
        best_poly[i] = b_poly
    end

    return DataFrame(; score=best_score, orientation=best_rotation, geometry=best_poly)
end
