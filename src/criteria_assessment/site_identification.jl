"""Methods for identifying potential deployment locations."""

"""
    proportion_suitable(subsection::BitMatrix, square_offset::Tuple=(-4,5))::Matrix{Int16}

Calculate the the proportion of the subsection that is suitable for deployments.
The `subsection` is the surrounding a rough hectare area centred on each cell of a raster
marked as being suitable according to user-selected criteria.

Cells on the edges of a raster object are assessed using a smaller surrounding area, rather
than shifting the window inward. In usual applications, there will be no target pixel
close to the edge due to the use of buffer areas.

# Arguments
- `x` : Matrix of boolean pixels after filtering with user criteria.
- `square_offset` : The number of pixels +/- around a center "target" pixel to assess as the
                    moving window. Defaults to (-4, 5).
                    Assuming a 10m² pixel, the default `square_offset` resolves to a one
                    hectare area.

# Returns
Matrix of values 0 - 100 indicating the percentage of the area around the target pixel that
meet suitability criteria.
"""
function proportion_suitable(x::BitMatrix; square_offset::Tuple=(-4, 5))::Matrix{Int16}
    subsection_dims = size(x)
    target_area = zeros(Int16, subsection_dims)

    @floop for row_col in ThreadsX.findall(x)
        (row, col) = Tuple(row_col)
        x_left = max(col + square_offset[1], 1)
        x_right = min(col + square_offset[2], subsection_dims[2])

        y_top = max(row + square_offset[1], 1)
        y_bottom = min(row + square_offset[2], subsection_dims[1])

        target_area[row, col] = Int16(sum(@views x[y_top:y_bottom, x_left:x_right]))
    end

    return target_area
end

"""
    filter_distances(
        target_rast::Raster,
        gdf::DataFrame,
        dist_nm
    )::Raster

Exclude pixels in `target_rast` that are beyond `dist_nm` (nautical miles) from a geometry
in `gdf`.  `target_rast` and `gdf` should be in the same CRS (EPSG:7844 / GDA2020 for GBR-reef-guidance-assessment).

# Arguments
- `target_rast` : Boolean raster of suitable pixels to filter.
- `gdf` : DataFrame with `geometry` column that contains vector objects of interest.
- `dist_nm` : Filtering distance from geometry object in nautical miles.

# Returns
- `tmp_areas` : Raster of filtered pixels containing only pixels within target distance
from a geometry centroid.
"""
function filter_distances(
    target_rast::Raster, gdf::DataFrame, dist; units::String="NM"
)::Raster
    tmp_areas = copy(target_rast)

    # First dimension is the rows (latitude)
    # Second dimension is the cols (longitude)
    raster_lat = Vector{Float64}(tmp_areas.dims[1].val)
    raster_lon = Vector{Float64}(tmp_areas.dims[2].val)

    @floop for row_col in ThreadsX.findall(tmp_areas)
        (lat_ind, lon_ind) = Tuple(row_col)
        point = AG.createpoint()

        lon = raster_lon[lon_ind]
        lat = raster_lat[lat_ind]
        AG.addpoint!(point, lon, lat)

        pixel_dists = AG.distance.([point], port_locs.geometry)
        geom_point = gdf[argmin(pixel_dists), :geometry]
        geom_point = (AG.getx(geom_point, 0), AG.gety(geom_point, 0))

        dist_nearest = Distances.haversine(geom_point, (lon, lat))

        # Convert from meters to nautical miles
        if units == "NM"
            dist_nearest = dist_nearest / 1852
        end

        # Convert from meters to kilometers
        if units == "km"
            dist_nearest = dist_nearest / 1000
        end

        tmp_areas.data[lon_ind, lat_ind] = dist_nearest < dist ? 1 : 0
    end

    return tmp_areas
end

"""
    mask_region(reg_assess_data, reg, qp, rtype)

# Arguments
- `reg_assess_data` : Regional assessment data
- `reg` : The region name to assess
- `qp` : query parameters
- `rtype` : region type (one of `:slopes` or `:flats`)

# Returns
Raster of region with locations that meet criteria masked.
"""
function mask_region(reg_assess_data, reg, qp, rtype)
    criteria_names, lbs, ubs = remove_rugosity(reg, parse_criteria_query(qp)...)

    # Otherwise, create the file
    @debug "$(now()) : Assessing criteria"
    assess = reg_assess_data[reg]
    mask_data = threshold_mask(
        assess,
        Symbol(rtype),
        CriteriaBounds.(criteria_names, lbs, ubs)
    )

    return mask_data
end

"""
    assess_region(reg_assess_data, reg, qp, rtype)

Perform raster suitability assessment based on user defined criteria.

# Arguments
- `reg_assess_data` : Dictionary containing the regional data paths, reef outlines and full region names.
- `reg` : Name of the region being assessed (format `Cairns-Cooktown` rather than `Cairns/Cooktown Management Area`).
- `qp` : Dict containing bounds for each variable being filtered.
- `rtype` : Type of zone to assess (flats or slopes).

# Returns
GeoTiff file of surrounding hectare suitability (1-100%) based on the criteria bounds input
by a user.
"""
function assess_region(reg_assess_data, reg, qp, rtype)
    @debug "Assessing region's suitability score"

    # Make mask of suitable locations
    mask_data = mask_region(reg_assess_data, reg, qp, rtype)

    # Assess remaining pixels for their suitability
    @debug "Calculating proportional suitability score"
    suitability_scores = proportion_suitable(mask_data.data)

    return rebuild(mask_data, suitability_scores)
end

"""
    assess_sites(
        reg_assess_data::OrderedDict,
        reg::String,
        rtype::String,
        pixel_criteria::Dict,
        site_criteria::Dict,
        assess_locs::Raster
    )

# Arguments
- `reg_assess_data` : Regional assessment data
- `reg` : Short region name
- `rtype` : Slopes or Flats assessment type
- `pixel_criteria` : parameters to assess specific locations with
- `site_criteria` : parameters to assess sites based on their polygonal representation
- `assess_locs` : Raster of suitability scores for each valid pixel

# Returns
GeoDataFrame of all potential sites
"""
function assess_sites(
    reg_assess_data::OrderedDict,
    reg::String,
    rtype::String,
    pixel_criteria::Dict,
    site_criteria::Dict,
    assess_locs::Raster
)
    criteria_names, lbs, ubs = remove_rugosity(reg, parse_criteria_query(pixel_criteria)...)

    # Otherwise, create the file
    @debug "$(now()) : Assessing criteria table"
    assess = reg_assess_data[reg]
    crit_pixels = apply_criteria_lookup(
        assess,
        Symbol(rtype),
        CriteriaBounds.(criteria_names, lbs, ubs)
    )

    res = abs(step(dims(assess_locs, X)))
    target_crs = convert(EPSG, crs(assess_locs))

    suitability_threshold = parse(Int64, (site_criteria["SuitabilityThreshold"]))
    assess_locs = identify_search_pixels(assess_locs, x -> x .> suitability_threshold)

    # Need reef outlines to indicate direction of the reef edge
    gdf = REGIONAL_DATA["reef_outlines"]
    reef_outlines = buffer_simplify(gdf)
    reef_outlines = polygon_to_lines.(reef_outlines)

    x_dist = parse(Int64, site_criteria["xdist"])
    y_dist = parse(Int64, site_criteria["ydist"])
    @debug "Assessing site polygons for $(size(assess_locs, 1)) locations."
    initial_polygons = identify_edge_aligned_sites(
        crit_pixels,
        assess_locs,
        res,
        gdf,
        x_dist,
        y_dist,
        target_crs,
        reef_outlines,
        reg
    )

    return initial_polygons
end
