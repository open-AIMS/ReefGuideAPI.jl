"""Methods for identifying potential deployment locations."""

using FLoops, ThreadsX

"""
    proportion_suitable(subsection::BitMatrix, window::Tuple=(-4,5))::Matrix{Int16}

Calculate the the proportion of the subsection that is suitable for deployments.
Subsection is the surrounding a rough hectare area centred on each cell of a raster marked
as being suitable according to user-selected criteria.

# Arguments
- `x` : Matrix of boolean pixels after filtering with user criteria.
- `window` : Window size to assess. Default window (-4,5) assesses a square hectare around each target pixel where the resolution of pixels is 10m.
"""
function proportion_suitable(x::BitMatrix; window::Tuple=(-4, 5))::Matrix{Int16}
    x′ = zeros(Int16, size(x))

    @floop for row_col in ThreadsX.findall(x)
        (row, col) = Tuple(row_col)
        x_left = max(col + window[1], 1)
        x_right = min(col + window[2], size(x, 1))

        y_top = max(row + window[1], 1)
        y_bottom = min(row + window[2], size(x, 2))

        x′[row, col] = Int16(sum(@views x[y_top:y_bottom, x_left:x_right]))
    end

    return x′
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
- `target_rast` : Raster of suitable pixels (Bool) to filter pixels from.
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

    @floop for row_col in findall(tmp_areas)
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
# TODO: Better name, and address duplication in "/assess/{reg}/{rtype}"
"""
function _temp_assess_region(reg_assess_data, reg, qp, rtype, config)
    criteria_names, lbs, ubs = remove_rugosity(reg, parse_criteria_query(qp)...)

    # Otherwise, create the file
    @debug "$(now()) : Assessing criteria"
    assess = reg_assess_data[reg]
    mask_data = make_threshold_mask(
        assess,
        Symbol(rtype),
        CriteriaBounds.(criteria_names, lbs, ubs)
    )

    return mask_data
end

function _temp_filter_lookup(reg_assess_data, reg, qp, rtype, config)
    criteria_names, lbs, ubs = remove_rugosity(reg, parse_criteria_query(qp)...)

    # Otherwise, create the file
    @debug "$(now()) : Assessing criteria table"
    assess = reg_assess_data[reg]
    crit_lookup = apply_criteria_lookup(
        assess,
        Symbol(rtype),
        CriteriaBounds.(criteria_names, lbs, ubs)
    )

    return crit_lookup
end

"""
    assess_region(reg_assess_data, reg, qp, rtype, config)

Perform raster suitability assessment based on user defined criteria.

# Arguments
- `reg_assess_data` : Dictionary containing the regional data paths, reef outlines and full region names.
- `reg` : Name of the region being assessed (format `Cairns-Cooktown` rather than `Cairns/Cooktown Management Area`).
- `qp` : Dict containing bounds for each variable being filtered.
- `rtype` : Type of zone to assess (flats or slopes).
- `config` : Information from `.config.toml` file.

# Returns
GeoTiff file of surrounding hectare suitability (1-100%) based on the criteria bounds input
by a user.
"""
function assess_region(reg_assess_data, reg, qp, rtype, config)
    @debug "Assessing region's suitability score"

    file_id = string(hash(qp))
    assessed_tmp_path = _cache_location(config)
    assessed_path_tif = joinpath(assessed_tmp_path, file_id * "_suitable.tiff")

    if isfile(assessed_path_tif)
        return file(assessed_path_tif)
    end

    # Make mask of suitable locations
    mask_data = _temp_assess_region(reg_assess_data, reg, qp, rtype, config)

    # Assess remaining pixels for their suitability
    @debug "Calculating proportional suitability score"
    suitability_scores = proportion_suitable(mask_data.data)

    @debug "$(now()) : Running on thread $(threadid())"
    @debug "Writing to $(assessed_path)"
    Rasters.write(
        assessed_path,
        rebuild(mask_data, suitability_scores);
        ext=".tiff",
        source="gdal",
        driver="COG",
        options=Dict{String,String}(
            "COMPRESS" => "DEFLATE",
            "SPARSE_OK" => "TRUE",
            "OVERVIEW_COUNT" => "5",
            "BLOCKSIZE" => "256",
            "NUM_THREADS" => n_gdal_threads(config)
        ),
        force=true
    )

    return file(assessed_path)
end

"""
    site_assess_region(reg_assess_data, reg, criteria_qp, assessment_qp, rtype, config)

Perform site suitability assessment using polygon searches with user defined criteria.

# Arguments
- `reg_assess_data` : Dictionary containing the regional data paths, reef outlines and full region names.
- `reg` : Name of the region being assessed (format `Cairns-Cooktown` rather than `Cairns/Cooktown Management Area`).
- `criteria_qp` : Dict containing bounds for each variable being filtered.
- `assessment_qp` : Dict containing the dimensions of the search polygon (`xdist`, `ydist`) and the `SuitabilityThreshold` to identify search pixels.
- `rtype` : Type of zone to assess (flats or slopes).
- `config` : Information from `.config.toml` file.

# Returns
GeoJSON file containing result site polygons after filtering out polygons < 0.33 score and
keeping the highest scoring polygon where intersecting polygons occur.
"""
function site_assess_region(reg_assess_data, reg, criteria_qp, assessment_qp, rtype, config)
    @debug "Assessing region's suitability score"

    file_id = string(hash(assessment_qp))
    assessed_tmp_path = _cache_location(config)
    assessed_path_geojson = joinpath(
        assessed_tmp_path, file_id * "_potential_sites.geojson"
    )

    if isfile(assessed_path_geojson)
        return file(assessed_path_geojson)
    end

    file_id = string(hash(criteria_qp))
    assessed_tmp_path = _cache_location(config)
    assessed_path_tif = joinpath(assessed_tmp_path, file_id * "_suitable.tiff")

    if isfile(assessed_path_tif)
        scan_locs = Raster(assessed_path_tif)
    else
        # Make mask of suitable locations
        mask_data = _temp_assess_region(reg_assess_data, reg, criteria_qp, rtype, config)

        # Assess remaining pixels for their suitability
        @debug "Calculating proportional suitability score"
        suitability_scores = proportion_suitable(mask_data.data)

        # Need rebuild(mask_data, suitability_scores) as input to identify_potential_sites_edges
        scan_locs = rebuild(mask_data, suitability_scores)
    end

    # Need dataframe of valid_lookup pixels
    @debug "Pre-processing assessment inputs."
    crit_pixels = _temp_filter_lookup(reg_assess_data, reg, criteria_qp, rtype, config)

    res = abs(step(dims(scan_locs, X)))
    target_crs = convert(EPSG, crs(scan_locs))

    suitability_threshold = parse(Int64, (assessment_qp["SuitabilityThreshold"]))
    scan_locs = identify_search_pixels(scan_locs, x -> x .> suitability_threshold)

    # Need reef outlines
    gdf = REGIONAL_DATA["reef_outlines"]
    reef_outlines = buffer_simplify(gdf)
    reef_outlines = polygon_to_lines.(reef_outlines)

    x_dist = parse(Int64, assessment_qp["xdist"])
    y_dist = parse(Int64, assessment_qp["ydist"])

    @debug "Performing polygon site assessment."
    initial_polygons = identify_potential_sites_edges(
        crit_pixels,
        scan_locs,
        res,
        gdf,
        x_dist,
        y_dist,
        target_crs,
        reef_outlines,
        reg
    )
    output_geojson(assessed_path_geojson, filter_sites(initial_polygons))

    return file(assessed_path_geojson)
end
