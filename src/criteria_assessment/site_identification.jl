"""Methods for identifying potential deployment locations."""

"""
    proportion_suitable(
        x::Union{BitMatrix,SparseMatrixCSC{Bool,Int64}}; square_offset::Tuple=(-4, 5)
    )::SparseMatrixCSC{UInt8,Int64}

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
function proportion_suitable(
    x::Union{BitMatrix,SparseMatrixCSC{Bool,Int64}}; square_offset::Tuple=(-4, 5)
)::SparseMatrixCSC{UInt8,Int64}
    subsection_dims = size(x)
    target_area = spzeros(UInt8, subsection_dims)

    for row_col in findall(x)
        (row, col) = Tuple(row_col)
        x_left = max(col + square_offset[1], 1)
        x_right = min(col + square_offset[2], subsection_dims[2])

        y_top = max(row + square_offset[1], 1)
        y_bottom = min(row + square_offset[2], subsection_dims[1])

        target_area[row, col] = UInt8(sum(@views x[y_top:y_bottom, x_left:x_right]))
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
    lookup_assess_region(reg_assess_data, reg, qp, rtype; x_dist=100.0, y_dist=100.0)

Perform suitability assessment with the lookup table based on user-defined criteria.
This is currently orders of magnitude slower than the raster-based approach, although
it uses significantly less memory.

# Arguments
- `reg_assess_data` : Dictionary containing the regional data paths, reef outlines and full region names.
- `reg` : Name of the region being assessed (format `Cairns-Cooktown` rather than `Cairns/Cooktown Management Area`).
- `qp` : Dict containing bounds for each variable being filtered.
- `rtype` : Type of zone to assess (flats or slopes).
- `x_dist` : width of search polygon
- `y_dist` : height of search polygon

# Returns
Raster of surrounding hectare suitability (1-100%) based on the criteria bounds input
by a user.
"""
function lookup_assess_region(reg_assess_data, reg, qp, rtype; x_dist=100.0, y_dist=100.0)
    criteria_names, lbs, ubs = remove_rugosity(reg, parse_criteria_query(qp)...)

    assess_locs = getfield(reg_assess_data[reg], Symbol("valid_$(rtype)"))

    # Filter look up table down to locations that meet criteria
    sel = true
    for (crit, lb, ub) in zip(criteria_names, lbs, ubs)
        lb = parse(Float32, lb)
        ub = parse(Float32, ub)

        sel = sel .& (lb .<= assess_locs[:, crit] .<= ub)
    end

    assess_locs = copy(assess_locs[sel, :])

    target_crs = crs(assess_locs.geometry[1])
    if isnothing(target_crs)
        target_crs = EPSG(4326)
    end

    # Estimate maximum number of pixels in the target area
    res = degrees_to_meters(step(dims(reg_assess_data[reg].stack)[1]), assess_locs.lats[1])
    max_count = floor(x_dist * y_dist / (res * res))

    assess_locs[:, :suitability_score] .= 0.0

    @debug "Assessing target area for $(nrow(assess_locs)) locations"
    search_box = initial_search_box(
        (assess_locs[1, :lons], assess_locs[1, :lats]),  # center point
        x_dist,                # x distance in meters
        y_dist,                # y distance in meters
        target_crs,
        0.0
    )

    # Create KD-tree to identify `n` nearest pixels
    # Using the lookup method for very large numbers of locations is prohibitively slow
    # hence why we use a KD-tree
    lon_lats = Matrix(assess_locs[:, [:lons, :lats]])'
    kdtree = KDTree(lon_lats; leafsize=25)
    @debug "$(now()) : Assessing suitability for $(nrow(assess_locs)) locations"
    for (i, coords) in enumerate(eachcol(lon_lats))
        # Retrieve the closest valid pixels
        idx, _ = knn(kdtree, coords, ceil(Int64, max_count))

        sb = move_geom(search_box, Tuple(coords))
        assess_locs[i, :suitability_score] = floor(
            Int64,
            (count(GO.contains.(Ref(sb), assess_locs[idx, :geometry])) / max_count) * 100
        )
    end
    @debug "$(now()) : Finished suitability assessment"

    return assess_locs
end

"""
    assess_region(reg_assess_data, reg, qp, rtype)

Perform raster suitability assessment based on user-defined criteria.

# Arguments
- `reg_assess_data` : Dictionary containing the regional data paths, reef outlines and \
                      full region names.
- `reg` : Name of the region being assessed (format `Cairns-Cooktown` rather than \
          `Cairns/Cooktown Management Area`).
- `qp` : Dict containing bounds for each variable being filtered.
- `rtype` : Type of zone to assess (flats or slopes).

# Returns
GeoTiff file of surrounding hectare suitability (1-100%) based on the criteria bounds input
by a user.
"""
function assess_region(reg_assess_data, reg, qp, rtype)::Raster
    # Make mask of suitable locations
    @debug "$(now()) : Creating mask for region"
    mask_data = mask_region(reg_assess_data, reg, qp, rtype)

    # Assess remaining pixels for their suitability
    @debug "$(now()) : Calculating proportional suitability score"
    suitability_scores = proportion_suitable(mask_data.data)

    @debug "$(now()) : Rebuilding raster and returning results"
    return rebuild(mask_data, suitability_scores)
end

"""
    assess_region(config, qp::Dict, reg::String, rtype::String, reg_assess_data::OrderedDict)

Convenience method wrapping around the analysis conducted by `assess_region()`.
Checks for previous assessment of indicated region and returns filename of cache if found.

"""
function assess_region(
    config, qp::Dict, reg::String, rtype::String, reg_assess_data::OrderedDict
)
    assessed_fn = cache_filename(
        extract_criteria(qp, suitability_criteria()), config, "$(reg)_suitable", "tiff"
    )
    if isfile(assessed_fn)
        return assessed_fn
    end

    @debug "$(now()) : Assessing region $(reg)"
    assessed = assess_region(reg_assess_data, reg, qp, rtype)

    @debug "$(now()) : Writing to $(assessed_fn)"
    _write_tiff(assessed_fn, assessed)

    return assessed_fn
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

    @debug "$(now()) : Identifying search pixels"
    target_locs = identify_search_pixels(assess_locs, x -> x .> suitability_threshold)

    if size(target_locs, 1) == 0
        # No viable set of locations, return empty dataframe
        return DataFrame(;
            score=[],
            orientation=[],
            qc_flag=[],
            geometry=[]
        )
    end

    # Need reef outlines to indicate direction of the reef edge
    gdf = REGIONAL_DATA["reef_outlines"]

    x_dist = parse(Int64, site_criteria["xdist"])
    y_dist = parse(Int64, site_criteria["ydist"])
    @debug "$(now()) : Assessing site polygons for $(size(target_locs, 1)) locations."
    initial_polygons = identify_edge_aligned_sites(
        crit_pixels,
        target_locs,
        res,
        gdf,
        x_dist,
        y_dist,
        target_crs,
        reg
    )

    return initial_polygons
end
