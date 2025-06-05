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
    assess_region(
        params::RegionalAssessmentParameters
    )::Raster

Perform raster suitability assessment based on user-defined criteria.

# Arguments
- params :: RegionalAssessmentParameters

# Returns
GeoTiff file of surrounding hectare suitability (1-100%) based on the criteria bounds input
by a user.
"""
function assess_region(
    params::RegionalAssessmentParameters
)::Raster
    # Make mask of suitable locations
    @debug "$(now()) : Creating mask for region"

    # Estimate size of mask (in MBs)
    lookup_tbl = params.region_data.slope_table
    region_dims = maximum(lookup_tbl.lon_idx), maximum(lookup_tbl.lat_idx)
    mask_size_MB = prod(region_dims) / 1024^2

    # Arbitrary threshold - use sparse matrices if likely to exceed memory
    if mask_size_MB < 700
        @debug "$(now()) : Creating mask as a regular matrix (est. size: $(mask_size_MB))"
        indicator = zeros(
            Int8, maximum(lookup_tbl.lon_idx), maximum(lookup_tbl.lat_idx)
        )
    else
        @debug "$(now()) : Creating mask as a sparse matrix (est. size: $(mask_size_MB))"
        indicator = ExtendableSparseMatrix{Int8,Int64}(
            maximum(lookup_tbl.lon_idx), maximum(lookup_tbl.lat_idx)
        )
    end

    # Builds out a set of criteria filters using the regional criteria.
    # NOTE this will only filter over available criteria
    filters = build_criteria_bounds_from_regional_criteria(params.regional_criteria)

    @debug "$(now()) : Applying criteria thresholds to generate mask layer"
    matching_idx = filter_lookup_table_by_criteria(lookup_tbl, filters)
    for r in eachrow(lookup_tbl[matching_idx, :])
        indicator[r.lon_idx, r.lat_idx] = true
    end

    return Raster(params.region_data.valid_extent; data=indicator)
end

function assess_sites(params::SuitabilityAssessmentParameters)
    target_crs = convert(EPSG, crs(regional_raster))
    suitability_threshold = params.suitability_threshold
    region = params.region

    @debug "$(now()) : Assessing criteria table for $(region)"
    # Get criteria bounds list from criteria
    filters = build_criteria_bounds_from_regional_criteria(params.regional_criteria)

    slope_table = params.region_data.slope_table
    matching_idx::BitVector = filter_lookup_table_by_criteria(
        slope_table,
        filters
    )

    res = abs(step(dims(regional_raster, X)))
    @debug "$(now()) : Assessing $(count(matching_idx)) candidate locations in $(region)."
    @debug "Finding optimal site alignment"
    initial_polygons = find_optimal_site_alignment(
        slope_table[matching_idx, :],
        res,
        params.x_dist,
        params.y_dist;
        target_crs=target_crs,
        degree_step=20.0,
        suit_threshold=suitability_threshold
    )

    return initial_polygons
end
