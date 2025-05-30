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
    assess_region(reg_assess_data, reg, qp, rtype)

Perform raster suitability assessment based on user-defined criteria.

# Arguments
- params :: RegionalAssessmentParameters

# Returns
GeoTiff file of surrounding hectare suitability (1-100%) based on the criteria bounds input
by a user.
"""
function assess_region(params::RegionalAssessmentParameters)::Raster
    # Make mask of suitable locations
    @debug "$(now()) : Creating mask for region"

    # Builds out a set of criteria filters using the regional criteria.
    # NOTE this will only filter over available criteria
    filters = build_criteria_bounds_from_regional_criteria(params.regional_criteria)

    # map our regional criteria 
    @debug "Applying criteria thresholds to generate mask layer"
    mask_data = filter_raster_by_criteria(
        # This is the raster stack
        params.region_data.raster_stack,
        # The slope table dataframe
        params.region_data.slope_table,
        # The list of criteria bounds
        filters
    )

    # Assess remaining pixels for their suitability
    @debug "$(now()) : Calculating proportional suitability score"
    suitability_scores = proportion_suitable(mask_data.data)

    @debug "$(now()) : Rebuilding raster and returning results"
    return rebuild(mask_data, suitability_scores)
end

function assess_sites(
    params::SuitabilityAssessmentParameters,
    regional_raster::Raster
)
    target_crs = convert(EPSG, crs(regional_raster))
    suitability_threshold = params.suitability_threshold
    region = params.region

    @debug "$(now()) : Identifying search pixels for $(region)"
    target_locs = lookup_df_from_raster(regional_raster, suitability_threshold)

    if size(target_locs, 1) == 0
        # No viable set of locations, return empty dataframe
        return DataFrame(;
            score=[],
            orientation=[],
            qc_flag=[],
            geometry=[]
        )
    end

    # Otherwise, create the file
    @debug "$(now()) : Assessing criteria table for $(region)"
    # Get criteria bounds list from criteria
    filters = build_criteria_bounds_from_regional_criteria(params.regional_criteria)

    crit_pixels::DataFrame = filter_lookup_table_by_criteria(
        # Slope table
        params.region_data.slope_table,
        filters
    )

    res = abs(step(dims(regional_raster, X)))
    @debug "$(now()) : Assessing $(size(target_locs, 1)) candidate locations in $(region)."
    @debug "Finding optimal site alignment"
    initial_polygons = find_optimal_site_alignment(
        crit_pixels,
        target_locs,
        res,
        params.x_dist,
        params.y_dist,
        target_crs
    )

    return initial_polygons
end
