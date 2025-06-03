"""
Helpers for job handlers which interrupt main workflow. 

For example, converting between job system interfaces and assessment interfaces.
"""

"""
Build regional assessment parameters from user input and regional data.

Creates a complete parameter set for regional assessment by merging user-specified
criteria bounds with regional defaults. Validates that the specified region exists.

# Arguments
- `input::RegionalAssessmentInput` : User input containing assessment parameters
- `regional_data::RegionalData` : Complete regional data for validation and defaults

# Returns
`RegionalAssessmentParameters` struct ready for assessment execution.

# Throws
- `ErrorException` : If specified region is not found in regional data
"""
function build_regional_assessment_parameters(
    input::RegionalAssessmentInput,
    regional_data::RegionalData
)::RegionalAssessmentParameters
    @info "Building regional assessment parameters" region = input.region

    # Validate region exists
    if !haskey(regional_data.regions, input.region)
        available_regions = collect(keys(regional_data.regions))
        @error "Region not found in regional data" region = input.region available_regions
        throw(
            ErrorException(
                "Regional data did not have data for region $(input.region). Available regions: $(join(available_regions, ", "))"
            )
        )
    end

    region_data = regional_data.regions[input.region]
    regional_criteria::BoundedCriteriaDict = Dict()
    regional_bounds::BoundedCriteriaDict = region_data.criteria

    for (criteria_id, possible_symbols) in PARAM_MAP
        bounds = get(regional_bounds, criteria_id, nothing)
        user_min =
            isnothing(possible_symbols) ? nothing :
            getproperty(input, first(possible_symbols))
        user_max =
            isnothing(possible_symbols) ? nothing :
            getproperty(input, last(possible_symbols))

        merged = merge_bounds(
            user_min,
            user_max,
            bounds
        )
        if !isnothing(merged)
            regional_criteria[criteria_id] = BoundedCriteria(;
                metadata=ASSESSMENT_CRITERIA[criteria_id],
                bounds=merged
            )
        end
    end

    return RegionalAssessmentParameters(;
        region=input.region,
        regional_criteria,
        region_data
    )
end

"""
Build suitability assessment parameters from user input and regional data.

Creates a complete parameter set for suitability assessment by merging user-specified
criteria bounds with regional defaults. Includes spatial dimensions for polygon analysis.

# Arguments
- `input::SuitabilityAssessmentInput` : User input containing assessment parameters and spatial dimensions
- `regional_data::RegionalData` : Complete regional data for validation and defaults

# Returns
`SuitabilityAssessmentParameters` struct ready for assessment execution.

# Throws
- `ErrorException` : If specified region is not found in regional data
"""
function build_suitability_assessment_parameters(
    input::SuitabilityAssessmentInput,
    regional_data::RegionalData
)::SuitabilityAssessmentParameters
    @info "Building suitability assessment parameters" region = input.region x_dist =
        input.x_dist y_dist = input.y_dist

    @debug "Building regional parameters first"
    regional_input = regional_job_from_suitability_job(input)
    regional_parameters = build_regional_assessment_parameters(
        regional_input,
        regional_data
    )
    # Extract threshold with default fallback
    threshold =
        !isnothing(input.threshold) ? input.threshold : DEFAULT_SUITABILITY_THRESHOLD
    @debug "Extending regional parameters with suitability inputs x_dist and ydist" x =
        input.x_dist y = input.y_dist
    return SuitabilityAssessmentParameters(;
        region=regional_parameters.region,
        regional_criteria=regional_parameters.regional_criteria,
        region_data=regional_parameters.region_data,
        suitability_threshold=Int64(threshold),
        x_dist=input.x_dist,
        y_dist=input.y_dist
    )
end

"""
Converts parameters from a suitability job into a regional job
"""
function regional_job_from_suitability_job(
    suitability_job::SuitabilityAssessmentInput
)::RegionalAssessmentInput
    return RegionalAssessmentInput(
        suitability_job.region,
        suitability_job.reef_type,
        suitability_job.depth_min,
        suitability_job.depth_max,
        suitability_job.slope_min,
        suitability_job.slope_max,
        suitability_job.rugosity_min,
        suitability_job.rugosity_max,
        suitability_job.waves_period_min,
        suitability_job.waves_period_max,
        suitability_job.waves_height_min,
        suitability_job.waves_height_max
    )
end
