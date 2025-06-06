# =============================================================================
# Assessment Parameters Constants
# =============================================================================

const DEFAULT_SUITABILITY_THRESHOLD::Int64 = 80

# =============================================================================
# Assessment Parameters Data Structures
# =============================================================================

"""
Payload for regional assessment actions - this includes all merged bounds and
regional data.

# Fields
- `region::String` : The region that is being assessed
- `regional_criteria::BoundedCriteriaDict` : The criteria to assess, including user provided bounds
- `region_data::RegionalDataEntry` : The data to consider for this region
"""
struct RegionalAssessmentParameters
    region::String
    regional_criteria::BoundedCriteriaDict
    region_data::RegionalDataEntry

    function RegionalAssessmentParameters(;
        region::String,
        regional_criteria::BoundedCriteriaDict,
        region_data::RegionalDataEntry,
    )
        return new(region, regional_criteria, region_data)
    end
end

"""
Payload for suitability assessment actions - this includes all merged bounds and
regional data plus spatial dimensions.

# Fields
- `region::String` : The region that is being assessed
- `regional_criteria::BoundedCriteriaDict` : The criteria to assess, including user provided bounds
- `region_data::RegionalDataEntry` : The data to consider for this region
- `suitability_threshold::Int64` : The cutoff to consider a site suitable
- `x_dist::Int64` : X dimension of polygon (metres)
- `y_dist::Int64` : Y dimension of polygon (metres)
"""
struct SuitabilityAssessmentParameters
    # Regional criteria
    region::String
    regional_criteria::BoundedCriteriaDict
    region_data::RegionalDataEntry
    suitability_threshold::Int64

    # Additional criteria
    x_dist::Int64
    y_dist::Int64

    function SuitabilityAssessmentParameters(;
        region::String,
        regional_criteria::BoundedCriteriaDict,
        region_data::RegionalDataEntry,
        suitability_threshold::Int64,
        x_dist::Int64,
        y_dist::Int64
    )
        @debug "Created SuitabilityAssessmentParameters" region suitability_threshold x_dist y_dist
        return new(
            region, regional_criteria, region_data, suitability_threshold, x_dist, y_dist
        )
    end
end

# =============================================================================
# Assessment Parameters Helper Functions
# =============================================================================

"""
Merge user-specified bounds with regional defaults.

Creates bounds using user values where provided, falling back to regional 
bounds for unspecified values. Returns nothing if regional criteria is not available.

# Arguments
- `user_min::OptionalValue{Float64}` : User-specified minimum value (optional)
- `user_max::OptionalValue{Float64}` : User-specified maximum value (optional)  
- `regional_criteria::OptionalValue{RegionalCriteriaEntry}` : Regional criteria with default bounds (optional)

# Returns
`Bounds` struct with merged values, or `nothing` if regional criteria unavailable.
"""
function merge_bounds(
    user_min::OptionalValue{Float64},
    user_max::OptionalValue{Float64},
    criteria::OptionalValue{BoundedCriteria}
)::OptionalValue{Bounds}
    if isnothing(criteria)
        return nothing
    end

    bounds = Bounds(;
        min=!isnothing(user_min) ? user_min : criteria.bounds.min,
        max=!isnothing(user_max) ? user_max : criteria.bounds.max
    )

    @debug "Merged bounds" min_val = bounds.min max_val = bounds.max user_specified_min =
        !isnothing(user_min) user_specified_max = !isnothing(user_max)

    return bounds
end

# Parameter mapping: criteria_id => (min_field, max_field) or nothing
const PARAM_MAP::Dict{String,OptionalValue{Tuple{Symbol,Symbol}}} = Dict(
    "Depth" => (:depth_min, :depth_max),
    "Slope" => (:slope_min, :slope_max),
    "Turbidity" => nothing,  # Not user-configurable
    "WavesHs" => (:waves_height_min, :waves_height_max),
    "WavesTp" => (:waves_period_min, :waves_period_max),
    "Rugosity" => (:rugosity_min, :rugosity_max)
)


"""
Generate a deterministic hash string for RegionalAssessmentParameters.

Creates a consistent hash based on assessment parameters that can be used
for cache file naming. Same parameters will always produce the same hash.

# Arguments
- `params::RegionalAssessmentParameters` : Assessment parameters to hash

# Returns
String hash suitable for use in cache file names.
"""
function regional_assessment_params_hash(params::RegionalAssessmentParameters)::String
    @debug "Generating hash for regional assessment parameters" region = params.region

    # Create hash input from key parameters
    hash_components = [
        params.region,
        # spread result list of components from regional criteria
        get_hash_components_from_regional_criteria(params.regional_criteria)...
    ]

    # Create deterministic hash
    hash_string = build_hash_from_components(hash_components)

    @debug "Generated assessment parameters hash" hash = hash_string components_count = length(
        hash_components
    )

    return hash_string
end

"""
Generate a deterministic hash string for SuitabilityAssessmentParameters.

Creates a consistent hash based on assessment parameters that can be used
for cache file naming. Same parameters will always produce the same hash.

# Arguments
- `params::SuitabilityAssessmentParameters` : Assessment parameters to hash

# Returns
String hash suitable for use in cache file names.
"""
function suitability_assessment_params_hash(params::SuitabilityAssessmentParameters)::String
    @debug "Generating hash for suitability assessment parameters" region = params.region threshold =
        params.suitability_threshold x_dist = params.x_dist y_dist = params.y_dist

    # Create hash input from key parameters including spatial dimensions
    hash_components = [
        params.region,
        string(params.suitability_threshold),
        string(params.x_dist),
        string(params.y_dist)
    ]

    # Add criteria bounds to hash (only non-nothing criteria)
    hash_components::Vector{String} = [
        hash_components;
        get_hash_components_from_regional_criteria(params.regional_criteria)
    ]

    # Create deterministic hash
    hash_string = build_hash_from_components(hash_components)

    @debug "Generated suitability parameters hash" hash = hash_string components_count = length(
        hash_components
    )

    return hash_string
end

"""
Build predictable file path for regional assessment results in configured cache
location.

Creates a complete file path for caching regional assessment results using the
configured cache directory and deterministic parameter-based naming.

# Arguments
- `params::RegionalAssessmentParameters` : Regional assessment parameters
- `ext::String` : File extension for the cache file
- `config::Dict` : Configuration dictionary containing cache settings

# Returns
String path to cache file location.
"""
function build_regional_assessment_file_path(
    params::RegionalAssessmentParameters;
    ext::String,
    config::Dict
)::String
    @debug "Building file path for regional assessment cache" region = params.region ext

    cache_path = _cache_location(config)
    param_hash = regional_assessment_params_hash(params)
    filename = "$(param_hash)_$(params.region)_regional_assessment.$(ext)"
    file_path = joinpath(cache_path, filename)

    @debug "Built regional assessment file path" file_path region = params.region hash =
        param_hash

    return file_path
end
