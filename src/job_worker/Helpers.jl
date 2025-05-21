"""
Helper methods for criteria parsing etc.
"""

# Default threshold when not provided in inputs
const DEFAULT_SUITABILITY_THRESHOLD = 50

"""
Min/max storage for environmental criteria
"""
mutable struct Range
    min::Float32
    max::Float32
    label::String
end

"""
Serialises the value of the range to min:max format
"""
function serialise_range(range::Range)::String
    return "$(range.min):$(range.max)"
end

"""
Builds a dictionary key-value pair from a Range
"""
function range_entry_to_kvp(range::Range)::Tuple{String,String}
    return (range.label, serialise_range(range))
end

"""
Typed/structured ranges for all notable criteria
"""
mutable struct RelevantRanges
    depth::Range
    slope::Range
    turbidity::Range
    waves_height::Range
    waves_period::Range
    rugosity::Range

    function RelevantRanges(;
        depth::Range,
        slope::Range,
        turbidity::Range,
        waves_height::Range,
        waves_period::Range,
        rugosity::Range
    )
        return new(depth, slope, turbidity, waves_height, waves_period, rugosity)
    end
end

"""
Dumps the relevant ranges into a dictionary following the appropriate style for
other methods
"""
function relevant_ranges_to_dict(ranges::RelevantRanges)::Dict{String,String}
    return Dict{String,String}(
        range_entry_to_kvp.([
            ranges.depth,
            ranges.slope,
            ranges.turbidity,
            ranges.waves_height,
            ranges.waves_period,
            ranges.rugosity
        ])
    )
end

"""
Converts DataFrame criteria ranges to a structured RelevantRanges object
"""
function structured_ranges_from_criteria_ranges(criteria_ranges::DataFrame)::RelevantRanges
    # Extract min/max values from the criteria_ranges DataFrame
    depth_range = Range(
        criteria_ranges[1, "Depth"],
        criteria_ranges[2, "Depth"],
        "Depth"
    )

    slope_range = Range(
        criteria_ranges[1, "Slope"],
        criteria_ranges[2, "Slope"],
        "Slope"
    )

    turbidity_range = Range(
        criteria_ranges[1, "Turbidity"],
        criteria_ranges[2, "Turbidity"],
        "Turbidity"
    )

    waves_height_range = Range(
        criteria_ranges[1, "WavesHs"],
        criteria_ranges[2, "WavesHs"],
        "WavesHs"
    )

    waves_period_range = Range(
        criteria_ranges[1, "WavesTp"],
        criteria_ranges[2, "WavesTp"],
        "WavesTp"
    )

    rugosity_range = Range(
        criteria_ranges[1, "Rugosity"],
        criteria_ranges[2, "Rugosity"],
        "Rugosity"
    )

    return RelevantRanges(
        depth=depth_range,
        slope=slope_range,
        turbidity=turbidity_range,
        waves_height=waves_height_range,
        waves_period=waves_period_range,
        rugosity=rugosity_range
    )
end

"""
Applies optional min/max overrides to a Range object
"""
function apply_optional_overrides!(range::Range, min_value::OptionalValue{Float64}, max_value::OptionalValue{Float64})
    if !isnothing(min_value)
        range.min = min_value
    end
    if !isnothing(max_value)
        range.max = max_value
    end
end

"""
Applies all criteria overrides from input to a RelevantRanges object
"""
function apply_criteria_overrides!(ranges::RelevantRanges, criteria::Union{RegionalAssessmentInput, SuitabilityAssessmentInput})
    # Apply overrides for each range
    apply_optional_overrides!(ranges.depth, criteria.depth_min, criteria.depth_max)
    apply_optional_overrides!(ranges.slope, criteria.slope_min, criteria.slope_max)
    apply_optional_overrides!(ranges.rugosity, criteria.rugosity_min, criteria.rugosity_max)
    apply_optional_overrides!(ranges.waves_period, criteria.waves_period_min, criteria.waves_period_max)
    apply_optional_overrides!(ranges.waves_height, criteria.waves_height_min, criteria.waves_height_max)
end

"""
Builds a parameters dictionary from RegionalAssessmentInput, applying overrides as needed
"""
function build_params_dictionary_from_regional_input(
    # The regional criteria job input
    criteria::RegionalAssessmentInput,
    # Criteria ranges as [[min,max], name] i.e. [2, name] = max of name
    criteria_ranges::DataFrame
)::Dict{String,String}
    # Get the structured ranges
    default_ranges = structured_ranges_from_criteria_ranges(criteria_ranges)
    
    # Apply all overrides
    apply_criteria_overrides!(default_ranges, criteria)
    
    # Base dictionary of ranges
    ranges_dict = relevant_ranges_to_dict(default_ranges)
    
    # Add in suitability threshold
    threshold_value = isnothing(criteria.threshold) ? DEFAULT_SUITABILITY_THRESHOLD : criteria.threshold
    ranges_dict["SuitabilityThreshold"] = "$(threshold_value)"
    
    return ranges_dict
end

"""
Builds a parameters dictionary from SuitabilityAssessmentInput, applying overrides and adding suitability-specific parameters
"""
function build_params_dictionary_from_suitability_input(
    # The suitability criteria job input
    criteria::SuitabilityAssessmentInput,
    # Criteria ranges as [[min,max], name] i.e. [2, name] = max of name
    criteria_ranges::DataFrame
)::Dict{String,String}
    # Get the structured ranges
    default_ranges = structured_ranges_from_criteria_ranges(criteria_ranges)
    
    # Apply all overrides
    apply_criteria_overrides!(default_ranges, criteria)
    
    # Base dictionary of ranges
    ranges_dict = relevant_ranges_to_dict(default_ranges)
    
    # Add in suitability threshold
    threshold_value = isnothing(criteria.threshold) ? DEFAULT_SUITABILITY_THRESHOLD : criteria.threshold
    ranges_dict["SuitabilityThreshold"] = "$(threshold_value)"
    
    # Suitability specific entries
    ranges_dict["xdist"] = "$(criteria.x_dist)"
    ranges_dict["ydist"] = "$(criteria.y_dist)"
    
    return ranges_dict
end
