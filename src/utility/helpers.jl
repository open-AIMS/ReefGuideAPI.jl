"""
Builds a hash by combining strings and hashing result
"""
function build_hash_from_components(components::Vector{String})::String
    return string(hash(join(components, "|")))
end

"""
Combines present regional criteria including bounds into hash components
"""
function get_hash_components_from_regional_criteria(
    criteria::RegionalCriteria
)::Vector{String}
    hash_components::Vector{String} = []
    for field in REGIONAL_CRITERIA_SYMBOLS
        criteria_entry::OptionalValue{RegionalCriteriaEntry} = getfield(
            criteria, field
        )
        if !isnothing(criteria_entry)
            push!(
                hash_components,
                "$(field)_$(criteria_entry.bounds.min)_$(criteria_entry.bounds.max)"
            )
        else
            push!(hash_components, "$(field)_null")
        end
    end
    return hash_components
end

"""
Convert RegionalCriteria to a vector of CriteriaBounds for assessment processing.

Transforms the RegionalCriteria struct into CriteriaBounds objects that include
evaluation functions. Only includes criteria that are available (non-nothing).

# Arguments
- `regional_criteria::RegionalCriteria` : Regional criteria with bounds to convert

# Returns
Vector of `CriteriaBounds` objects for available criteria.
"""
function build_criteria_bounds_from_regional_criteria(
    regional_criteria::RegionalCriteria
)::Vector{CriteriaBounds}
    @debug "Converting RegionalCriteria to CriteriaBounds vector"

    criteria_bounds = CriteriaBounds[]

    for field_symbol in REGIONAL_CRITERIA_SYMBOLS
        criteria_entry = getfield(regional_criteria, field_symbol)

        if !isnothing(criteria_entry)
            bounds = CriteriaBounds(
                # Field to get in the data
                criteria_entry.metadata.id,
                # Min/max bounds
                criteria_entry.bounds.min,
                criteria_entry.bounds.max
            )
            push!(criteria_bounds, bounds)
        else
            @debug "Skipped criteria - not available" criteria_id = String(field_symbol)
        end
    end

    @debug "Built CriteriaBounds vector" total_criteria = length(criteria_bounds) criteria_ids = [
        String(cb.name) for cb in criteria_bounds
    ]

    return criteria_bounds
end

"""
Convert a min/max tuple to a Bounds struct.

# Arguments
- `min_max::Tuple{Number,Number}` : Tuple containing (minimum, maximum) values

# Returns
`Bounds` struct with converted float values.
"""
function bounds_from_tuple(min_max::Tuple{Number,Number})::Bounds
    return Bounds(; min=min_max[1], max=min_max[2])
end

"""
Generate the filename for slope lookup data for a given region.

# Arguments  
- `region::RegionMetadata` : Region metadata containing ID

# Returns
String filename in format "{region_id}_slope_lookup.parq"
"""
function get_slope_parquet_filename(region::RegionMetadata)::String
    filename = "$(region.id)$(SLOPES_LOOKUP_SUFFIX)"
    @debug "Generated slope parquet filename" region_id = region.id filename
    return filename
end

"""
Create a dictionary mapping criteria IDs to regional criteria entries.

NOTE: Only includes criteria that are available for the region, as specified in the
region metadata and actually instantiated in the RegionalCriteria struct.

Uses the defined set of symbols on the regional criteria struct to iterate through

# Arguments
- `region_data::RegionalDataEntry` : Regional data containing criteria information

# Returns
Dictionary with criteria ID strings as keys and RegionalCriteriaEntry as values.
Only includes criteria that are both listed in region metadata and available as non-nothing.
"""
function build_regional_criteria_dictionary(
    region_data::RegionalDataEntry
)::Dict{String,RegionalCriteriaEntry}
    @debug "Building criteria dictionary for region" region_id = region_data.region_id available_in_metadata =
        region_data.region_metadata.available_criteria

    regional_criteria = region_data.criteria
    criteria_dict = Dict{String,RegionalCriteriaEntry}()

    # Only process criteria that are listed as available in the region metadata
    available_criteria_set = Set(region_data.region_metadata.available_criteria)

    for symbol in REGIONAL_CRITERIA_SYMBOLS
        possible_value::OptionalValue{RegionalCriteriaEntry} = getfield(
            regional_criteria, symbol
        )
        if (
            !isnothing(possible_value) &&
            possible_value.metadata.id ∈ available_criteria_set
        )
            criteria_dict[possible_value.metadata.id] = possible_value
        end
    end

    @debug "Built criteria dictionary" region_id = region_data.region_id available_in_metadata = length(
        available_criteria_set
    ) actually_available = length(criteria_dict) criteria_ids = collect(keys(criteria_dict))

    return criteria_dict
end

"""
Given a dictionary mapping criteria ID -> optional bounds, builds out a
RegionalCriteria object.
"""
function build_regional_criteria_from_criteria_dictionary(
    criteria::Dict{String,OptionalValue{Bounds}}
)
    function check_criteria(metadata::CriteriaMetadata)::OptionalValue{Bounds}
        if haskey(criteria, metadata.id) && !isnothing(criteria[metadata.id])
            return criteria[metadata.id]
        end
        return nothing
    end

    @debug "Creating RegionalCriteria by assessing each entry of criteria dictionary"
    return RegionalCriteria(;
        depth_bounds=check_criteria(ASSESSMENT_CRITERIA.depth),
        slope_bounds=check_criteria(ASSESSMENT_CRITERIA.slope),
        turbidity_bounds=check_criteria(ASSESSMENT_CRITERIA.turbidity),
        waves_height_bounds=check_criteria(ASSESSMENT_CRITERIA.waves_height),
        waves_period_bounds=check_criteria(ASSESSMENT_CRITERIA.waves_period),
        rugosity_bounds=check_criteria(ASSESSMENT_CRITERIA.rugosity)
    )
end

"""
Add longitude and latitude columns to a DataFrame based on geometry centroids.

Modifies the input DataFrame by adding 'lons' and 'lats' columns extracted
from the centroid coordinates of each geometry feature.

# Arguments
- `df::DataFrame` : DataFrame with geometry column containing spatial features
"""
function add_lat_long_columns_to_dataframe(df::DataFrame)::Nothing
    @debug "Adding lat/long columns to DataFrame" num_rows = nrow(df)
    # Extract coordinate tuples from geometry centroids
    coords = GI.coordinates.(df.geometry)
    # Add longitude column (first coordinate)
    df[!, :lons] .= first.(coords)
    # Add latitude column (second coordinate) 
    df[!, :lats] .= last.(coords)
    @debug "Successfully added coordinate columns"
    return nothing
end
