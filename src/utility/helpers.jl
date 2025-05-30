
"""
Builds a hash by combining strings and hashing result
"""
function build_hash_from_components(components::Vector{String})::String
    return string(hash(join(components, "|")))
end

"""
Returns a hash component for bounded criteria 

Carefully orders hash predictably and only includes present criteria
"""
function get_hash_components_from_regional_criteria(
    criteria::BoundedCriteriaDict
)::Vector{String}
    @debug "Hashing criteria..." criteria
    components::Vector{String} = []
    for id in keys(ASSESSMENT_CRITERIA)
        components = vcat(components,
            haskey(criteria, id) ?
            [
                id, string(criteria[id].bounds.min), string(criteria[id].bounds.max)
            ] : []
        )
    end
    return components
end

"""
Convert BoundedCriteriaDict to a vector of CriteriaBounds for assessment processing.
Transforms the BoundedCriteriaDict into CriteriaBounds objects that include
evaluation functions. Only includes criteria that are available in the dictionary.

# Arguments 
- `bounded_criteria_dict::BoundedCriteriaDict` : Dictionary of bounded criteria to convert

# Returns
Vector of `CriteriaBounds` objects for available criteria.
"""
function build_criteria_bounds_from_regional_criteria(
    bounded_criteria_dict::BoundedCriteriaDict
)::Vector{CriteriaBounds}
    @debug "Converting BoundedCriteriaDict to CriteriaBounds vector"
    criteria_bounds = CriteriaBounds[]

    for (criteria_id, bounded_criteria) in bounded_criteria_dict
        bounds = CriteriaBounds(
            # Field to get in the data
            bounded_criteria.metadata.id,
            # Min/max bounds  
            bounded_criteria.bounds.min,
            bounded_criteria.bounds.max
        )
        push!(criteria_bounds, bounds)
        @debug "Added criteria bounds" criteria_id = criteria_id min_val =
            bounded_criteria.bounds.min max_val = bounded_criteria.bounds.max
    end

    @debug "Built CriteriaBounds vector" total_criteria = length(criteria_bounds) criteria_ids = [
        cb.name for cb in criteria_bounds
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
