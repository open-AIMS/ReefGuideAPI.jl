"""Methods to filter criteria bounds over rasters and lookup tables"""

"""
CriteriaBounds combine lookup information for a given criteria, bounds, and a
rule (function) which enforces it for a given value
"""
struct CriteriaBounds{F<:Function}
    "The field ID of the criteria"
    name::Symbol
    "min"
    lower_bound::Float32
    "max"
    upper_bound::Float32
    "A function which takes a value and returns if matches the criteria"
    rule::F

    function CriteriaBounds(name::S, lb::S, ub::S)::CriteriaBounds where {S<:String}
        lower_bound::Float32 = parse(Float32, lb)
        upper_bound::Float32 = parse(Float32, ub)
        func = (x) -> lower_bound .<= x .<= upper_bound

        return new{Function}(Symbol(name), lower_bound, upper_bound, func)
    end

    function CriteriaBounds(
        name::String, lb::Float32, ub::Float32
    )::CriteriaBounds
        func = (x) -> lb .<= x .<= ub
        return new{Function}(Symbol(name), lb, ub, func)
    end
end

"""
Apply thresholds for each criteria.

# Arguments
- `criteria_stack` : RasterStack of criteria data for a given region
- `lookup` : Lookup dataframe for the region
- `criteria_bounds` : A vector of CriteriaBounds which contains named criteria
  with min/max ranges and a function to apply.

# Returns
BitMatrix indicating locations within desired thresholds
"""
function filter_raster_by_criteria(
    criteria_stack::RasterStack,
    lookup::DataFrame,
    criteria_bounds::Vector{CriteriaBounds}
)::Raster
    # Result store
    data = falses(size(criteria_stack))

    # Apply criteria
    res_lookup = trues(nrow(lookup))
    for filter::CriteriaBounds in criteria_bounds
        res_lookup .= res_lookup .& filter.rule(lookup[!, filter.name])
    end

    tmp = lookup[res_lookup, [:lon_idx, :lat_idx]]
    data[CartesianIndex.(tmp.lon_idx, tmp.lat_idx)] .= true

    res = Raster(criteria_stack.Depth; data=sparse(data), missingval=0)
    return res
end

"""
Filters the slope table (which contains raster param values too) by building a
bit mask AND'd for all thresholds
"""
function filter_lookup_table_by_criteria(
    slope_table::DataFrame,
    ruleset::Vector{CriteriaBounds}
)::DataFrame
    slope_table.all_crit .= 1

    for threshold in ruleset
        slope_table.all_crit =
            slope_table.all_crit .& threshold.rule(slope_table[!, threshold.name])
    end

    return slope_table[BitVector(slope_table.all_crit), :]
end

"""
    lookup_df_from_raster(raster::Raster, threshold::Union{Int64,Float64})::DataFrame

Build a look up table identifying all pixels in a raster that meet a suitability threshold.

# Arguments
- `raster` : Raster of regional data
- `threshold` : Suitability threshold value (greater or equal than)

# Returns
DataFrame containing indices, lon and lat for each pixel that is intended for further
analysis.
"""
function lookup_df_from_raster(raster::Raster, threshold::Union{Int64,Float64})::DataFrame
    criteria_matches::SparseMatrixCSC{Bool,Int64} = sparse(falses(size(raster)))
    Rasters.read!(raster .>= threshold, criteria_matches)
    indices::Vector{CartesianIndex{2}} = findall(criteria_matches)
    indices_lon::Vector{Float64} = lookup(raster, X)[first.(Tuple.(indices))]
    indices_lat::Vector{Float64} = lookup(raster, Y)[last.(Tuple.(indices))]

    return DataFrame(; indices=indices, lons=indices_lon, lats=indices_lat)
end
