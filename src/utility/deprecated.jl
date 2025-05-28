"""
==========
DEPRECATED 
==========
"""

function criteria_data_map()
    # TODO: Load from config?
    return OrderedDict(
        :Depth => "_bathy",
        :Benthic => "_benthic",
        :Geomorphic => "_geomorphic",
        :Slope => "_slope",
        :Turbidity => "_turbid",
        :WavesHs => "_waves_Hs",
        :WavesTp => "_waves_Tp",
        :Rugosity => "_rugosity",
        :ValidSlopes => "_valid_slopes",
        :ValidFlats => "_valid_flats"
    )
end

function search_criteria()::Vector{String}
    return string.(keys(criteria_data_map()))
end

function site_criteria()::Vector{String}
    return ["SuitabilityThreshold", "xdist", "ydist"]
end

function suitability_criteria()::Vector{String}
    return vcat(search_criteria(), ["SuitabilityThreshold"])
end

function criteria_data_map()
    # TODO: Load from config?
    return OrderedDict(
        :Depth => "_bathy",
        :Benthic => "_benthic",
        :Geomorphic => "_geomorphic",
        :Slope => "_slope",
        :Turbidity => "_turbid",
        :WavesHs => "_waves_Hs",
        :WavesTp => "_waves_Tp",
        :Rugosity => "_rugosity",
        :ValidSlopes => "_valid_slopes",
        :ValidFlats => "_valid_flats"

        # Unused datasets
        # :PortDistSlopes => "_PortDistSlopes",
        # :PortDistFlats => "_PortDistFlats"
    )
end

function search_criteria()::Vector{String}
    return string.(keys(criteria_data_map()))
end

function site_criteria()::Vector{String}
    return ["SuitabilityThreshold", "xdist", "ydist"]
end

function suitability_criteria()::Vector{String}
    return vcat(search_criteria(), ["SuitabilityThreshold"])
end

function extract_criteria(qp::T, criteria::Vector{String})::T where {T<:Dict{String,String}}
    return filter(
        k -> string(k.first) ∈ criteria, qp
    )
end

struct OldRegionalCriteria{T}
    stack::RasterStack
    valid_slopes::T
    valid_flats::T
end

function valid_slope_lon_inds(reg::OldRegionalCriteria)
    return reg.valid_slopes.lon_idx
end
function valid_slope_lat_inds(reg::OldRegionalCriteria)
    return reg.valid_slopes.lat_idx
end
function valid_flat_lon_inds(reg::OldRegionalCriteria)
    return reg.valid_flats.lon_idx
end
function valid_flat_lat_inds(reg::OldRegionalCriteria)
    return reg.valid_flats.lat_idx
end

struct CriteriaBounds{F<:Function}
    name::Symbol
    lower_bound::Float32
    upper_bound::Float32
    rule::F

    function CriteriaBounds(name::String, lb::S, ub::S)::CriteriaBounds where {S<:String}
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
