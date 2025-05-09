using Base.Threads

using Dates
using StructTypes

import Rasters: Between
using Oxygen: json, Request

include("query_parser.jl")
include("tiles.jl")
include("site_identification.jl")

const REEF_TYPE = [:slopes, :flats]

# HTTP response headers for COG files
const COG_HEADERS = [
    "Cache-Control" => "max-age=86400, no-transform"
]

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

struct RegionalCriteria{T}
    stack::RasterStack
    valid_slopes::T
    valid_flats::T
end

function valid_slope_lon_inds(reg::RegionalCriteria)
    return reg.valid_slopes.lon_idx
end
function valid_slope_lat_inds(reg::RegionalCriteria)
    return reg.valid_slopes.lat_idx
end
function valid_flat_lon_inds(reg::RegionalCriteria)
    return reg.valid_flats.lon_idx
end
function valid_flat_lat_inds(reg::RegionalCriteria)
    return reg.valid_flats.lat_idx
end

function Base.show(io::IO, ::MIME"text/plain", z::RegionalCriteria)
    # TODO: Include the extent
    println("""
    Criteria: $(names(z.stack))
    Number of valid slope locations: $(nrow(z.valid_slopes))
    Number of valid flat locations: $(nrow(z.valid_flats))
    """)
    return nothing
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
end

# Define struct type definition to auto-serialize/deserialize to JSON
StructTypes.StructType(::Type{CriteriaBounds}) = StructTypes.Struct()

"""
    criteria_middleware(handle)

Creates middleware that parses a criteria query before reaching an endpoint

# Example
`https://somewhere:8000/suitability/assess/region-name/reeftype?criteria_names=Depth,Slope&lb=-9.0,0.0&ub=-2.0,40`
"""
function criteria_middleware(handle)
    function (req)
        fd = queryparams(req)

        criteria_names = string.(split(fd["criteria_names"], ","))
        lbs = string.(split(fd["lb"], ","))
        ubs = string.(split(fd["ub"], ","))

        return handle(CriteriaBounds.(criteria_names, lbs, ubs))
    end
end

function setup_criteria_routes(config, auth)
    reg_assess_data = setup_regional_data(config)

    @get auth("/criteria/ranges") function (req::Request)
        # Transform cached criteria ranges to a dictionary for return as json.
        criteria_ranges = reg_assess_data["criteria_ranges"]
        criteria_names = names(criteria_ranges)

        @debug "Transforming criteria dataframe to JSON"
        ret = OrderedDict()
        for cn in criteria_names
            ret[cn] = OrderedDict(
                :min_val => criteria_ranges[1, cn],
                :max_val => criteria_ranges[2, cn]
            )
        end

        return json(ret)
    end
end
