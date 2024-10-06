using Base.Threads

using Dates
using StructTypes

import Rasters: Between
using Oxygen: json

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
        :ValidFlats => "_valid_flats",
        :PortDistSlopes => "_PortDistSlopes",
        :PortDistFlats => "_PortDistFlats"
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
    _write_cog(file_path::String, data::Raster)::Nothing

Write out a COG using common options.

# Arguments
- `file_path` : Path to write data out to
- `data` : Raster data to write out
"""
function _write_cog(file_path::String, data::Raster)::Nothing
    Rasters.write(
        file_path,
        data;
        ext=".tiff",
        source="gdal",
        driver="COG",
        options=Dict{String,String}(
            "COMPRESS" => "DEFLATE",
            "SPARSE_OK" => "TRUE",
            "OVERVIEW_COUNT" => "5",
            "BLOCKSIZE" => tile_size(config),
            "NUM_THREADS" => n_gdal_threads(config)
        ),
        force=true
    )

    return nothing
end

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

# Not sure why, but this ain't working
# reeftype_router = router("/suitability", middleware=[criteria_middleware], tags=["suitability"])

function setup_region_routes(config, auth)
    reg_assess_data = setup_regional_data(config)

    @get auth("/assess/{reg}/{rtype}") function (req::Request, reg::String, rtype::String)
        qp = queryparams(req)
        file_id = string(hash(qp))
        mask_temp_path = _cache_location(config)
        mask_path = joinpath(mask_temp_path, file_id * ".tiff")

        if isfile(mask_path)
            return file(mask_path; headers=COG_HEADERS)
        end

        # Otherwise, create the file
        @debug "$(now()) : Assessing criteria"
        mask_data = mask_region(reg_assess_data, reg, qp, rtype)

        @debug "$(now()) : Running on thread $(threadid())"
        @debug "Writing to $(mask_path)"
        # Writing time: ~10-25 seconds
        Rasters.write(
            mask_path,
            UInt8.(mask_data);
            ext=".tiff",
            source="gdal",
            driver="COG",
            options=Dict{String,String}(
                "COMPRESS" => "DEFLATE",
                "SPARSE_OK" => "TRUE",
                "OVERVIEW_COUNT" => "5",
                "BLOCKSIZE" => "256",
                "NUM_THREADS" => n_gdal_threads(config)
            ),
            force=true
        )

        return file(mask_path; headers=COG_HEADERS)
    end

    @get auth("/suitability/assess/{reg}/{rtype}") function (
        req::Request, reg::String, rtype::String
    )
        # somewhere:8000/suitability/assess/region-name/reeftype?criteria_names=Depth,Slope&lb=-9.0,0.0&ub=-2.0,40
        # 127.0.0.1:8000/suitability/assess/Cairns-Cooktown/slopes?Depth=-4.0:-2.0&Slope=0.0:40.0&Rugosity=0.0:6.0

        qp = queryparams(req)
        return assess_region(reg_assess_data, reg, qp, rtype, config)
    end

    @get auth("/suitability/site-suitability/{reg}/{rtype}") function (
        req::Request, reg::String, rtype::String
    )
        # 127.0.0.1:8000/suitability/site-suitability/Cairns-Cooktown/slopes?Depth=-4.0:-2.0&Slope=0.0:40.0&Rugosity=0.0:6.0&SuitabilityThreshold=95&xdist=450&ydist=50
        qp = queryparams(req)
        criteria_names = string.(keys(criteria_data_map()))
        criteria_qp = filter(k -> k.first ∈ criteria_names, qp)

        assessment_qp = filter(
            k -> string(k.first) ∈ ["SuitabilityThreshold", "xdist", "ydist"], qp
        )

        return site_assess_region(
            reg_assess_data, reg, criteria_qp, assessment_qp, rtype, config
        )
    end

    @get auth("/bounds/{reg}") function (req::Request, reg::String)
        rst_stack = reg_assess_data[reg].stack

        return json(Rasters.bounds(rst_stack))
    end

    # Form for testing/dev
    # https:://somewhere:8000/suitability/assess/region-name/reeftype?criteria_names=Depth,Slope&lb=-9.0,0.0&ub=-2.0,40
    @get "/" function ()
        return html("""
               <form action="/assess/Cairns-Cooktown/slopes" method="post">
                   <label for="criteria_names">Criteria Names:</label><br>
                   <input type="text" id="criteria_names" name="criteria"><br>

                   <label for="lb">Lower Bound:</label><br>
                   <input type="text" id="lb" name="lower_bound"><br><br>

                   <label for="ub">Upper Bound:</label><br>
                   <input type="text" id="ub" name="upper_bound"><br><br>

                   <input type="submit" value="Submit">
               </form>
               """)
    end

    # Parse the form data and return it
    @post auth("/form") function (req)
        data = formdata(req)
        return data
    end
end
