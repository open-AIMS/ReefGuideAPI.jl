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
    _write_cog(file_path::String, data::Raster, config::Dict)::Nothing

Write out a COG using common options.

# Arguments
- `file_path` : Path to write data out to
- `data` : Raster data to write out
"""
function _write_cog(file_path::String, data::Raster, config::Dict)::Nothing
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
            "BLOCKSIZE" => string(first(tile_size(config))),
            "NUM_THREADS" => n_gdal_threads(config)
        ),
        force=true
    )

    return nothing
end

function _write_tiff(file_path::String, data::Raster)::Nothing
    Rasters.write(
        file_path,
        data;
        ext=".tiff",
        source="gdal",
        driver="gtiff",
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
        mask_path = cache_filename(qp, config, "", "tiff")
        if isfile(mask_path)
            return file(mask_path; headers=COG_HEADERS)
        end

        # Otherwise, create the file
        @debug "$(now()) : Assessing criteria"
        mask_data = mask_region(reg_assess_data, reg, qp, rtype)

        @debug "$(now()) : Running on thread $(threadid())"
        @debug "Writing to $(mask_path)"
        # Writing time: ~10-25 seconds
        _write_cog(mask_path, UInt8.(mask_data), config)

        return file(mask_path; headers=COG_HEADERS)
    end

    @get auth("/suitability/assess/{reg}/{rtype}") function (
        req::Request, reg::String, rtype::String
    )
        # somewhere:8000/suitability/assess/region-name/reeftype?criteria_names=Depth,Slope&lb=-9.0,0.0&ub=-2.0,40
        # 127.0.0.1:8000/suitability/assess/Cairns-Cooktown/slopes?Depth=-4.0:-2.0&Slope=0.0:40.0&Rugosity=0.0:6.0
        # 127.0.0.1:8000/suitability/assess/Cairns-Cooktown/slopes?Depth=-4.0:-2.0&Slope=0.0:40.0&Rugosity=0.0:6.0&SuitabilityThreshold=95

        qp = queryparams(req)
        assessed_fn = cache_filename(qp, config, "$(reg)_suitable", "tiff")
        if isfile(assessed_fn)
            return file(assessed_fn; headers=COG_HEADERS)
        end

        @debug "$(now()) : Assessing region $(reg)"
        assessed = assess_region(reg_assess_data, reg, qp, rtype)

        @debug "$(now()) : Writing to $(assessed_fn)"
        _write_tiff(assessed_fn, assessed)

        return file(assessed_fn; headers=COG_HEADERS)
    end

    @get auth("/suitability/site-suitability/{reg}/{rtype}") function (
        req::Request, reg::String, rtype::String
    )
        # 127.0.0.1:8000/suitability/site-suitability/Cairns-Cooktown/slopes?Depth=-4.0:-2.0&Slope=0.0:40.0&Rugosity=0.0:6.0&SuitabilityThreshold=95&xdist=450&ydist=50
        qp = queryparams(req)
        suitable_sites_fn = cache_filename(
            qp, config, "$(reg)_potential_sites", "geojson"
        )
        if isfile(suitable_sites_fn)
            return file(suitable_sites_fn)
        end

        # Assess location suitability if needed
        assessed_fn = cache_filename(qp, config, "$(reg)_suitable", "tiff")
        if isfile(assessed_fn)
            assessed = Raster(assessed_fn)
        else
            assessed = assess_region(reg_assess_data, reg, qp, rtype)
            _write_tiff(assessed_fn, assessed)
        end

        # Extract criteria and assessment
        criteria_names = string.(keys(criteria_data_map()))
        pixel_criteria = filter(k -> k.first ∈ criteria_names, qp)
        site_criteria = filter(
            k -> string(k.first) ∈ ["SuitabilityThreshold", "xdist", "ydist"], qp
        )

        best_sites = filter_sites(
            assess_sites(
                reg_assess_data, reg, rtype, pixel_criteria, site_criteria, assessed
            )
        )

        # Specifically clear from memory to invoke garbage collector
        assessed = nothing

        if size(best_sites, 1) == 0
            return json(nothing)
        end

        output_geojson(suitable_sites_fn, best_sites)
        return file(suitable_sites_fn)
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
