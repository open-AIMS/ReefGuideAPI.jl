using Base.Threads

using Dates
using StructTypes

import Rasters: Between
using CairoMakie, ImageIO, Images


const REEF_TYPE = [:slopes, :flats]

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
https:://somewhere:8000/suitability/assess/region-name/reeftype?criteria_names=Depth,Slope&lb=-9.0,0.0&ub=-2.0,40
"""
function criteria_middleware(handle)
    function(req)
        fd = queryparams(req)

        criteria_names = string.(split(fd["criteria_names"], ","))
        lbs = string.(split(fd["lb"], ","))
        ubs = string.(split(fd["ub"], ","))

        return handle(CriteriaBounds.(criteria_names, lbs, ubs))
    end
end

# Not sure why, but this ain't working
# reeftype_router = router("/suitability", middleware=[criteria_middleware], tags=["suitability"])

"""
    _tile_to_lon_lat(z, x, y)

Obtain lon/lat bounds of a requested tile.

Note: Zoom levels are capped to 14
"""
function _tile_to_lon_lat(z, x, y)
    z = min(14, z)  # Cap zoom level to something reasonable

    # Calculate the boundaries of the tile
    n = 2.0^z
    lat_min = atan(sinh(π * (1 - 2 * (y + 1) / n))) * 180.0 / π
    lat_max = atan(sinh(π * (1 - 2 * y / n))) * 180.0 / π
    lon_min = x / n * 360.0 - 180.0
    lon_max = (x + 1) / n * 360.0 - 180.0

    # West, East, North, South
    return lon_min, lon_max, lat_min, lat_max
end

"""
    _lon_lat_to_tile(lon, lat, zoom)

Identify the corresponding tile coordinates for a given lon/lat.
"""
function _lon_lat_to_tile(lon, lat, zoom)
    n = 2.0^zoom
    x = floor(Int64, (lon + 180.0) / 360.0 * n)

    lat_rad = lat * π / 180.0
    y = floor(Int64, (1.0 - log(tan(lat_rad) + 1.0 / cos(lat_rad)) / π) / 2.0 * n)

    return x, y
end

"""
    fast_resample(rst::Raster, new_size::Tuple{Int, Int})::Matrix

Quickly resample a Raster to a new size using nearest neighbor interpolation.

Applies a simple nearest neighbor algorithm, prioritising performance over accuracy.


# Arguments
- `rst::Raster`: The input raster to be resampled.
- `new_size::Tuple{Int, Int}`: The desired dimensions of the output raster as (width, height).

# Returns
Matrix with the resampled data.

# Notes
- The approach prioritizes speed over accuracy and is particularly effective for downsampling large rasters.
- It uses nearest neighbor interpolation, which may result in aliasing artifacts, especially when upsampling.
- The function is multi-threaded. Ensure Julia is started with multiple threads for optimal performance.

# Example
```julia
large_raster = Raster(rand(UInt8, 14756, 14838); dims=(X(1:1:14756), Y(1:1:14838)))
small_matrix = fast_resample(large_raster, (256, 256))
```
"""
function fast_resample(rst::Raster, new_size::Tuple{Int, Int})::Matrix
    old_size = size(rst)
    x_ratio = old_size[1] / new_size[1]
    y_ratio = old_size[2] / new_size[2]

    resampled = Array{eltype(rst)}(undef, new_size)

    Threads.@threads for y in 1:new_size[2]
        for x in 1:new_size[1]
            src_x = round(Int, x * x_ratio)
            src_y = round(Int, y * y_ratio)
            resampled[x, y] = rst[src_x, src_y]
        end
    end

    return resampled
end

function setup_region_routes(config)
    reg_assess_data = setup_regional_data(config)

    @get "/assess/{reg}/{rtype}" function (req::Request, reg::String, rtype::String)
        qp = queryparams(req)
        file_id = string(hash(qp))
        mask_temp_path = _cache_location(config)
        mask_path = joinpath(mask_temp_path, file_id*".tiff")

        if isfile(mask_path)
            return file(mask_path)
        end

        # Otherwise, create the file
        @debug "$(now()) : Assessing criteria"
        # Filtering time: 0.6 - 7.0 seconds
        criteria_names = string.(split(qp["criteria_names"], ","))
        lbs = string.(split(qp["lb"], ","))
        ubs = string.(split(qp["ub"], ","))

        if !contains(reg, "Townsville")
            # Remove rugosity layer from consideration as it doesn't exist for regions
            # outside of Townsville.
            pos = findfirst(lowercase.(criteria_names) .== "rugosity")
            criteria_names = [cname for (i, cname) in enumerate(criteria_names) if i != pos]
            lbs = [lb for (i, lb) in enumerate(lbs) if i != pos]
            ubs = [ub for (i, ub) in enumerate(ubs) if i != pos]
        end

        assess = reg_assess_data[reg]
        mask_data = make_threshold_mask(
            assess,
            Symbol(rtype),
            CriteriaBounds.(criteria_names, lbs, ubs)
        )

        @debug "$(now()) : Running on thread $(threadid())"
        @debug "Writing to $(mask_path)"
        # Writing time: ~10-25 seconds
        rst_layer = assess.stack[names(assess.stack)[1]]
        Rasters.write(
            mask_path,
            Raster(rst_layer; data=UInt8.(mask_data), missingval=0);
            ext=".tiff",
            source="gdal",
            driver="COG",
            options=Dict{String,String}(
                "COMPRESS"=>"DEFLATE",
                "SPARSE_OK"=>"TRUE",
                "OVERVIEW_COUNT"=>"5",
                "BLOCKSIZE"=>"256",
                "NUM_THREADS"=>n_gdal_threads(config)
            ),
            force=true
        )

        return file(mask_path)
    end

    # Form for testing/dev
    # https:://somewhere:8000/suitability/assess/region-name/reeftype?criteria_names=Depth,Slope&lb=-9.0,0.0&ub=-2.0,40
    @get "/" function()
        html("""
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
    @post "/form" function(req)
        data = formdata(req)
        return data
    end
end


# function setup_criteria_routes()
#     crit_list = criteria_data_map()

#     for (c, data_pattern) in crit_list
#         @get "{region}/$c" function (req, region::S, lb::T, ub::T) where {S,T}
#             return Image(make_threshold_mask(region::String, c::Symbol, crit_map::Dict))

#             loaded_data = retrieve_data(PREPPED_DATA_DIR, data_pattern)
#             return lb .<= loaded_data .<= ub
#         end
#     end
# end
