using Base.Threads

using Dates
using StructTypes


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

function setup_region_routes(config)
    reg_assess_data = regional_assessment_data(config)

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
        Rasters.write(
            mask_path,
            Raster(assess.stack[names(assess.stack)[1]]; data=mask_data, missingval=0);
            ext=".tiff",
            source="gdal",
            driver="COG",
            options=Dict{String,String}(
                "COMPRESS"=>"DEFLATE",
                "SPARSE_OK"=>"TRUE",
                "OVERVIEW_COUNT"=>"5"
            )
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
