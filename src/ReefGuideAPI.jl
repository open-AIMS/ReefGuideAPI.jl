module ReefGuideAPI

using Base.Threads
using
    Glob,
    TOML

using DataFrames
using OrderedCollections
using Memoization
using SparseArrays

import GeoDataFrames as GDF
using
    GeoParquet,
    Rasters

using
    FLoops,
    HTTP,
    Oxygen

include("assessment/criteria.jl")
include("geom_handlers/site_assessment.jl")
include("assessment/query_thresholds.jl")


function get_regions()
    # TODO: Comes from config?
    regions = String[
        "Townsville-Whitsunday",
        "Cairns-Cooktown",
        "Mackay-Capricorn",
        "FarNorthern"
    ]

    return regions
end

function setup_regional_data(reef_data_path::String)
    regional_assessment_data = OrderedDict{String, RegionalCriteria}()
    for reg in get_regions()
        data_paths = String[]
        data_names = String[]

        slope_table = nothing
        flat_table = nothing

        for (k, dp) in criteria_data_map()
            g = glob("$reg*$dp.tif", reef_data_path)
            if length(g) == 0
                continue
            end

            push!(data_paths, first(g))
            push!(data_names, string(k))
            if occursin("valid", string(dp))
                # Load up Parquet files
                parq_file = replace(first(g), ".tif"=>"_lookup.parq")

                if occursin("slope", string(dp))
                    slope_table = GeoParquet.read(parq_file)
                elseif occursin("flat", string(dp))
                    flat_table = GeoParquet.read(parq_file)
                else
                    error("Unknown lookup found: $(parq_file)")
                end
            end
        end

        rst_stack = RasterStack(data_paths; name=data_names, lazy=true)
        regional_assessment_data[reg] = RegionalCriteria(
            rst_stack,
            slope_table,
            flat_table
        )
    end

    return regional_assessment_data
end

function regional_assessment_data(config)
    return setup_regional_data(config["prepped_data"]["PREPPED_DATA_DIR"])
end

function _cache_location(config)
    cache_loc = try
        in_debug = "DEBUG_MODE" in config["server_config"]
        if in_debug && lowercase(config["server_config"]["DEBUG_MODE"]) == "true"
            mktempdir()
        else
            config["server_config"]["CACHE_DIR"]
        end
    catch
        mktempdir()
    end

    return cache_loc
end

function n_gdal_threads(config)::String
    n_cog_threads = try
        config["server_config"]["COG_THREADS"]
    catch
        "1"
    end

    return n_cog_threads
end

function start_server(config_path)
    config = TOML.parsefile(config_path)
    setup_region_routes(config)

    if Threads.nthreads() > 1
        serveparallel(port=8000)
    else
        serve(port=8000)
    end
end

export
    RegionalCriteria,
    criteria_data_map,
    regional_assessment_data

# Methods to assess/identify deployment "plots" of reef.
export
    assess_reef_site,
    identify_potential_sites

# Geometry handling
export
    create_poly,
    create_bbox,
    port_buffer_mask

# Raster->Index interactions (defunct?)
export
    valid_slope_lon_inds,
    valid_slope_lat_inds,
    valid_flat_lon_inds,
    valid_flat_lat_inds

end
