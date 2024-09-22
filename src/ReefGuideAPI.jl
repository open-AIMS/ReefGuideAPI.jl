module ReefGuideAPI

using Base.Threads
using
    Glob,
    TOML

using Serialization

using DataFrames
using OrderedCollections
using Memoization
using SparseArrays

import GeoDataFrames as GDF
using
    ArchGDAL,
    GeoParquet,
    Rasters

using
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

function setup_regional_data(config::Dict)
    return setup_regional_data(config["prepped_data"]["PREPPED_DATA_DIR"])
end
function setup_regional_data(reef_data_path::String)
    if @isdefined(REGIONAL_DATA)
        @debug "Using previously generated regional data store."
        sleep(1)  # Pause so message is noticeably visible
        return REGIONAL_DATA
    end

    @debug "Setting up regional data store..."

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

        # Pre-extract long/lat coordinates
        coords = GI.coordinates.(slope_table.geometry)
        slope_table[!, :lons] .= first.(coords)
        slope_table[!, :lats] .= last.(coords)

        coords = GI.coordinates.(flat_table.geometry)
        flat_table[!, :lons] .= first.(coords)
        flat_table[!, :lats] .= last.(coords)

        rst_stack = RasterStack(data_paths; name=data_names, lazy=true)
        regional_assessment_data[reg] = RegionalCriteria(
            rst_stack,
            slope_table,
            flat_table
        )
    end

    # Remember, `@eval` runs in global scope.
    @eval const REGIONAL_DATA = $(regional_assessment_data)

    return REGIONAL_DATA
end

"""
    _cache_location(config::Dict)::String

Retrieve cache location for geotiffs.
"""
function _cache_location(config::Dict)::String
    cache_loc = try
        in_debug = "DEBUG_MODE" in config["server_config"]
        if in_debug && lowercase(config["server_config"]["DEBUG_MODE"]) == "true"
            mktempdir()
        else
            config["server_config"]["TIFF_CACHE_DIR"]
        end
    catch
        mktempdir()
    end

    return cache_loc
end

"""
    n_gdal_threads(config::Dict)::String

Retrieve the configured number of threads to use when writing COGs with GDAL.
"""
function n_gdal_threads(config::Dict)::String
    n_cog_threads = try
        config["server_config"]["COG_THREADS"]
    catch
        "1"  # Default to using a single thread for GDAL write
    end

    return n_cog_threads
end

"""
    tile_size(config::Dict)::Tuple

Retrieve the configured size of map tiles in pixels (width and height / lon and lat).
"""
function tile_size(config::Dict)::Tuple
    tile_dims = try
        res = parse(Int, config["server_config"]["TILE_SIZE"])
        (res, res)
    catch
        (256, 256)  # 256x256
    end

    return tile_dims
end

function start_server(config_path)
    @info "Launching server... please wait"

    @info "Parsing configuration from $(config_path)..."
    config = TOML.parsefile(config_path)
    @info "Successfully parsed configuration."

    @info "Setting up region routes..."
    setup_region_routes(config)
    @info "Completed region routes setup."

    @info "Setting up tile routes..."
    setup_tile_routes()
    @info "Completed tile routes setup."

    @info "Initialisation complete, starting server on port 8000."
    @info "Starting with $(Threads.nthreads()) threads..."
    if Threads.nthreads() > 1
        serveparallel(host="0.0.0.0", port=8000)
    else
        serve(host="0.0.0.0", port=8000)
    end
end

export
    RegionalCriteria,
    criteria_data_map

# Methods to assess/identify deployment "plots" of reef.
export
    assess_reef_site,
    identify_potential_sites_edges,
    filter_intersecting_sites

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
