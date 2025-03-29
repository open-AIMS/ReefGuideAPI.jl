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

using FLoops, ThreadsX

import GeoDataFrames as GDF
using
    ArchGDAL,
    GeoParquet,
    Rasters

using
    HTTP,
    Oxygen

include("criteria_assessment/criteria.jl")
include("criteria_assessment/query_thresholds.jl")

include("site_assessment/common_functions.jl")
include("site_assessment/best_fit_polygons.jl")

include("Middleware.jl")
include("admin.jl")

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

"""
    initialize_regional_data_cache(reef_data_path::String, reg_cache_fn::String)

Create initial regional data store with data from `reef_data_path`, excluding geospatial
data and save to `reg_cache_fn` path.
"""
function initialize_regional_data_cache(reef_data_path::String, reg_cache_fn::String)
    regional_assessment_data = OrderedDict{
        String,Union{RegionalCriteria,DataFrame,Dict}
    }()
    for reg in get_regions()
        @debug "$(now()) : Initializing cache for $reg"
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
                parq_file = replace(first(g), ".tif" => "_lookup.parq")

                if occursin("slope", string(dp))
                    slope_table = GeoParquet.read(parq_file)
                elseif occursin("flat", string(dp))
                    @warn "Skipping data for reef flats as it is currently unused"
                    # flat_table = GeoParquet.read(parq_file)
                else
                    msg = "Unknown lookup found: $(parq_file). Must be 'slope' or 'flat'"
                    throw(ArgumentError(msg))
                end
            end
        end

        # Pre-extract long/lat coordinates
        coords = GI.coordinates.(slope_table.geometry)
        slope_table[!, :lons] .= first.(coords)
        slope_table[!, :lats] .= last.(coords)

        # coords = GI.coordinates.(flat_table.geometry)
        # flat_table[!, :lons] .= first.(coords)
        # flat_table[!, :lats] .= last.(coords)

        rst_stack = RasterStack(data_paths; name=data_names, lazy=true)

        # Constrain to just the areas with valid data (with a 0.05 degree buffer)
        # min_lon = min(minimum(slope_table.lons), minimum(flat_table.lons)) - 0.05
        # max_lon = max(maximum(slope_table.lons), maximum(flat_table.lons)) + 0.05
        # min_lat = min(minimum(slope_table.lats), minimum(flat_table.lats)) - 0.05
        # max_lat = max(maximum(slope_table.lats), maximum(flat_table.lats)) + 0.05
        # rst_stack = view(rst_stack, X(min_lon .. max_lon), Y(min_lat .. max_lat))

        regional_assessment_data[reg] = RegionalCriteria(
            rst_stack,
            slope_table,
            slope_table[[1], :]  # flat_table
        )

        @debug "$(now()) : Finished initialization for $reg"
    end

    regional_assessment_data["region_long_names"] = Dict(
        "FarNorthern" => "Far Northern Management Area",
        "Cairns-Cooktown" => "Cairns/Cooktown Management Area",
        "Townsville-Whitsunday" => "Townsville/Whitsunday Management Area",
        "Mackay-Capricorn" => "Mackay/Capricorn Management Area"
    )

    # Store cache on disk to avoid excessive cold startup times
    @debug "Saving regional data cache to disk"
    serialize(reg_cache_fn, regional_assessment_data)

    return regional_assessment_data
end

"""
    setup_regional_data(config::Dict)

Load regional data to act as an in-memory cache.

# Arguments
- `config` : Configuration settings, typically loaded from a TOML file.
- `reef_data_path` : Path to pre-prepared reef data

# Returns
OrderedDict of `RegionalCriteria` for each region.
"""
function setup_regional_data(config::Dict)
    reef_data_path = config["prepped_data"]["PREPPED_DATA_DIR"]
    reg_cache_dir = config["server_config"]["REGIONAL_CACHE_DIR"]
    reg_cache_fn = joinpath(reg_cache_dir, "regional_cache.dat")

    if @isdefined(REGIONAL_DATA)
        @debug "Using previously generated regional data store."
    elseif isfile(reg_cache_fn)
        @debug "Loading regional data cache from disk"
        # Updates to packages like DiskArrays can break deserialization
        try
            @eval const REGIONAL_DATA = deserialize($(reg_cache_fn))
        catch err
            @warn "Failed to deserialize $(reg_cache_fn) with error:" err
            rm(reg_cache_fn)
        end
    end

    if !@isdefined(REGIONAL_DATA)
        @debug "Setting up regional data store..."
        regional_assessment_data = initialize_regional_data_cache(
            reef_data_path,
            reg_cache_fn
        )
        # Remember, `@eval` runs in global scope.
        @eval const REGIONAL_DATA = $(regional_assessment_data)
    end

    # If REGIONAL_DATA is defined, but failed to load supporting data that cannot be
    # cached to disk, such as the reef outlines, (e.g., due to incorrect config), then it
    # will cause errors later on.
    # Then there's no way to address this, even between web server sessions, as `const`
    # values cannot be modified.
    # Here, we check for existence and try to load again if needed.
    if !haskey(REGIONAL_DATA, "reef_outlines")
        reef_outline_path = joinpath(reef_data_path, "rrap_canonical_outlines.gpkg")
        REGIONAL_DATA["reef_outlines"] = GDF.read(reef_outline_path)
    end

    return REGIONAL_DATA
end

"""
    _cache_location(config::Dict)::String

Retrieve cache location for geotiffs.
"""
function _cache_location(config::Dict)::String
    cache_loc = try
        in_debug = haskey(config["server_config"], "DEBUG_MODE")
        if in_debug && lowercase(config["server_config"]["DEBUG_MODE"]) == "true"
            if "DEBUG_CACHE_DIR" ∉ keys(ENV)
                ENV["DEBUG_CACHE_DIR"] = mktempdir()
            end

            ENV["DEBUG_CACHE_DIR"]
        else
            config["server_config"]["TIFF_CACHE_DIR"]
        end
    catch err
        @warn "Encountered error:" err
        if "DEBUG_CACHE_DIR" ∉ keys(ENV)
            ENV["DEBUG_CACHE_DIR"] = mktempdir()
        end

        ENV["DEBUG_CACHE_DIR"]
    end

    return cache_loc
end

"""
    cache_filename(qp::Dict, config::Dict, suffix::String, ext::String)

Generate a filename for a cache.

# Arguments
- `qp` : Query parameters to hash
- `config` : app configuration (to extract cache parent directory from)
- `suffix` : a suffix to use in the filename (pass `""` if none required)
- `ext` : file extension to use
"""
function cache_filename(qp::Dict, config::Dict, suffix::String, ext::String)
    file_id = string(hash(qp))
    temp_path = _cache_location(config)
    cache_file_path = joinpath(temp_path, "$(file_id)$(suffix).$(ext)")

    return cache_file_path
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

function get_auth_router(config::Dict)
    # Setup auth middleware - depends on config.toml - can return identity func
    auth = setup_jwt_middleware(config)
    return router(""; middleware=[auth])
end

"""
    warmup_cache(config_path::String)

Invokes warm up of regional data cache to reduce later spin up times.
"""
function warmup_cache(config_path::String)
    config = TOML.parsefile(config_path)

    # Create re-usable empty tile
    no_data_path = cache_filename(Dict("no_data" => "none"), config, "no_data", "png")
    if !isfile(no_data_path)
        save(no_data_path, zeros(RGBA, tile_size(config)))
    end

    return setup_regional_data(config)
end

function start_server(config_path)
    @info "Launching server... please wait"

    warmup_cache(config_path)

    @info "Parsing configuration from $(config_path)..."
    config = TOML.parsefile(config_path)

    @info "Setting up auth middleware and router."
    auth = get_auth_router(config)

    @info "Setting up region routes..."
    setup_region_routes(config, auth)

    @info "Setting up tile routes..."
    setup_tile_routes(config, auth)

    @info "Setting up admin routes..."
    setup_admin_routes(config)

    port = 8000
    @info "Initialisation complete, starting server on port $(port) with $(Threads.nthreads()) threads."

    return serve(;
        middleware=[CorsMiddleware],
        host="0.0.0.0",
        port=port,
        parallel=Threads.nthreads() > 1,
        is_prioritized=(req::HTTP.Request) -> req.target == "/health"
    )
end

export
    RegionalCriteria,
    criteria_data_map

# Methods to assess/identify deployment "plots" of reef.
export
    assess_reef_site,
    identify_edge_aligned_sites,
    filter_sites,
    output_geojson

# Geometry handling
export
    create_poly,
    create_bbox,
    port_buffer_mask,
    meters_to_degrees,
    polygon_to_lines

# Raster->Index interactions (defunct?)
export
    valid_slope_lon_inds,
    valid_slope_lat_inds,
    valid_flat_lon_inds,
    valid_flat_lat_inds

# Ruleset thresholds
export
    within_thresholds

end
