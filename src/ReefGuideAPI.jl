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

include("Middleware.jl")
include("admin.jl")
include("file_io.jl")
include("server_cache.jl")

include("job_management/JobInterface.jl")
include("job_management/DiskService.jl")

include("criteria_assessment/criteria.jl")
include("criteria_assessment/query_thresholds.jl")
include("criteria_assessment/regional_assessment.jl")

include("site_assessment/common_functions.jl")
include("site_assessment/best_fit_polygons.jl")

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

function get_auth_router(config::Dict)
    # Setup auth middleware - depends on config.toml - can return identity func
    auth = setup_jwt_middleware(config)
    return router(""; middleware=[auth])
end

function start_server(config_path)
    @info "Launching server... please wait"

    warmup_cache(config_path)

    @info "Parsing configuration from $(config_path)..."
    config = TOML.parsefile(config_path)

    @info "Setting up auth middleware and router."
    auth = get_auth_router(config)

    @info "Setting up criteria routes..."
    setup_criteria_routes(config, auth)

    @info "Setting up region routes..."
    setup_region_routes(config, auth)

    @info "Setting up tile routes..."
    setup_tile_routes(config, auth)

    @info "Setting up job routes..."
    setup_job_routes(config, auth)

    @info "Setting up admin routes..."
    setup_admin_routes(config)

    # Which host should we listen on?
    host = config["server_config"]["HOST"]
    # Which port should we listen on?
    port = 8000

    @info "Initialisation complete, starting server listening on host: $(host) at port $(port) with $(Threads.nthreads()) threads."

    return serve(;
        middleware=[CorsMiddleware],
        host=host,
        port=port,
        parallel=Threads.nthreads() > 1,
        is_prioritized=(req::HTTP.Request) -> req.target == "/health"
    )
end

function runner_loop(config)
    # What endpoint are we using?
    job_endpoint = config["job_runner"]["JOB_API_ENDPOINT"]

    @info "Job runner awaiting jobs at job_endpoint $job_endpoint"
end

function start_runner(config_path)
    @info "Launching job runner... please wait"

    @info "Warming up cache"
    warmup_cache(config_path)

    @info "Parsing configuration from $(config_path)..."
    config = TOML.parsefile(config_path)

    @info "Setting up regional assessment data"
    reg_assess_data = setup_regional_data(config)

    @info "Initialisation complete. Initiating job runner."
    runner_loop(config)
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

end
