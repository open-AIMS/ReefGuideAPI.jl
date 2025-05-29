module ReefGuideAPI

import GeoDataFrames as GDF
using
    Base.Threads,
    Glob,
    TOML,
    ArchGDAL,
    GeoParquet,
    Rasters,
    HTTP,
    Oxygen,
    Serialization,
    DataFrames,
    OrderedCollections,
    Memoization,
    SparseArrays,
    FLoops, ThreadsX

# Utilities and helpers for assessments
include("utility/index.jl")

# Assessment logic
include("assessment_methods/index.jl")

# Worker system
include("job_worker/index.jl")

function start_server(config_path)
    @info "Launching server... please wait"

    @info "Parsing configuration from $(config_path)..."
    config = TOML.parsefile(config_path)

    @info "Initialising regional data and setting up tile cache"
    initialise_data_with_cache(config)

    @info "Setting up auth middleware and router."
    auth = get_auth_router(config)

    @info "Setting up job routes..."
    setup_job_routes(config, auth)

    @info "Setting up utility routes..."
    setup_utility_routes(config, auth)

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

"""
Create and initialize a worker from the environment.

This is a blocking operation until the worker times out.
"""
function start_worker()
    @info "Initializing worker from environment variables..."
    worker = create_worker_from_env()
    @info "Parsing TOML config"
    config = TOML.parsefile(worker.config.config_path)
    @info "Loading regional data"
    initialise_data_with_cache(config)
    @info "Starting worker loop from ReefGuideAPI.jl"
    start(worker)
    @info "Worker closed itself..."
end

export start_worker, start_server

end
