struct Job
    id::Int64
    type::String
    # We don't really know what type of data this will be 
    input_payload::Any
end

struct JobAssignment
    # ID of the assignment
    id:Int64
    # ID of the tasked job
    job_id:Int64
    # Path to where the data can be stored (in s3)
    storage_uri:string
end

struct WorkerConfig
    # API connection settings
    api_endpoint::String

    # Worker behavior
    job_types::Vector{String}
    poll_interval_ms::Int
    max_concurrent_jobs::Int
    idle_timeout_ms::Int

    # Auth settings
    username::String
    password::String

    # Constructor with defaults to handle optional fields
    function Config(
        api_endpoint::String,
        job_types::Vector{String},
        username::String,
        password::String;
        poll_interval_ms::Int=1000,
        max_concurrent_jobs::Int=1,
        idle_timeout_ms::Int=2 * 60 * 1000,
    )
        return new(
            api_endpoint,
            job_types,
            poll_interval_ms,
            max_concurrent_jobs,
            idle_timeout_ms,
            username,
            password
        )
    end
end

mutable struct Worker
    config::WorkerConfig
    activeJobs::Dict{Int,Timer}
    isPolling::Bool
    httpClient::AuthApiClient
    metadata::TaskIdentifiers
    idleTimeout::Union{Timer,Nothing}
    lastActivityTimestamp::DateTime

    function TestWorker(
        config::WorkerConfig, httpClient::AuthApiClient, metadata::TaskIdentifiers
    )
        return new(
            config,
            Dict{Int,Timer}(),
            false,
            httpClient,
            metadata,
            nothing,
            now()
        )
    end
end
