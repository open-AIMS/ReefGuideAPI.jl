using Dates

"""
Configuration for the worker
"""
struct WorkerConfig
    # API connection settings
    api_endpoint::String

    # Worker behavior
    job_types::Vector{String}
    poll_interval_ms::Int
    idle_timeout_ms::Int

    # Auth settings - these target the reefguide-web-api providing necessary
    # credentials for creating jobs etc
    username::String
    password::String

    # Constructor with defaults to handle optional fields
    function WorkerConfig(
        api_endpoint::String,
        job_types::Vector{String},
        username::String,
        password::String;
        poll_interval_ms::Int=1000,
        idle_timeout_ms::Int=2 * 60 * 1000
    )
        return new(
            api_endpoint,
            job_types,
            poll_interval_ms,
            idle_timeout_ms,
            username,
            password
        )
    end
end

"""
Validation error for config values
"""
struct ConfigValidationError <: Exception
    field::String
    message::String
end

Base.showerror(io::IO, e::ConfigValidationError) =
    print(io, "ConfigValidationError: $(e.field) - $(e.message)")

"""
Validates a URL string
"""
function validate_url(url::String, field_name::String)
    # Basic URL validation - could be enhanced with regex
    if !startswith(url, "http://") && !startswith(url, "https://")
        throw(ConfigValidationError(field_name, "Invalid URL format: $url"))
    end
    return url
end

"""
Validates and parses an integer from string
"""
function parse_int(value::String, field_name::String, min_value::Int=nothing)
    try
        parsed = parse(Int, value)
        if !isnothing(min_value) && parsed < min_value
            throw(ConfigValidationError(field_name, "Value must be at least $min_value"))
        end
        return parsed
    catch e
        if isa(e, ArgumentError)
            throw(ConfigValidationError(field_name, "Cannot parse '$value' as integer"))
        end
        rethrow(e)
    end
end

"""
Gets an environment variable with validation
"""
function get_env(key::String, required::Bool=true)
    value = get(ENV, key, nothing)
    if isnothing(value) && required
        throw(ConfigValidationError(key, "Required environment variable not set"))
    end
    return value
end

"""
Load configuration from environment variables
"""
function load_config_from_env()::WorkerConfig
    # Fetch and validate required environment variables
    api_endpoint = get_env("API_ENDPOINT")
    validate_url(api_endpoint, "API_ENDPOINT")

    job_types_str = get_env("JOB_TYPES")
    job_types = strip.(split(job_types_str, ','))
    if isempty(job_types) || any(isempty, job_types)
        throw(
            ConfigValidationError(
                "JOB_TYPES", "At least one non-empty job type must be specified"
            )
        )
    end

    username = get_env("USERNAME")
    if isempty(username)
        throw(ConfigValidationError("USERNAME", "Username cannot be empty"))
    end

    password = get_env("PASSWORD")
    if isempty(password)
        throw(ConfigValidationError("PASSWORD", "Password cannot be empty"))
    end

    # Optional environment variables with defaults
    poll_interval_ms = parse_int(
        something(get_env("POLL_INTERVAL_MS", false), "1000"),
        "POLL_INTERVAL_MS",
        1000
    )

    idle_timeout_ms = parse_int(
        something(get_env("IDLE_TIMEOUT_MS", false), string(2 * 60 * 1000)),
        "IDLE_TIMEOUT_MS",
        1
    )

    # Create and return the config object
    return WorkerConfig(
        api_endpoint,
        job_types,
        username,
        password,
        poll_interval_ms,
        idle_timeout_ms,
    )
end

"""
Create a Worker instance from environment variables
"""
function create_worker_from_env()::WorkerService
    # Load configuration from environment variables
    config = load_config_from_env()

    # Create credentials and API client
    credentials = Credentials(config.username, config.password)
    http_client = AuthApiClient(config.api_endpoint, credentials)

    # Create and return worker instance
    return WorkerService(config, http_client, task_metadata)
end