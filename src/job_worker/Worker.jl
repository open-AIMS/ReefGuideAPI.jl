using Dates
using Random
using Logging
using JSON3

"""
Represents a job that needs to be processed
"""
struct Job
    id::Int64
    type::String
    # We don't really know what type of data this will be 
    input_payload::Any
end

"""
Represents an assignment for a job
"""
struct JobAssignment
    # ID of the assignment
    id::Int64
    # ID of the tasked job
    job_id::Int64
    # Path to where the data can be stored (in s3)
    storage_uri::String
end

"""
Context provided to job handlers with all necessary information
"""
struct JobContext
    # The job to be processed
    job::Job
    # The job assignment details
    assignment::JobAssignment
    # The API client for making HTTP requests
    http_client::Any
    # Task metadata
    task_metadata::Any

    # Constructor that takes all fields
    function JobContext(job, assignment, http_client, task_metadata)
        return new(job, assignment, http_client, task_metadata)
    end
end

"""
Job handler type definition - a function that processes jobs of a specific type
"""
abstract type JobHandler end

"""
Default job handler that does minimal processing
"""
struct DefaultJobHandler <: JobHandler end

"""
Process a job with the default handler
Return a tuple of (success::Bool, result_payload::Any)
"""
function process(::DefaultJobHandler, context::JobContext)
    # Default implementation - 90% success rate
    success = rand() > 0.1
    sleep(rand(5:15))  # Simulate processing time between 5-15 seconds
    return (success, success ? Dict() : nothing)
end

"""
Handler that processes jobs using the Jobs module handlers
"""
struct TypedJobHandler <: JobHandler end

"""
Process method that uses the Jobs module to handle the job
"""
function process(::TypedJobHandler, context::JobContext)
    try
        # Extract job type from the job
        job_type_str = context.job.type

        # Convert string to JobType enum (safely)
        job_type = try
            getproperty(Jobs, Symbol(job_type_str))
        catch e
            @error "Unknown job type: $job_type_str" exception = (e, catch_backtrace())
            return (false, nothing)
        end

        # Get storage URI from the assignment
        storage_uri = context.assignment.storage_uri

        # Process the job using the Jobs framework
        @info "Processing job $(context.job.id) with type $(job_type_str)"
        output = Jobs.process_job(job_type, context.job.input_payload, storage_uri)

        # Convert output to a dictionary for the worker framework
        result_payload = JSON3.read(JSON3.write(output), Dict)

        # Return success and the result
        return (true, result_payload)
    catch e
        @error "Error processing job: $e" exception = (e, catch_backtrace())
        return (false, nothing)
    end
end

"""
Worker configuration
"""
struct WorkerConfig
    # API connection settings
    api_endpoint::String

    # Worker behavior
    job_types::Vector{String}
    poll_interval_ms::Int
    idle_timeout_ms::Int

    # Auth settings
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
The Worker Service that manages job processing
"""
mutable struct WorkerService
    # Configuration
    config::WorkerConfig

    # Whether the worker is currently running
    is_running::Bool

    # HTTP client for API calls
    http_client::Any

    # Task metadata
    metadata::Any

    # Last activity timestamp
    last_activity_timestamp::DateTime

    # Job handlers registry - maps job types to handler functions
    job_handlers::Dict{String,JobHandler}

    # Constructor
    function WorkerService(config::WorkerConfig, http_client, metadata)
        worker = new(
            config,
            false,
            http_client,
            metadata,
            now(),
            Dict{String,JobHandler}()
        )

        # Automatically register the TypedJobHandler for all supported job types
        register_typed_handlers!(worker)

        return worker
    end
end

"""
Register the TypedJobHandler for all job types supported by the Jobs module
"""
function register_typed_handlers!(worker::WorkerService)
    # Loop through all values in the JobType enum
    for job_type in instances(Jobs.JobType)
        # Convert enum to string
        job_type_str = string(job_type)

        # Register the TypedJobHandler for this job type
        register_handler!(worker, job_type_str, TypedJobHandler())
    end
end

"""
Register a job handler for a specific job type
"""
function register_handler!(worker::WorkerService, job_type::String, handler::JobHandler)
    worker.job_handlers[job_type] = handler
    @info "Registered handler for job type: $job_type"
    return worker
end

"""
Get the appropriate handler for a job type
"""
function get_handler(worker::WorkerService, job_type::String)::JobHandler
    return get(worker.job_handlers, job_type, DefaultJobHandler())
end

"""
Update the last activity timestamp
"""
function update_last_activity!(worker::WorkerService)
    worker.last_activity_timestamp = now()
    return worker
end

"""
Start the worker
"""
function start(worker::WorkerService)
    @info "Starting worker with config:" job_types = worker.config.job_types poll_interval_ms =
        worker.config.poll_interval_ms idle_timeout_ms = worker.config.idle_timeout_ms

    worker.is_running = true

    # Run the main loop in the current thread
    run_worker_loop(worker)

    return worker
end

"""
Stop the worker
"""
function stop(worker::WorkerService)
    @info "Stopping worker..."
    worker.is_running = false
    return worker
end

"""
Main worker loop
"""
function run_worker_loop(worker::WorkerService)
    @info "Starting worker loop"

    while worker.is_running
        try
            # Poll for a job
            job = poll_for_job(worker)

            # Process job if found
            if !isnothing(job)
                process_job_completely(worker, job)
            end

            # Check for idle timeout
            check_idle_timeout(worker)

            # Sleep before next poll
            sleep(worker.config.poll_interval_ms / 1000)
        catch e
            @error "Error in worker loop: $e" exception = (e, catch_backtrace())
            # Sleep briefly before retrying to avoid hammering the API on errors
            sleep(1.0)
        end
    end

    @info "Worker loop ended"
end

"""
Check if worker has been idle too long and should shut down
"""
function check_idle_timeout(worker::WorkerService)
    if worker.config.idle_timeout_ms > 0
        idle_time_ms = Dates.value(now() - worker.last_activity_timestamp) ÷ 1_000_000
        if idle_time_ms >= worker.config.idle_timeout_ms
            @info "Worker idle for $(idle_time_ms)ms, shutting down..."
            worker.is_running = false
        end
    end
end

"""
Poll for a single available job
"""
function poll_for_job(worker::WorkerService)::Union{Job,Nothing}
    try
        # Get available jobs
        response = worker.http_client.get(
            "/jobs/poll"; params=Dict("jobType" => worker.config.job_types[1])
        )

        jobs = response.jobs

        if isempty(jobs)
            return nothing
        end

        # Update activity timestamp when we find a job
        update_last_activity!(worker)

        # Return the first available job
        return jobs[1]
    catch e
        @error "Error polling for jobs: $e" exception = (e, catch_backtrace())
        return nothing
    end
end

"""
Process a job completely (claim, process, complete)
"""
function process_job_completely(worker::WorkerService, job::Job)
    try
        # Try to claim the job
        assignment = claim_job(worker, job)

        if isnothing(assignment)
            @warn "Failed to claim job $(job.id)"
            return nothing
        end

        # Get the appropriate handler for this job type
        handler = get_handler(worker, job.type)

        # Process the job synchronously
        @info "Processing job $(job.id) with handler for type $(job.type)"

        # Create context for the handler
        context = JobContext(job, assignment, worker.http_client, worker.metadata)

        # Process the job with the handler
        success, result_payload = process(handler, context)

        # Complete the job
        complete_job(worker, assignment.id, job, success, result_payload)
    catch e
        @error "Error processing job $(job.id): $e" exception = (e, catch_backtrace())
    end
end

"""
Claim a job and get an assignment
"""
function claim_job(worker::WorkerService, job::Job)::Union{JobAssignment,Nothing}
    try
        # Try to claim the job
        task_arn =
            isnothing(worker.metadata.task_arn) ? "Unknown - metadata lookup failure" :
            worker.metadata.task_arn
        cluster_arn =
            isnothing(worker.metadata.cluster_arn) ? "Unknown - metadata lookup failure" :
            worker.metadata.cluster_arn

        assignment_response = worker.http_client.post(
            "/jobs/assign",
            Dict(
                "jobId" => job.id,
                "ecsTaskArn" => task_arn,
                "ecsClusterArn" => cluster_arn
            )
        )

        assignment = assignment_response.assignment
        @info "Claimed job $(job.id), assignment $(assignment.id)"

        # Update activity timestamp when we claim a job
        update_last_activity!(worker)

        return assignment
    catch e
        @error "Error claiming job $(job.id): $e" exception = (e, catch_backtrace())
        return nothing
    end
end

"""
Complete a job
"""
function complete_job(
    worker::WorkerService,
    assignment_id::Int64,
    job::Job,
    success::Bool,
    result_payload::Any
)
    try
        @info "Completing job $(job.id)"

        worker.http_client.post(
            "/jobs/assignments/$(assignment_id)/result",
            Dict(
                "status" => success ? "SUCCEEDED" : "FAILED",
                "resultPayload" => result_payload
            )
        )

        @info "Job $(job.id) completed with status: $(success ? "SUCCESS" : "FAILURE")"

        # Update activity timestamp when we complete a job
        update_last_activity!(worker)
    catch e
        @error "Error completing job $(job.id): $e" exception = (e, catch_backtrace())
    end
end