"""
The worker service class which manages connecting to all the job worker
components to orchestrate consuming jobs. Principally polls for, completes and
reports back jobs done, on a loop, until idle for a configurable idle time.
"""

using Dates
using Random
using Logging
using JSON3

include("Config.jl")
include("ECS.jl")
include("HttpClient.jl")
include("Jobs.jl")
include("Storage.jl")

"""
Represents a job that needs to be processed
"""
struct Job
    "Job ID in the DB"
    id::Int64
    "What type of Job - strict set of enums"
    type::String
    "The payload defining the job parameters - should correspond to input payload type registered in Jobs.jl"
    input_payload::Any
end

"""
Represents an assignment for a job
"""
struct JobAssignment
    "ID of the assignment in the DB"
    id::Int64
    "ID of the tasked job"
    job_id::Int64
    "Path to where the data can be stored (in s3)"
    storage_uri::String
end

"""
API Response after an assignment
"""
struct JobAssignmentResponse
    assignment::JobAssignment
end

"""
Context provided to job handlers with all necessary information
"""
struct JobContext
    "The job to be processed"
    job::Job
    "The job assignment details"
    assignment::JobAssignment
    "The API client for making HTTP requests"
    http_client::AuthApiClient
    "Task metadata"
    task_metadata::Any
    "Path to the config file"
    config_path::String
    "AWS region for s3 storage"
    aws_region::String
    "Constructor that takes all fields"
    function JobContext(;
        job::Job,
        assignment::JobAssignment,
        http_client::AuthApiClient,
        task_metadata::Any,
        config_path::String,
        aws_region::String="ap-southeast-2"
    )
        return new(job, assignment, http_client, task_metadata, config_path, aws_region)
    end
end

"""
Job handler type definition - a function that processes jobs of a specific type
"""
abstract type JobHandler end

"""
Mock job handler that does minimal processing
"""
struct MockJobHandler <: JobHandler end

"""
Process a job with the mock handler
Return a tuple of (success::Bool, result_payload::Any)
"""
function process(::MockJobHandler, context::JobContext)
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
        job_type = symbol_to_job_type[Symbol(job_type_str)]
        if isnothing(job_type)
            @error "Unknown job type: $job_type_str" exception = (e, catch_backtrace())
            return (false, nothing)
        end

        # Get storage URI from the assignment
        storage_uri = context.assignment.storage_uri

        # Process the job using the Jobs framework
        @debug "Processing job $(context.job.id) with type $(job_type_str)"
        output::AbstractJobOutput = process_job(
            job_type, context.job.input_payload,
            HandlerContext(;
                storage_uri=storage_uri,
                config_path=context.config_path,
                aws_region=context.aws_region
            )
        )

        @debug "Result from process_job $(output)"

        # Convert output to a dictionary for the worker framework If not
        # possible then try to just send empty payload
        result_payload::Dict = Dict()
        try
            result_payload = JSON3.read(JSON3.write(output), Dict)
        catch
            @debug "Error occurred while trying to convert a task output payload into a dictionary - presumably this is due to an empty result payload struct."
        end

        @debug "Parsed payload $(result_payload)"

        # Return success and the result
        return (true, result_payload)
    catch e
        @error "Error processing job: $e" exception = (e, catch_backtrace())
        return (false, nothing)
    end
end

"""
The Worker Service that manages job processing
"""
mutable struct WorkerService
    "Configuration"
    config::WorkerConfig

    "Whether the worker is currently running"
    is_running::Bool

    "HTTP client for API calls"
    http_client::Any

    "Task identifiers and other metadata about this running service"
    metadata::TaskIdentifiers

    "When was the last time this worker did something (got a job/finished it?)"
    last_activity_timestamp::DateTime

    "Job handlers registry - maps job types to handler functions"
    job_handlers::Dict{String,JobHandler}

    function WorkerService(config::WorkerConfig, http_client, identifiers; mock::Bool=false)
        worker = new(
            config,
            false,
            http_client,
            identifiers,
            now(),
            Dict{String,JobHandler}()
        )

        # Automatically register the TypedJobHandler for all supported job types
        if mock
            register_mock_handlers!(worker)
        else
            register_typed_handlers!(worker)
        end

        return worker
    end
end

"""
Register the MockJobHandler for all job types supported by the Jobs module
"""
function register_mock_handlers!(worker::WorkerService)
    # Loop through all values in the JobType enum
    for job_type in instances(JobType)
        # Convert enum to string
        job_type_str = string(job_type)

        # Register the TypedJobHandler for this job type
        register_handler!(worker, job_type_str, MockJobHandler())
    end
end

"""
Register the TypedJobHandler for all job types supported by the Jobs module
"""
function register_typed_handlers!(worker::WorkerService)
    # Loop through all values in the JobType enum
    for job_type in instances(JobType)
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
    @debug "Registered handler for job type: $job_type"
    return worker
end

"""
Get the appropriate handler for a job type
"""
function get_handler(worker::WorkerService, job_type::String)::JobHandler
    return Base.get(worker.job_handlers, job_type, MockJobHandler())
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
    @debug "Starting worker with config:" job_types = worker.config.job_types poll_interval_ms =
        worker.config.poll_interval_ms idle_timeout_ms = worker.config.idle_timeout_ms

    worker.is_running = true

    # Run the main loop in the current thread
    return run_worker_loop(worker)
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
            @debug "Polling for a job"
            job = poll_for_job(worker)

            # Process job if found
            if !isnothing(job)
                @debug "Found a job"
                process_job_completely(worker, job)
            end

            @debug "Checking for idle timeout"
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
        idle_time_ms = floor(
            Int64, Dates.value(now() - worker.last_activity_timestamp)
        )
        @debug "Idle time (milliseconds) $(idle_time_ms). Configured idle timeout is $(worker.config.idle_timeout_ms)"
        if idle_time_ms >= worker.config.idle_timeout_ms
            @debug "Worker idle for $(idle_time_ms)ms, shutting down..."
            worker.is_running = false
        end
    end
end

function poll_for_job(worker::WorkerService)::Union{Job,Nothing}
    try
        # Get available jobs 
        # TODO if we only handle specific subsets, might want to filter more here
        response = HTTPGet(
            worker.http_client, "/jobs/poll";
        )

        @debug "Response from jobs poll: $(response)"

        # Check if we have jobs in the response
        if isempty(response.jobs)
            @debug "No jobs available in response"
            return nothing
        end

        jobs = response.jobs

        # Return the first available job which is of a type that we can handle
        for (i, job_data) in enumerate(jobs)
            @debug "Processing Job[$(i)] = $(job_data)"

            try
                # Parse the job data
                parsed = JSON3.read(JSON3.write(job_data), Job)
                @debug "Result of parsing: $(parsed)"

                # Check if this job type is one we can handle
                if parsed.type in worker.config.job_types
                    @info "Found suitable job of type $(parsed.type)"

                    # Update activity timestamp when we find potential jobs
                    update_last_activity!(worker)

                    return parsed
                else
                    @debug "Skipping job $(i) of type $(parsed.type) (not in our supported types)"
                end
            catch e
                # Handle malformed jobs by logging and continuing to the next one
                @warn "Skipping malformed job at index $(i): $e" exception = (
                    e, catch_backtrace()
                )
                continue
            end
        end

        # If we get here, we found no suitable jobs
        @debug "No suitable jobs found among $(length(jobs)) available jobs"
        return nothing
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
        context = JobContext(;
            job=job,
            assignment=assignment,
            http_client=worker.http_client,
            task_metadata=worker.metadata,
            config_path=worker.config.config_path,
            aws_region=worker.config.aws_region
        )

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

        assignment_response = HTTPPost(
            worker.http_client,
            "/jobs/assign",
            Dict(
                "jobId" => job.id,
                "ecsTaskArn" => task_arn,
                "ecsClusterArn" => cluster_arn
            )
        )

        @debug "Assignment response $(assignment_response)"

        if isnothing(assignment_response)
            @error "Failed job assignment, there was no response."
            return nothing
        end

        # parse into assignment
        assignment_data::JobAssignmentResponse = JSON3.read(
            JSON3.write(assignment_response), JobAssignmentResponse
        )

        assignment = assignment_data.assignment

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

        HTTPPost(
            worker.http_client,
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
