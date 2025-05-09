using JSON3
using Logging
using Dates

#
# Core type definitions
#

"""
Enum for job types matching the API definition
"""
@enum JobType begin
    CRITERIA_POLYGONS
end

symbol_to_job_type = Dict(zip(Symbol.(instances(JobType)), instances(JobType)))

"""
Enum for storage schemes matching the API definition
"""
@enum StorageScheme begin
    S3
    # Add more storage schemes as needed
end

"""
Abstract type for job input payloads
All concrete job input types should inherit from this
"""
abstract type AbstractJobInput end

"""
Abstract type for job output payloads
All concrete job output types should inherit from this
"""
abstract type AbstractJobOutput end

"""
Abstract type for job handler implementations
All concrete job handlers should inherit from this
"""
abstract type AbstractJobHandler end

"""
Registry mapping job types to handlers, input/output types, and validators
"""
struct JobRegistry
    handlers::Dict{JobType,AbstractJobHandler}
    input_types::Dict{JobType,Type{<:AbstractJobInput}}
    output_types::Dict{JobType,Type{<:AbstractJobOutput}}

    function JobRegistry()
        return new(
            Dict{JobType,AbstractJobHandler}(),
            Dict{JobType,Type{<:AbstractJobInput}}(),
            Dict{JobType,Type{<:AbstractJobOutput}}()
        )
    end
end

# Global registry instance
const JOB_REGISTRY = JobRegistry()

#
# Registration functions
#

"""
Register a job handler for a specific job type
"""
function register_job_handler!(
    job_type::JobType,
    handler::AbstractJobHandler,
    input_type::Type{<:AbstractJobInput},
    output_type::Type{<:AbstractJobOutput}
)
    JOB_REGISTRY.handlers[job_type] = handler
    JOB_REGISTRY.input_types[job_type] = input_type
    JOB_REGISTRY.output_types[job_type] = output_type

    @debug "Registered handler for job type: $job_type"
    return nothing
end

"""
Get the appropriate handler for a job type
"""
function get_job_handler(job_type::JobType)::AbstractJobHandler
    if !haskey(JOB_REGISTRY.handlers, job_type)
        error("No handler registered for job type: $job_type")
    end
    return JOB_REGISTRY.handlers[job_type]
end

#
# Validation functions
#

"""
Parse and validate a job input payload
"""
function validate_job_input(job_type::JobType, raw_payload::Any)
    if !haskey(JOB_REGISTRY.input_types, job_type)
        error("No input type registered for job type: $job_type")
    end

    input_type = JOB_REGISTRY.input_types[job_type]

    try
        # Parse the raw JSON payload into the appropriate type
        return JSON3.read(JSON3.write(raw_payload), input_type)
    catch e
        @error "Input validation failed for job type $job_type" exception = (
            e, catch_backtrace()
        )
        error("Invalid input payload for job type: $job_type")
    end
end

"""
Validate a job output payload
"""
function validate_job_output(job_type::JobType, output::AbstractJobOutput)
    if !haskey(JOB_REGISTRY.output_types, job_type)
        error("No output type registered for job type: $job_type")
    end

    expected_type = JOB_REGISTRY.output_types[job_type]

    if !isa(output, expected_type)
        error("Output payload is not of the correct type for job type: $job_type")
    end

    return output
end

#
# Job processing
#

"""
Process a job using the appropriate handler
"""
function process_job(
    job_type::JobType, input_payload::Any, storage_uri::String
)::AbstractJobOutput
    # Get the registered handler
    handler = get_job_handler(job_type)

    # Validate and convert input payload
    typed_input = validate_job_input(job_type, input_payload)

    # Process the job
    @debug "Processing job of type: $job_type"
    output = handle_job(handler, typed_input, storage_uri)

    # Validate output
    validate_job_output(job_type, output)

    return output
end

#
# ================= 
# CRITERIA_POLYGONS 
# ================= 
#

"""
Input payload for CRITERIA_POLYGONS job
"""
struct CriteriaPolygonsInput <: AbstractJobInput
    id::Int64
    # Add more fields as needed to match the TS schema
end

"""
Output payload for CRITERIA_POLYGONS job
"""
struct CriteriaPolygonsOutput <: AbstractJobOutput
end

"""
Handler for CRITERIA_POLYGONS jobs
"""
struct CriteriaPolygonsHandler <: AbstractJobHandler end

"""
Process a CRITERIA_POLYGONS job
"""
function handle_job(
    ::CriteriaPolygonsHandler, input::CriteriaPolygonsInput, storage_uri::String
)::CriteriaPolygonsOutput
    @debug "Processing criteria polygons job with id: $(input.id)"

    # Simulate processing time
    sleep(10)

    @debug "Finished criteria polygons job with id: $(input.id)"
    @debug "Could write something to $(storage_uri) if desired."

    # This is where the actual job processing would happen
    # For now, we just return a dummy output
    return CriteriaPolygonsOutput()
end

#
# ====
# INIT
# ====
#

#
# Register the job types when the module loads
#
function __init__()
    # Register the CRITERIA_POLYGONS job handler
    register_job_handler!(
        CRITERIA_POLYGONS,
        CriteriaPolygonsHandler(),
        CriteriaPolygonsInput,
        CriteriaPolygonsOutput
    )

    @debug "Jobs module initialized with handlers"
end
