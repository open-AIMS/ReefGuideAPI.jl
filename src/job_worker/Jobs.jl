"""
This is the file where handlers, input and output payloads are registered to
handle jobs for this worker.
"""

using JSON3
using Logging
using Dates

"""
    create_job_id(query_params::Dict)::String

Generate a job id based on query parameters.
"""
function create_job_id(query_params::Dict)::String
    return string(hash(query_params))
end

"""
Builds a predictable file name based on pulling out the regional assessment
relevant criteria - in the configured cache location.
"""
function build_regional_assessment_file_path(;
    query_params::Dict, region::String, reef_type::String, ext::String,
    config::Dict
)::String
    @debug "Ascertaining file name for regional assessment"
    cache_path = _cache_location(config)

    # Use only the criteria relevant to regional assessment
    regional_assess_criteria = extract_criteria(query_params, suitability_criteria())

    # NOTE: This is where we could add additional hash params to invalidate
    # cache
    job_id = create_job_id(regional_assess_criteria)

    return joinpath(
        cache_path, "$(job_id)_$(region)_$(reef_type)_regional_assessment.$(ext)"
    )
end

# ================
# Type Definitions
# ================

"""
Enum for job types matching the API definition
"""
@enum JobType begin
    SUITABILITY_ASSESSMENT
    REGIONAL_ASSESSMENT
    TEST
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
A context object passed through to a job handler
"""
struct HandlerContext
    "The path to the s3 storage location permitted for writing"
    storage_uri::String
    "The path to the config file"
    config_path::String
    aws_region::String

    function HandlerContext(;
        storage_uri::String, config_path::String, aws_region::String="ap-southeast-2"
    )
        return new(storage_uri, config_path, aws_region)
    end
end

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

# ======================
# Registration functions
# ======================

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
    job_type::JobType, input_payload::Any, context::HandlerContext
)::AbstractJobOutput
    # Get the registered handler
    handler = get_job_handler(job_type)

    # Validate and convert input payload
    typed_input = validate_job_input(job_type, input_payload)

    # Process the job
    @debug "Processing job of type: $job_type"
    output = handle_job(handler, typed_input, context)

    # Validate output
    validate_job_output(job_type, output)

    return output
end

#
# ================= 
# TEST 
# ================= 
#

"""
Input payload for TEST job
"""
struct TestInput <: AbstractJobInput
    id::Int64
end

"""
Output payload for TEST job
"""
struct TestOutput <: AbstractJobOutput
end

"""
Handler for TEST jobs
"""
struct TestHandler <: AbstractJobHandler end

"""
Process a TEST job
"""
function handle_job(
    ::TestHandler, input::TestInput, context::HandlerContext
)::TestOutput
    @debug "Processing test job with id: $(input.id)"

    # Simulate processing time
    sleep(10)

    @debug "Finished test job with id: $(input.id)"
    @debug "Could write something to $(context.storage_uri) if desired."

    # This is where the actual job processing would happen
    # For now, we just return a dummy output
    return TestOutput()
end

#
# ===================
# REGIONAL_ASSESSMENT 
# ===================
#

"""
Input payload for REGIONAL_ASSESSMENT job

Subset of CRITERIA_ASSESSMENT payload
"""
struct RegionalAssessmentInput <: AbstractJobInput
    # High level config
    "Region for assessment"
    region::String
    "The type of reef, slopes or flats"
    reef_type::String

    # Criteria
    "The depth range (min)"
    depth_min::Float64
    "The depth range (max)"
    depth_max::Float64
    "The slope range (min)"
    slope_min::Float64
    "The slope range (max)"
    slope_max::Float64
    "The rugosity range (min)"
    rugosity_min::Float64
    "The rugosity range (max)"
    rugosity_max::Float64
    "Suitability threshold (min)"
    threshold::Int64
end

"""
Output payload for REGIONAL_ASSESSMENT job
"""
struct RegionalAssessmentOutput <: AbstractJobOutput
    cog_path::String
end

"""
Handler for REGIONAL_ASSESSMENT jobs
"""
struct RegionalAssessmentHandler <: AbstractJobHandler end

"""
Handler for the regional assessment job. 
"""
function handle_job(
    ::RegionalAssessmentHandler, input::RegionalAssessmentInput,
    context::HandlerContext
)::RegionalAssessmentOutput
    @info "Initiating regional assessment task"

    @info "Parsing configuration from $(context.config_path)..."
    config = TOML.parsefile(context.config_path)
    @info "Configuration parsing complete."

    @info "Setting up regional assessment data"
    reg_assess_data = setup_regional_data(config)
    @info "Done setting up regional assessment data"

    @info "Performing regional assessment"

    # Pull out these parameters in the format previously expected 
    reg = input.region
    rtype = input.reef_type

    # This is the format expected by prior methods TODO improve the separation
    # of concerns between query string parameters and function inputs - the API
    # routing concerns should not be exposed to the application layer
    qp::Dict{String,String} = Dict{String,String}(
        "Depth" => "$(input.depth_min):$(input.depth_max)",
        "Slope" => "$(input.slope_min):$(input.slope_max)",
        "Rugosity" => "$(input.rugosity_min):$(input.rugosity_max)",
        "SuitabilityThreshold" => "$(input.threshold)"
    )

    assessed_fn = build_regional_assessment_file_path(;
        query_params=qp, region=reg, reef_type=rtype, ext="tiff", config
    )
    @debug "COG File name: $(assessed_fn)"

    if !isfile(assessed_fn)
        @debug "File system cache was not hit for this task"
        @debug "Assessing region $(reg)"
        assessed = assess_region(reg_assess_data, reg, qp, rtype)

        @debug now() "Writing COG of regional assessment to $(assessed_fn)"
        _write_cog(assessed_fn, assessed, config)
        @debug now() "Finished writing cog "
    else
        @info "Cache hit - skipping regional assessment process and re-uploading to output!"
    end

    # Now upload this to s3 
    client = S3StorageClient(; region=context.aws_region)

    # Output file names
    output_file_name_rel = "regional_assessment.tiff"
    full_s3_target = "$(context.storage_uri)/$(output_file_name_rel)"
    @debug "File paths:" relative = output_file_name_rel absolute = full_s3_target

    @debug now() "Initiating file upload"
    upload_file(client, assessed_fn, full_s3_target)
    @debug now() "File upload completed"

    @debug "Finished regional assessment job."
    return RegionalAssessmentOutput(
        output_file_name_rel
    )
end

#
# ======================
# SUITABILITY_ASSESSMENT 
# ======================
#

"""
Input payload for SUITABILITY_ASSESSMENT job
"""
struct SuitabilityAssessmentInput <: AbstractJobInput
    # High level config

    "Region for assessment"
    region::String
    "The type of reef, slopes or flats"
    reef_type::String

    # Criteria
    "The depth range (min)"
    depth_min::Float64
    "The depth range (max)"
    depth_max::Float64
    "The slope range (min)"
    slope_min::Float64
    "The slope range (max)"
    slope_max::Float64
    "The rugosity range (min)"
    rugosity_min::Float64
    "The rugosity range (max)"
    rugosity_max::Float64
    "Suitability threshold (min)"
    threshold::Int64
    "Length dimension of target polygon"
    x_dist::Int64
    "Width dimension of target polygon"
    y_dist::Int64
end

"""
Output payload for SUITABILITY_ASSESSMENT job
"""
struct SuitabilityAssessmentOutput <: AbstractJobOutput
    geojson_path::String
end

"""
Handler for SUITABILITY_ASSESSMENT jobs
"""
struct SuitabilityAssessmentHandler <: AbstractJobHandler end

"""
Handler for the suitability assessment job. 
"""
function handle_job(
    ::SuitabilityAssessmentHandler, input::SuitabilityAssessmentInput,
    context::HandlerContext
)::SuitabilityAssessmentOutput
    @info "Initiating site assessment task"

    @info "Parsing configuration from $(context.config_path)..."
    config = TOML.parsefile(context.config_path)
    @info "Configuration parsing complete."

    @info "Setting up regional assessment data"
    reg_assess_data = setup_regional_data(config)
    @info "Done setting up regional assessment data"

    @info "Performing regional assessment (dependency of site assessment)"

    # Pull out these parameters in the format previously expected 
    reg = input.region
    rtype = input.reef_type

    # This is the format expected by prior methods TODO improve the separation
    # of concerns between query string parameters and function inputs - the API
    # routing concerns should not be exposed to the application layer
    qp::Dict{String,String} = Dict{String,String}(
        "Depth" => "$(input.depth_min):$(input.depth_max)",
        "Slope" => "$(input.slope_min):$(input.slope_max)",
        "Rugosity" => "$(input.rugosity_min):$(input.rugosity_max)",
        "SuitabilityThreshold" => "$(input.threshold)",
        "xdist" => "$(input.x_dist)",
        "ydist" => "$(input.y_dist)"
    )

    assessed_fn = build_regional_assessment_file_path(;
        query_params=qp, region=reg, reef_type=rtype, ext="tiff", config
    )
    @debug "COG File name: $(assessed_fn)"

    if !isfile(assessed_fn)
        @debug "File system cache was not hit for this task"
        @debug "Assessing region $(reg)"
        assessed = assess_region(reg_assess_data, reg, qp, rtype)

        @debug "Writing COG to $(assessed_fn)"
        _write_cog(assessed_fn, assessed, config)
    else
        @info "Cache hit - skipping regional assessment process!"
        @debug "Pulling out raster from cache"
        assessed = Raster(assessed_fn; missingval=0, lazy=true)
    end

    # Extract criteria and assessment
    pixel_criteria = extract_criteria(qp, search_criteria())
    deploy_site_criteria = extract_criteria(qp, site_criteria())

    @debug "Performing site assessment"
    best_sites = filter_sites(
        assess_sites(
            reg_assess_data, reg, rtype, pixel_criteria, deploy_site_criteria,
            assessed
        )
    )

    # Specifically clear from memory to invoke garbage collector
    assessed = nothing

    @debug "Writing to temporary file"
    geojson_name = "$(tempname()).geojson"
    @debug "File name $(geojson_name)"

    if nrow(best_sites) == 0
        open(geojson_name, "w") do f
            JSON.print(f, nothing)
        end
    else
        output_geojson(geojson_name, best_sites)
    end

    # Now upload this to s3 
    client = S3StorageClient(; region=context.aws_region)

    # Output file names
    output_file_name_rel = "suitable.geojson"
    full_s3_target = "$(context.storage_uri)/$(output_file_name_rel)"
    @debug "File paths:" relative = output_file_name_rel absolute = full_s3_target

    upload_file(client, geojson_name, full_s3_target)

    # clean up temp file
    if isfile(geojson_name)
        @debug "Cleaned up temp file"
        rm(geojson_name)
    end

    @debug "Finished suitability assessment job."
    return SuitabilityAssessmentOutput(
        output_file_name_rel
    )
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
    # Register the TEST job handler
    register_job_handler!(
        TEST,
        TestHandler(),
        TestInput,
        TestOutput
    )

    # Register the SUITABILITY_ASSESSMENT job handler
    register_job_handler!(
        SUITABILITY_ASSESSMENT,
        SuitabilityAssessmentHandler(),
        SuitabilityAssessmentInput,
        SuitabilityAssessmentOutput
    )

    # Register the REGIONAL_ASSESSMENT job handler
    register_job_handler!(
        REGIONAL_ASSESSMENT,
        RegionalAssessmentHandler(),
        RegionalAssessmentInput,
        RegionalAssessmentOutput
    )

    @debug "Jobs module initialized with handlers"
end
