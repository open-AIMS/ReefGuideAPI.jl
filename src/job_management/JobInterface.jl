using StructTypes

abstract type JobService end

mutable struct JobAttributes
    const job_id::String

    # Should be an enum
    status::String  # One of "no job", "completed", "processing", "error"
    const result_loc::String
    const access_url::String

    # Init datetime, error message?
end

# Define struct type definition to auto-serialize/deserialize to JSON
StructTypes.StructType(::Type{JobAttributes}) = StructTypes.Struct()

"""
    JobAttributes(attr::Dict)::JobAttributes

Construct JobAttributes from a dictionary.
"""
function JobAttributes(attr::Dict)::JobAttributes
    details::Dict{String,String} = JSON.parse(JSON.json(attr))

    return JobAttributes(
        details["job_id"],
        details["status"],
        details["result_loc"],
        details["access_url"]
    )
end

"""
    create_job_id(query_params::Dict)::String

Generate a job id based on query parameters.
"""
function create_job_id(query_params::Dict)::String
    return string(hash(query_params))
end

"""
    create_job_url(job_id::String)::String

Generate URL endpoint to retrieve results of a job.
"""
function create_job_url(job_id::String)::String
    return "/job/result/$(job_id)"
end

function job_details end
function job_status end
function submit_job end
function update_job! end
function job_result end
