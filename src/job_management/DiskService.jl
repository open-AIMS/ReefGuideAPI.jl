struct DiskService <: JobService
    cache_dir::String
end

function _job_location(srv::DiskService, job_id::String)
    return joinpath(srv.cache_dir, job_id)
end

"""
    job_details(srv::DiskService, job_id::String)

Retrieve details of a job from disk.

# Arguments
- `srv` : DiskService
- `job_id` : ID of job, typically based on a hash of the search criteria

# Returns
Job details in a dictionary
"""
function job_details(srv::DiskService, job_id::String)::JobAttributes
    d = try
        JSON.parsefile(_job_location(srv, job_id))
    catch err
        if !(err isa SystemError)
            rethrow(err)
        end

        Dict(
            "job_id" => "invalid",
            "status" => "no job",
            "result_loc" => "",
            "access_url" => ""
        )
    end

    return JobAttributes(d)
end

"""
    submit_job(srv::DiskService, job_id::String, result_loc::String)::JobAttributes

Submit a job by writing state to disk.

# Arguments
- `srv` : Service type
- `job_id` : ID of job, typically based on a hash of search criteria
- `result_loc` : Expected location of results
"""
function submit_job(srv::DiskService, job_id::String, result_loc::String)::JobAttributes
    fn = _job_location(srv, job_id)
    attr = JobAttributes(job_id, "processing", result_loc, create_job_url(job_id))
    open(fn, "w") do f
        JSON.print(f, attr)
    end

    return attr
end

function update_job!(srv::DiskService, job_id::String, attrs::JobAttributes)::Nothing
    fn = _job_location(srv, job_id)
    open(fn, "w") do f
        JSON.print(f, attrs)
    end

    return nothing
end

"""
    job_status(srv::DiskService, job_id::String)::String

Retrieve status of a job.

# Arguments
- `srv` : Service type
- `job_id` : ID of job, typically based on a hash of search criteria
"""
function job_status(srv::DiskService, job_id::String)::String
    details = job_details(srv, job_id)

    return details.status
end

function job_status(details::JobAttributes)
    return details.status
end

"""
    job_result(srv::DiskService, job_id::String)::Union{Bool,String}

Retrieve the location of results for a job.
Raises ArgumentError in cases where the job does not exist.

# Arguments
- `srv` : Service type
- `job_id` : ID of job, typically based on a hash of search criteria
"""
function job_result(srv::DiskService, job_id::String)::Union{Bool,String}
    details = job_details(srv, job_id)
    if details.status == "no job"
        throw(ArgumentError("Job does not exist"))
    end

    if details.status != "completed"
        return false
    end

    return details.result_loc
end
