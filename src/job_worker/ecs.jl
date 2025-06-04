"""
Manages the ECS specific runtime environment for the worker node.
"""

"""
Core identifiers and metadata for an ECS Fargate task
"""
struct TaskIdentifiers
    task_id::Union{String,Nothing}
    task_arn::Union{String,Nothing}
    cluster_arn::Union{String,Nothing}
    task_family::Union{String,Nothing}
    task_revision::Union{Int,Nothing}
    availability_zone::Union{String,Nothing}

    # Default constructor with all fields as nothing
    function TaskIdentifiers(
        task_id=nothing,
        task_arn=nothing,
        cluster_arn=nothing,
        task_family=nothing,
        task_revision=nothing,
        availability_zone=nothing
    )
        return new(
            task_id,
            task_arn,
            cluster_arn,
            task_family,
            task_revision,
            availability_zone
        )
    end
end

"""
Shape of the metadata response from ECS Task Metadata Endpoint V4
"""
struct TaskMetadataResponse
    TaskARN::String
    Family::String
    Revision::Int
    Cluster::String
    AvailabilityZone::String
end

"""
Retrieves identifiers and metadata for the current ECS Fargate task.
Uses the ECS Task Metadata Endpoint V4 to fetch task information.

# Returns
- `TaskIdentifiers`: Object containing task metadata and identifiers

# Throws
- `Exception`: If metadata cannot be retrieved or parsed
"""
function get_task_metadata()::TaskIdentifiers
    # ECS Fargate runtime provides this endpoint
    metadata_uri = Base.get(ENV, "ECS_CONTAINER_METADATA_URI_V4", nothing)

    # We may not have it locally - but the web API is currently tolerant of this
    if isnothing(metadata_uri)
        throw(ErrorException("Not running in ECS environment - metadata URI not found"))
    end

    try
        # Get the result from the endpoint
        response = HTTP.get("$(metadata_uri)/task")
        if response.status != 200
            throw(ErrorException("HTTP error! status: $(response.status)"))
        end

        # Cast stream -> string and then parse 
        task_metadata::TaskMetadataResponse = JSON3.read(
            String(response.body), TaskMetadataResponse
        )

        # Extract task ID from ARN
        task_arn_parts = split(task_metadata.TaskARN, '/')
        task_id = length(task_arn_parts) > 0 ? task_arn_parts[end] : "unknown"

        return TaskIdentifiers(
            task_id,
            task_metadata.TaskARN,
            task_metadata.Cluster,
            task_metadata.Family,
            task_metadata.Revision,
            task_metadata.AvailabilityZone
        )
    catch e
        if e isa HTTP.ExceptionRequest.StatusError
            throw(ErrorException("Failed to retrieve ECS task metadata: $(e.status)"))
        else
            throw(ErrorException("Failed to retrieve ECS task metadata"))
        end
    end
end

"""
Safely attempts to get task metadata with fallback values.
Won't throw errors - returns empty object if metadata can't be retrieved.

# Returns
- `TaskIdentifiers`: Metadata object with potentially missing fields
"""
function get_task_metadata_safe()::TaskIdentifiers
    try
        return get_task_metadata()
    catch e
        @warn "Failed to get task metadata. Falling back to fake values."
        return TaskIdentifiers()
    end
end
