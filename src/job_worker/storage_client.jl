"""
An abstract storage client and implementation for S3 - which is untested/not
used yet.
"""

"""
Abstract type for storage clients
All concrete storage clients should inherit from this
"""
abstract type StorageClient end

"""
Parse a storage URI into its components
Returns a tuple of (scheme, bucket, path)
"""
function parse_storage_uri(uri::String)::Tuple{String,String,String}
    # Match s3://bucket-name/path or other schemes
    m = match(r"^([a-z0-9]+)://([^/]+)/(.*)$", uri)
    if isnothing(m)
        throw(ArgumentError("Invalid storage URI format: $uri"))
    end

    scheme = m.captures[1]
    bucket = m.captures[2]
    path = m.captures[3]

    return (scheme, bucket, path)
end

"""
S3 Storage Client implementation
"""
struct S3StorageClient <: StorageClient
    region::String
    s3_endpoint::OptionalValue{String}

    # Constructor with defaults
    function S3StorageClient(;
        region::String,
        s3_endpoint::OptionalValue{String}=nothing
    )
        return new(region, s3_endpoint)
    end
end

"""
Upload a file to the specified storage URI
"""
function upload_file(
    client::S3StorageClient,
    local_path::String,
    storage_path::String
)::String
    try
        # Parse the storage URI - concat the target path + base path
        scheme, bucket, key = parse_storage_uri(storage_path)

        if scheme != "s3"
            throw(ArgumentError("Expected S3 URI, got $scheme"))
        end

        @debug "Uploading file from $(local_path) to $(storage_path)"

        aws =
            !isnothing(client.s3_endpoint) ?
            # Use minio special config if we need to
            Minio.MinioConfig(
                client.s3_endpoint;
                # TODO would be nicer to pass through config tree
                username=ENV["MINIO_USERNAME"],
                password=ENV["MINIO_PASSWORD"]
            ) :
            # Otherwise use typical AWS config
            AWS.AWSConfig(; region=client.region)

        # Read the file content
        file_data = Base.read(local_path)

        # Upload to S3
        AWSS3.s3_put(aws, bucket, key, file_data)

        # For demonstration purposes only - simulate the upload
        @debug "Uploaded file to $storage_path"

        return storage_path
    catch e
        @error "Failed to upload file to S3: $e" exception = (e, catch_backtrace())
        rethrow(e)
    end
end

"""
Factory function to create an S3 storage client based on environment variables
"""
function create_s3_client()::S3StorageClient
    # Get configuration from environment
    region = Base.get(ENV, "AWS_REGION", "ap-southeast-2")
    return S3StorageClient(;
        region=region
    )
end
