"""
A HTTP client which wraps auth/headers to talk to the web API for the job
system.
"""

"""
Custom error type for API errors
"""
struct ApiError <: Exception
    message::String
    status_code::Int
    response::Any

    ApiError(message::String, status_code::Integer, response=nothing) = new(
        message, status_code, response
    )
end

"""
Represents JWT payload structure 

(This is what's in the token)
"""
struct JWTPayload
    id::String
    email::String
    roles::Vector{String}
    exp::Int
end

"""
Credentials for authentication
"""
struct Credentials
    email::String
    password::String
end

"""
Authentication tokens
"""
struct AuthTokens
    token::String
    refresh_token::Union{String,Nothing}
end

"""
Client for authenticated API requests
"""
mutable struct AuthApiClient
    base_url::String
    credentials::Credentials
    tokens::Union{AuthTokens,Nothing}
    http_headers::Dict{String,String}
    token_refresh_threshold_ms::Int

    function AuthApiClient(base_url::String, credentials::Credentials)
        return new(
            base_url,
            credentials,
            nothing,
            Dict("Content-Type" => "application/json"),
            60  # 1 minute in seconds
        )
    end
end

"""
Get a valid token, refreshing if necessary
"""
function get_valid_token(client::AuthApiClient)::Union{String,Nothing}
    if isnothing(client.tokens)
        @debug "Logging in"
        login!(client)
        return isnothing(client.tokens) ? nothing : client.tokens.token
    end

    # Decode token to check expiration (no need to validate so we specify no encoding)
    decoded_token = JSONWebTokens.decode(JSONWebTokens.None(), client.tokens.token)

    exp_time = decoded_token["exp"]

    # Get current unix timestamp 
    current_unix = floor(Int64, Dates.datetime2unix(Dates.now(Dates.UTC)))

    # convert to ms
    expires_in = exp_time - current_unix

    # Refresh if close to expiration
    @debug "Token expires in $(expires_in) seconds..."
    if expires_in <= client.token_refresh_threshold_ms
        @debug "Refreshing now"
        refresh_token!(client)
    end

    return isnothing(client.tokens) ? nothing : client.tokens.token
end

"""
Login to get new tokens
"""
function login!(client::AuthApiClient)
    try
        response = HTTP.post(
            "$(client.base_url)/auth/login",
            client.http_headers,
            JSON3.write(
                Dict(
                    "email" => client.credentials.email,
                    "password" => client.credentials.password
                )
            )
        )

        if response.status == 200
            response_data = JSON3.read(String(response.body))
            client.tokens = AuthTokens(
                response_data.token,
                haskey(response_data, :refreshToken) ? response_data.refreshToken : nothing
            )
        else
            throw(ApiError("Login failed", response.status, String(response.body)))
        end
    catch e
        if e isa HTTP.ExceptionRequest.StatusError
            throw(ApiError("Failed to login", e.status, String(e.response.body)))
        else
            throw(ApiError("Failed to login", 500, nothing))
        end
    end
end

"""
Refresh the authentication token
"""
function refresh_token!(client::AuthApiClient)
    @debug "Token refresh started at: " Dates.format(now(), "yyyy-mm-dd HH:MM:SS")

    try
        # If no refresh token available, try logging in
        if isnothing(client.tokens) || isnothing(client.tokens.refresh_token)
            @debug "No refresh token, logging in..."
            login!(client)
            return nothing
        end

        response = HTTP.post(
            "$(client.base_url)/auth/token",
            client.http_headers,
            JSON3.write(Dict(
                "refreshToken" => client.tokens.refresh_token
            ))
        )

        if response.status != 200
            @info "Non 200 response from refresh token endpoint: " response.status
            throw(ApiError("Non 200 response from refresh token", response.status))
        end

        response_data = JSON3.read(String(response.body))
        client.tokens = AuthTokens(
            response_data.token,
            client.tokens.refresh_token
        )
    catch e
        @info "Error caught during refresh"
        # If refresh fails, try logging in again
        client.tokens = nothing
        login!(client)
    end

    @debug "Token refresh completed at: " Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
end

"""
Get authorization headers with token
"""
function auth_headers(client::AuthApiClient)
    headers = copy(client.http_headers)
    @debug "Getting valid token"
    token = get_valid_token(client)
    if !isnothing(token)
        headers["Authorization"] = "Bearer $token"
    end
    return headers
end

"""
Determine if path requires authentication
"""
function requires_auth(path::String)::Bool
    return !(
        endswith(path, "/auth/login") ||
        endswith(path, "/auth/register") ||
        endswith(path, "/auth/token")
    )
end

# HTTP Methods

function setupHeaders(client::AuthApiClient, path::String)
    return (
        url="$(client.base_url)$path",
        headers=requires_auth(path) ? auth_headers(client) : client.http_headers
    )
end

function HTTPGet(client::AuthApiClient, path::String; params::Dict=Dict())::Any
    url, headers = setupHeaders(client, path)

    response = HTTP.get(url, headers; query=params)
    return JSON3.read(String(response.body))
end

function HTTPPost(
    client::AuthApiClient, path::String, data::Union{Dict,Nothing}=nothing
)::Union{Any,Nothing}
    url, headers = setupHeaders(client, path)

    body = isnothing(data) ? "" : JSON3.write(data)
    response = HTTP.post(url, headers, body)
    response_content = String(response.body)
    @debug "Response content = $(response_content)"

    if isempty(response_content)
        return nothing
    else
        return JSON3.read(response_content)
    end
end

function HTTPPut(
    client::AuthApiClient, path::String, data::Union{Dict,Nothing}=nothing
)::Any
    url, headers = setupHeaders(client, path)

    body = isnothing(data) ? "" : JSON3.write(data)
    response = HTTP.put(url, headers, body)

    if isempty(String(response.body))
        return nothing
    else
        return JSON3.read(String(response.body))
    end
end

function HTTPPatch(
    client::AuthApiClient, path::String, data::Union{Dict,Nothing}=nothing
)::Any
    url, headers = setupHeaders(client, path)

    body = isnothing(data) ? "" : JSON3.write(data)
    response = HTTP.patch(url, headers, body)

    if isempty(String(response.body))
        return nothing
    else
        return JSON3.read(String(response.body))
    end
end

function HTTPDelete(client::AuthApiClient, path::String)::Any
    url, headers = setupHeaders(client, path)

    response = HTTP.delete(url, headers)

    if isempty(String(response.body))
        return nothing
    else
        return JSON3.read(String(response.body))
    end
end
