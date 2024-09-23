using JSONWebTokens
using HTTP
using Dates
using JSON

const CORS_HEADERS = [
    "Access-Control-Allow-Origin" => "*",
    "Access-Control-Allow-Headers" => "*",
    "Access-Control-Allow-Methods" => "POST, GET, OPTIONS"
]

# https://juliaweb.github.io/HTTP.jl/stable/examples/#Cors-Server
function CorsMiddleware(handler)
    return function(req::HTTP.Request)
        @debug "CORS middleware"
        # determine if this is a pre-flight request from the browser
        if HTTP.method(req)=="OPTIONS"
            return HTTP.Response(200, CORS_HEADERS)  
        else 
            return handler(req) # passes the request to the AuthMiddleware
        end
    end
end


function setup_jwt_middleware(config::Dict)
    if !get(config["jwt_auth"], "JWT_ENABLED", false)
        return identity  # Return a pass-through middleware if JWT is not enabled
    end

    jwt_iss = config["jwt_auth"]["JWT_ISS"]
    # WKT endpoint for rest API
    public_key_url = config["jwt_auth"]["WKT_ENDPOINT"]

    # Fetch the public key
    response = HTTP.get(public_key_url)
    jwks_json = JSON.parse(String(response.body))
    # Assuming the first key in the JWKS is the one we want
    public_key = jwks_json["keys"][1]["n"]
    rsa_public = JSONWebTokens.RS256(public_key)

    function jwt_auth_middleware(handler)
        return function(req::HTTP.Request)
            auth_header = HTTP.header(req, "Authorization", "")
            if !startswith(auth_header, "Bearer ")
                return HTTP.Response(401, "Unauthorized: Missing or invalid Authorization header.")
            end

            token = strip(auth_header[8:end])  # Remove "Bearer " prefix

            try
                payload = JSONWebTokens.decode(rsa_public, token)
                
                # Check if the token is expired
                exp = get(payload, "exp", nothing)
                if exp !== nothing && exp < time()
                    return HTTP.Response(401, "Unauthorized: Token expired")
                end

                # Check if the issuer matches
                if get(payload, "iss", "") != jwt_iss
                    return HTTP.Response(401, "Unauthorized: Invalid token issuer")
                end

                # If we've made it this far, the token is valid
                return handler(req)
            catch e
                if e isa JSONWebTokens.InvalidSignatureError
                    return HTTP.Response(401, "Unauthorized: Invalid token signature")
                elseif e isa JSONWebTokens.MalformedJWTError
                    return HTTP.Response(401, "Unauthorized: Invalid token format")
                else
                    # Log the error for debugging
                    @error "Unexpected error during JWT validation" exception=(e, catch_backtrace())
                    return HTTP.Response(500, "Internal Server Error")
                end
            end
        end
    end

    return jwt_auth_middleware
end