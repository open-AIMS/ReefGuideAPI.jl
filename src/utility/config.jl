"""Methods for accessing properties of config."""

"""
    _cache_location(config::Dict)::String

Retrieve cache location for geotiffs.
"""
function _cache_location(config::Dict)::String
    cache_loc = try
        in_debug = haskey(config["server_config"], "DEBUG_MODE")
        if in_debug && lowercase(config["server_config"]["DEBUG_MODE"]) == "true"
            if "DEBUG_CACHE_DIR" ∉ keys(ENV)
                ENV["DEBUG_CACHE_DIR"] = mktempdir()
            end

            ENV["DEBUG_CACHE_DIR"]
        else
            config["server_config"]["TIFF_CACHE_DIR"]
        end
    catch err
        @warn "Encountered error:" err
        if "DEBUG_CACHE_DIR" ∉ keys(ENV)
            ENV["DEBUG_CACHE_DIR"] = mktempdir()
        end

        ENV["DEBUG_CACHE_DIR"]
    end

    return cache_loc
end

"""
    n_gdal_threads(config::Dict)::String

Retrieve the configured number of threads to use when writing COGs with GDAL.
"""
function _n_gdal_threads(config::Dict)::String
    n_cog_threads = try
        config["server_config"]["COG_THREADS"]
    catch
        "1"  # Default to using a single thread for GDAL write
    end

    return n_cog_threads
end
