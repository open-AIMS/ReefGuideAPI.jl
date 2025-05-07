"""Methods for common file I/O."""

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
    cache_filename(qp::Dict, config::Dict, suffix::String, ext::String)

Generate a filename for a cache.

# Arguments
- `qp` : Query parameters to hash
- `config` : app configuration (to extract cache parent directory from)
- `suffix` : a suffix to use in the filename (pass `""` if none required)
- `ext` : file extension to use
"""
function cache_filename(qp::Dict, config::Dict, suffix::String, ext::String)
    file_id = create_job_id(qp)
    temp_path = _cache_location(config)
    cache_file_path = joinpath(temp_path, "$(file_id)$(suffix).$(ext)")

    return cache_file_path
end

"""
    n_gdal_threads(config::Dict)::String

Retrieve the configured number of threads to use when writing COGs with GDAL.
"""
function n_gdal_threads(config::Dict)::String
    n_cog_threads = try
        config["server_config"]["COG_THREADS"]
    catch
        "1"  # Default to using a single thread for GDAL write
    end

    return n_cog_threads
end

"""
    _write_cog(file_path::String, data::Raster, config::Dict)::Nothing

Write out a COG using common options.

# Arguments
- `file_path` : Path to write data out to
- `data` : Raster data to write out
"""
function _write_cog(file_path::String, data::Raster, config::Dict)::Nothing
    Rasters.write(
        file_path,
        data;
        ext=".tiff",
        source="gdal",
        driver="COG",
        options=Dict{String,String}(
            "COMPRESS" => "DEFLATE",
            "SPARSE_OK" => "TRUE",
            "OVERVIEW_COUNT" => "5",
            "BLOCKSIZE" => string(first(tile_size(config))),
            "NUM_THREADS" => n_gdal_threads(config)
        ),
        force=true
    )

    return nothing
end

"""
    _write_tiff(file_path::String, data::Raster)::Nothing

Write out a geotiff using common options.

# Arguments
- `file_path` : Path to write data out to
- `data` : Raster data to write out
"""
function _write_tiff(file_path::String, data::Raster)::Nothing
    Rasters.write(
        file_path,
        data;
        ext=".tiff",
        source="gdal",
        driver="gtiff",
        force=true
    )

    return nothing
end
