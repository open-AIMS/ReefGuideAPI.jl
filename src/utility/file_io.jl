"""Methods for common file I/O."""

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
            "NUM_THREADS" => _n_gdal_threads(config)
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
