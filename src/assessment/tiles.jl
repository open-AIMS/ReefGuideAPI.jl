"""
Helper methods to support tiling
"""

"""
    _tile_to_lon_lat(z::T, x::T, y::T) where {T<:Int64}

Obtain lon/lat of top-left corner of a requested tile.

# Returns
lon, lat
"""
function _tile_to_lon_lat(z::T, x::T, y::T) where {T<:Int64}
    n = 2.0^z
    lon = x / n * 360.0 - 180.0
    lat_rad = atan(sinh(π * (1 - 2 * y / n)))
    lat = rad2deg(lat_rad)

    return (lon, lat)
end

"""
    _tile_bounds(z::T, x::T, y::T) where {T<:Int64}

Obtain lon/lat bounds of a requested tile.

# Returns
West, East, North South (min lon, max lon, lat max, lat min)
"""
function _tile_bounds(z::T, x::T, y::T) where {T<:Int64}
    # Calculate the boundaries of the tile
    n = 2.0^z
    lat_min = atan(sinh(π * (1 - 2 * (y + 1) / n))) * 180.0 / π
    lat_max = atan(sinh(π * (1 - 2 * y / n))) * 180.0 / π
    lon_min = x / n * 360.0 - 180.0
    lon_max = (x + 1) / n * 360.0 - 180.0

    # West, East, North, South
    return lon_min, lon_max, lat_max, lat_min
end

"""
    _lon_lat_to_tile(zoom, lon, lat)

Identify the corresponding tile coordinates for a given lon/lat.

# Returns
x and y tile coordinates
"""
function _lon_lat_to_tile(zoom, lon, lat)
    n = 2.0^zoom
    x = floor(Int64, (lon + 180.0) / 360.0 * n)

    lat_rad = lat * π / 180.0
    y = floor(Int64, (1.0 - log(tan(lat_rad) + 1.0 / cos(lat_rad)) / π) / 2.0 * n)

    return x, y
end

"""
    nearest(rst::Raster, tile_size::Tuple{Int, Int})::Matrix

Resample a raster to a tile size using nearest neighbor interpolation.
This approach prioritising performance over accuracy.

# Arguments
- `rst`: The input raster to be resampled.
- `tile_size`: The desired dimensions of the tile (lat, long).

# Returns
Matrix with the resampled data.

# Examples
```julia
large_raster = Raster(rand(UInt8, 14756, 14838); dims=(X(1:1:14756), Y(1:1:14838)))
small_matrix = nearest(large_raster, (256, 256))
```
"""
function nearest(rst::Raster, tile_size::Tuple{Int, Int})::Matrix
    old_size = size(rst)

    # Important: must flip axes!
    # Rasters.jl stores longitude along first dimension (rows) by default.
    x_ratio = old_size[1] / tile_size[2]
    y_ratio = old_size[2] / tile_size[1]

    resampled = zeros(eltype(rst), tile_size)

    Threads.@threads for lat in 1:tile_size[1]
        for lon in 1:tile_size[2]
            # Use area averaging for downsampling
            x_start = max(1, floor(Int, (lon - 1) * x_ratio) + 1)
            x_end = min(old_size[1], ceil(Int, lon * x_ratio))
            y_start = max(1, floor(Int, (lat - 1) * y_ratio) + 1)
            y_end = min(old_size[2], ceil(Int, lat * y_ratio))

            count_val = count(rst[x_start:x_end, y_start:y_end].data .> 0)
            if count_val == 0
                continue
            end

            sum_val = sum(rst[x_start:x_end, y_start:y_end].data)
            resampled[lat, lon] = ceil(sum_val / count_val)
        end
    end

    return resampled
end

"""
    masked_nearest(rst::Raster, z::Int, x::Int, y::Int, tile_size::Tuple{Int,Int}, orig_rst_size::Tuple{Int,Int})::Matrix

Resample a raster using nearest neighbor interpolation when the tile includes area outside
where data exists (e.g., viewing the globe where the data may appear in a small corner of
the tile). This approach prioritising performance over accuracy.

# Arguments
- `rst`: The input raster to be resampled.
- `z`: Tile zoom level requested.
- `x`: x coordinate for requested tile.
- `y`: y coordinate for the requested tile.
- `tile_size`: The desired dimensions of the tile (lat, long).

# Returns
Matrix with the resampled data.
"""
function masked_nearest(
    rst::Raster,
    z::Int,
    x::Int,
    y::Int,
    tile_size::Tuple{Int,Int}
)::Matrix
    # Bounds for the requested tile
    (t_lon_min, t_lon_max, t_lat_max, t_lat_min) = _tile_bounds(z, x, y)

    # Bounds for the area of interest (AOI; where we have data)
    ((aoi_lon_min, aoi_lon_max), (aoi_lat_min, aoi_lat_max)) = Rasters.bounds(rst)

    # Create an empty tile (lat/long)
    lat_size = tile_size[1]
    long_size = tile_size[2]
    tile = fill(0.0, lat_size, long_size)

    lons = @. t_lon_min + (t_lon_max - t_lon_min) * ((1:long_size) - 1) / (long_size - 1)
    lats = @. t_lat_max - (t_lat_max - t_lat_min) * ((1:lat_size) - 1) / (lat_size - 1)
    in_lons = aoi_lon_min .<= lons .<= aoi_lon_max
    in_lats = aoi_lat_min .<= lats .<= aoi_lat_max

    # Sample data that is within area of interest
    data_x = round.(Int, (lons[in_lons] .- aoi_lon_min) / (aoi_lon_max - aoi_lon_min) * size(rst, 1))
    data_y = round.(Int, (aoi_lat_max .- lats[in_lats]) / (aoi_lat_max - aoi_lat_min) * size(rst, 2))

    tile[in_lats, in_lons] .= rst[data_x, data_y]'

    return tile
end

function setup_tile_routes(auth)
    @get("/to-tile/{zoom}/{lon}/{lat}", middleware=[auth]) function (req::Request, zoom::Int64, lon::Float64, lat::Float64)
        x, y = _lon_lat_to_tile(zoom, lon, lat)
        return json(Dict(:x=>x, :y=>y))
    end

    @get("/to-lonlat/{zoom}/{x}/{y}", middleware=[auth]) function (req::Request, zoom::Int64, x::Int64, y::Int64)
        lon_min, lon_max, lat_max, lat_min = _tile_bounds(zoom, x, y)
        return json(
                Dict(
                :lon_min=>lon_min,
                :lon_max=>lon_max,
                :lat_max=>lat_max,
                :lat_min=>lat_min
            )
        )
    end

end
