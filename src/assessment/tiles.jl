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
    fast_resample(rst::Raster, new_size::Tuple{Int, Int})::Matrix

Quickly resample a Raster to a new size using nearest neighbor interpolation.
Applies a simple nearest neighbor algorithm, prioritising performance over accuracy.

# Arguments
- `rst::Raster`: The input raster to be resampled.
- `new_size::Tuple{Int, Int}`: The desired dimensions of the output raster as (width, height).

# Returns
Matrix with the resampled data.

# Examples
```julia
large_raster = Raster(rand(UInt8, 14756, 14838); dims=(X(1:1:14756), Y(1:1:14838)))
small_matrix = fast_resample(large_raster, (256, 256))
```
"""
function fast_resample(rst::Raster, new_size::Tuple{Int, Int})::Matrix
    old_size = size(rst)
    x_ratio = old_size[1] / new_size[1]
    y_ratio = old_size[2] / new_size[2]

    resampled = Array{eltype(rst)}(undef, new_size)

    Threads.@threads for y in 1:new_size[2]
        for x in 1:new_size[1]
            src_x = round(Int, x * x_ratio)
            src_y = round(Int, y * y_ratio)
            resampled[x, y] = rst[src_x, src_y]
        end
    end

    return resampled
end
