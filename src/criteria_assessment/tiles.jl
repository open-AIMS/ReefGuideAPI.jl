"""
Helper methods to support tiling
"""

using Images, ImageIO, Interpolations

# HTTP response headers for tile images
const TILE_HEADERS = [
    "Cache-Control" => "max-age=86400, no-transform"
]

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

# Helper functions for Web Mercator projection
_lat_to_y(lat::Float64) = log(tan(π / 4 + lat * π / 360))
_y_to_lat(y::Float64) = 360 / π * atan(exp(y)) - 90

"""
    adjusted_nearest(rst::Raster, z::Int, x::Int, y::Int, tile_size::Tuple{Int,Int}, orig_rst_size::Tuple{Int,Int})::Matrix

Resample a raster using nearest neighbor interpolation when the tile includes area outside
where data exists (e.g., viewing the globe where the data may appear in a small corner of
the tile). This approach attempts to account for planetary curvature while still maintaining
some performance.

# Arguments
- `rst`: The input raster to be resampled.
- `z`: Tile zoom level requested.
- `x`: x coordinate for requested tile.
- `y`: y coordinate for the requested tile.
- `tile_size`: The desired dimensions of the tile (long, lat).

# Returns
Matrix with the resampled data.
"""
function adjusted_nearest(
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

    # Create an empty tile (long/lat)
    long_size, lat_size = tile_size
    tile = zeros(lat_size, long_size)

    # Generate longitude and latitude arrays for the tile
    lons = @. mod(
        t_lon_min + (t_lon_max - t_lon_min) * ((1:long_size) - 1) / (long_size - 1) + 180,
        360
    ) - 180
    lats = @. _y_to_lat(
        _lat_to_y(t_lat_max) -
        (_lat_to_y(t_lat_max) - _lat_to_y(t_lat_min)) * ((1:lat_size) - 1) / (lat_size - 1)
    )

    # Determine which points are within the area of interest
    in_lons = aoi_lon_min .<= lons .<= aoi_lon_max
    in_lats = aoi_lat_min .<= lats .<= aoi_lat_max

    # Sample data that is within area of interest
    for (i, lon) in enumerate(lons), (j, lat) in enumerate(lats)
        if in_lons[i] && in_lats[j]
            x_idx =
                round(
                    Int,
                    (lon - aoi_lon_min) / (aoi_lon_max - aoi_lon_min) * (size(rst, 1) - 1)
                ) + 1
            y_idx =
                round(
                    Int,
                    (_lat_to_y(aoi_lat_max) - _lat_to_y(lat)) /
                    (_lat_to_y(aoi_lat_max) - _lat_to_y(aoi_lat_min)) * (size(rst, 2) - 1)
                ) + 1
            x_idx = clamp(x_idx, 1, size(rst, 1))
            y_idx = clamp(y_idx, 1, size(rst, 2))
            tile[j, i] = rst[x_idx, y_idx]
        end
    end

    return tile
end

function setup_tile_routes(config, auth)
    @get auth("/to-tile/{zoom}/{lon}/{lat}") function (
        req::Request, zoom::Int64, lon::Float64, lat::Float64
    )
        x, y = _lon_lat_to_tile(zoom, lon, lat)
        return json(Dict(:x => x, :y => y))
    end

    @get auth("/to-lonlat/{zoom}/{x}/{y}") function (
        req::Request, zoom::Int64, x::Int64, y::Int64
    )
        lon_min, lon_max, lat_max, lat_min = _tile_bounds(zoom, x, y)
        return json(
            Dict(
                :lon_min => lon_min,
                :lon_max => lon_max,
                :lat_max => lat_max,
                :lat_min => lat_min
            )
        )
    end

    reg_assess_data = setup_regional_data(config)
    @get auth("/tile/{z}/{x}/{y}") function (req::Request, z::Int64, x::Int64, y::Int64)
        # http://127.0.0.1:8000/tile/{z}/{x}/{y}?region=Cairns-Cooktown&rtype=slopes&Depth=-9.0:0.0&Slope=0.0:40.0&Rugosity=0.0:3.0
        # http://127.0.0.1:8000/tile/8/231/139?region=Cairns-Cooktown&rtype=slopes&Depth=-9.0:0.0&Slope=0.0:40.0&Rugosity=0.0:3.0
        # http://127.0.0.1:8000/tile/7/115/69?region=Cairns-Cooktown&rtype=slopes&Depth=-9.0:0.0&Slope=0.0:40.0&Rugosity=0.0:3.0
        # http://127.0.0.1:8000/tile/8/231/139?region=Cairns-Cooktown&rtype=slopes&Depth=-9.0:0.0&Slope=0.0:40.0&Rugosity=0.0:3.0

        qp = queryparams(req)
        mask_path = cache_filename(qp, config, "", "png")
        if isfile(mask_path)
            return file(mask_path; headers=TILE_HEADERS)
        end

        # Otherwise, create the file
        thread_id = Threads.threadid()
        @debug "Thread $(thread_id) - $(now()) : Assessing criteria"
        # Filtering time: 0.6 - 7.0 seconds
        reg = qp["region"]
        rtype = qp["rtype"]

        criteria_names, lbs, ubs = remove_rugosity(reg, parse_criteria_query(qp)...)

        # Calculate tile bounds
        lon_min, lon_max, lat_max, lat_min = _tile_bounds(z, x, y)
        @debug "Thread $(thread_id) - $(now()) : Calculated bounds (z/x/y, lon bounds, lat bounds): $z $x $y | $(_tile_to_lon_lat(z, x, y)) | ($(lon_min), $(lon_max)), ($(lat_min), $(lat_max))"

        # Extract relevant data based on tile coordinates
        @debug "Thread $(thread_id) - $(now()) : Extracting tile data"
        mask_data = threshold_mask(
            reg_assess_data[reg],
            Symbol(rtype),
            CriteriaBounds.(criteria_names, lbs, ubs),
            (lon_min, lon_max),
            (lat_min, lat_max)
        )

        if any(size(mask_data) .== 0)
            no_data_path = cache_filename(
                Dict("no_data" => "none"), config, "no_data", "png"
            )

            @debug "Thread $(thread_id) - No data for $reg ($rtype) at $z/$x/$y"
            return file(no_data_path; headers=TILE_HEADERS)
        end

        @debug "Thread $(thread_id) - Extracted data size: $(size(mask_data))"

        @debug "Thread $(thread_id) - $(now()) : Creating PNG (with transparency)"
        img = zeros(RGBA, tile_size(config))
        if (z < 12)
            # Account for geographic positioning when zoomed out further than
            # raster area
            resampled = adjusted_nearest(mask_data, z, x, y, tile_size(config))
        else
            # Zoomed in close so less need to account for curvature
            # BSpline(Constant()) is equivalent to nearest neighbor.
            # See details in: https://juliaimages.org/ImageTransformations.jl/stable/reference/#Low-level-warping-API
            resampled = imresize(
                mask_data.data', tile_size(config); method=BSpline(Constant())
            )
        end

        img[resampled .== 1] .= RGBA(0, 0, 0, 1)

        @debug "Thread $(thread_id) - $(now()) : Saving and serving file"
        save(mask_path, img)
        return file(mask_path; headers=TILE_HEADERS)
    end
end
