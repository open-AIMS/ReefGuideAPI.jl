"""
"""

"""
    within_thresholds(ctype::Val, data::Raster, lb::T, ub::T) where {T}

Apply in-bound constraints.

# Notes
Why is this a simple one line function?
Because we want to be able to cache results for each constraint type.
"""
function within_thresholds(reg::Val, ctype::Val, data::Raster, lb::T, ub::T) where {T}
    return within_thresholds(data, lb, ub)
end
function within_thresholds(req, data::N, lb::T, ub::T) where {N,T}
    return within_thresholds(data, lb, ub)
end
@memoize function within_thresholds(data::Raster, lb::T, ub::T) where {T}
    return (lb .<= data .<= ub)
end
@memoize function within_thresholds(data::Vector, lb::T, ub::T) where {T}
    return (lb .<= data .<= ub)
end

function within_thresholds(data, lb::T, ub::T) where {T}
    return (lb .<= data .<= ub)
end

"""
    port_buffer_mask(gdf::DataFrame, dist::Float64; unit::String="NM")

Create a masking buffer around indicated port locations.

# Arguments
- `gdf` : GeoDataFrame of port locations (given as long/lat points)
- `dist` : distance from port in degrees (deg), kilometers (km), or nautical miles (NM; default)
- `unit` : unit `dist` is in
"""
function port_buffer_mask(gdf::DataFrame, dist::Float64; unit::String="NM")
    # Determine conversion factor (nautical miles or kilometers)
    conv_factor = 1.0
    if unit == "NM"
        conv_factor = 60.0  # 60 NM = 1 degree
    elseif unit == "km"
        conv_factor = 111.0  # 111 km = 1 degree
    elseif unit != "deg"
        error("Unknown distance unit requested. Can only be one of `NM` or `km` or `deg`")
    end

    ports = gdf.geometry  # TODO: Replace with `GI.geometrycolumns()`

    # Make buffer around ports
    buffered_ports = GO.buffer.(ports, dist / conv_factor)

    # Combine all geoms into one
    port_mask = reduce((x1, x2) -> LibGEOS.union(x1, x2), buffered_ports)

    return port_mask
end

"""
    filter_distances(
        target_rast::Raster,
        dist_buffer
    )::Raster

Apply a mask to exclude pixels that are outside the indicated distance buffer(s).

`target_rast` and the `dist_buffer` should be in the same CRS (e.g., EPSG:7844 / GDA2020).

# Arguments
- `target_rast` : Raster of suitable pixels (Bool) to filter pixels from.
- `dist_buffer` : Buffer geometry to use as the mask.

# Returns
- Masked boolean raster indicating pixels that are within the target distance.
"""
function filter_distances(target_rast::Raster, dist_buffer)::Raster
    # Mask out areas outside considered distance from port
    return mask(Raster(target_rast; missingval=0); with=dist_buffer)
end

"""
    valid_lonlat_inds(data::DataFrame, criteria::Symbol, lb::T, ub::T) where {T}

Retrieve the indices of valid data for a region.

# Arguments
- `data` :
- `criteria` :
- `lb` :
- `ub` :

# Returns
Tuple{Vector{Int64}, Vector{Int64}}, of lon and lat indices.
"""
function valid_lonlat_inds(data::DataFrame, criteria::Symbol, lb::T, ub::T) where {T}
    valid_locs = within_thresholds(data[!, criteria], lb, ub)

    lon_pos = data[valid_locs, :lon_idx]
    lat_pos = data[valid_locs, :lat_idx]

    return lon_pos, lat_pos
end

"""
    valid_pixel_positions(data::DataFrame, criteria::Symbol, lb::T, ub::T) where {T}

Obtain the pixel positions of valid data.

Intended for use in applications similar to [ImageryLayer - client side pixel filter](https://developers.arcgis.com/javascript/latest/sample-code/layers-imagery-pixelvalues/).

# Arguments
- `data` :
- `criteria` :
- `lb` : lower bound
- `ub` : upper bound
"""
function valid_pixel_positions(data::DataFrame, criteria::Symbol, lb::T, ub::T) where {T}
    lon_pos, lat_pos = valid_lonlat_inds(data, criteria, lb, ub)
    pixel_pos = lon_pos .* lat_pos

    return pixel_pos
end

function _create_filter(bounds::Tuple)
    return (x) -> bounds[1] .< x .<= bounds[2]
end

"""
    apply_criteria_thresholds(criteria_stack::RasterStack, lookup::DataFrame, ruleset::Vector{CriteriaBounds{Function}})::Raster
    apply_criteria_thresholds(criteria_stack::RasterStack, lookup::DataFrame, ruleset::Dict)::Raster
    apply_criteria_thresholds(criteria_stack::RasterStack, lookup::DataFrame, ruleset::NamedTuple)::Raster

Apply thresholds for each criteria.

# Arguments
- `criteria_stack` : RasterStack of criteria data for a given region
- `lookup` : Lookup dataframe for the region
- `ruleset` : A set of CriteriaBounds, Dictionary or NamedTuple indicating a mapping of
              criteria names to their lower and upper thresholds.

# Returns
BitMatrix indicating locations within desired thresholds
"""
function apply_criteria_thresholds(
    criteria_stack::RasterStack,
    lookup::DataFrame,
    ruleset::Dict
)::Raster
    ruleset = NamedTuple{(keys(ruleset)...,)}(
        Tuple(_create_filter.(values(ruleset)))
    )

    return apply_criteria_thresholds(criteria_stack, lookup, ruleset)
end
function apply_criteria_thresholds(
    criteria_stack::RasterStack,
    lookup::DataFrame,
    ruleset::NamedTuple
)::Raster
    # Result store
    res = Raster(criteria_stack.Depth; data=falses(size(criteria_stack)), missingval=0)

    res_lookup = trues(nrow(lookup))
    for rule_name in keys(ruleset)
        res_lookup .= res_lookup .& ruleset[rule_name](lookup[!, rule_name])
    end

    tmp = lookup[res_lookup, [:lon_idx, :lat_idx]]
    res[CartesianIndex.(tmp.lon_idx, tmp.lat_idx)] .= true

    return res
end
function apply_criteria_thresholds(
    criteria_stack::T,
    lookup::DataFrame,
    ruleset::Vector{CriteriaBounds{Function}}
)::Raster where {T}
    # Result store
    res = Raster(criteria_stack.Depth; data=falses(size(criteria_stack)), missingval=0)

    res_lookup = trues(nrow(lookup))
    for threshold in ruleset
        res_lookup .= res_lookup .& threshold.rule(lookup[!, threshold.name])
    end

    tmp = lookup[res_lookup, [:lon_idx, :lat_idx]]
    res[CartesianIndex.(tmp.lon_idx, tmp.lat_idx)] .= true

    return res
end

"""
    apply_criteria_lookup(
        reg_criteria,
        rtype::Symbol,
        ruleset::Vector{CriteriaBounds{Function}}
    )

Filter valid_lookup table by applying user defined `ruleset` criteria.

# Arguments
- `reg_criteria` : RegionalCriteria containing valid_rtype lookup table for filtering.
- `rtype` : Flats or slope category for assessment.
- `ruleset` : User defined ruleset for upper and lower bounds.

# Returns
Filtered lookup table containing points that meet all criteria in `ruleset`.
"""
function apply_criteria_lookup(
    reg_criteria::RegionalCriteria,
    rtype::Symbol,
    ruleset::Vector{CriteriaBounds{Function}}
)
    lookup = getfield(reg_criteria, Symbol(:valid_, rtype))
    lookup.all_crit .= 1

    for threshold in ruleset
        lookup.all_crit = lookup.all_crit .& threshold.rule(lookup.name)
    end

    lookup = lookup[lookup.all_crit, :]
    lookup.lon = first.(GI.coordinates.(lookup.geometry))
    lookup.lat = last.(GI.coordinates.(lookup.geometry))

    return lookup
end

"""
    make_threshold_mask(reg::String, rtype::Symbol, crit_map)

Generate mask for a given region and reef type (slopes or flats) according to thresholds
applied to a set of criteria.

# Notes
- Zeros indicate locations to mask **out**.
- Ones indicate locations to **keep**.

# Arguments
- `reg_criteria` : RegionalCriteria to assess
- `rtype` : reef type to assess (`:slopes` or `:flats`)
- `crit_map` : List of criteria thresholds to apply (see `apply_criteria_thresholds()`)
- `lons` : Longitudinal extent (min and max, required when generating masks for tiles)
- `lats` : Latitudinal extent (min and max, required when generating masks for tiles)

# Returns
True/false mask indicating locations within desired thresholds.
"""
function make_threshold_mask(reg_criteria, rtype::Symbol, crit_map)::Raster
    valid_lookup = getfield(reg_criteria, Symbol(:valid_, rtype))
    mask_layer = apply_criteria_thresholds(
        reg_criteria.stack,
        valid_lookup,
        crit_map
    )

    return mask_layer
end
function make_threshold_mask(
    reg_criteria,
    rtype::Symbol,
    crit_map,
    lons::Tuple,
    lats::Tuple
)::Raster
    lookup = getfield(reg_criteria, Symbol(:valid_, rtype))

    (lat1, lat2) = lats[1] > lats[2] ? (lats[2], lats[1]) : (lats[1], lats[2])

    within_search = (
        (lons[1] .<= lookup.lons .<= lons[2]) .&
        (lat1 .<= lookup.lats .<= lat2)
    )
    lookup = lookup[within_search, :]

    # Need to pass in full representation of the raster as the lookup table relies on
    # the original Cartesian indices.
    res = apply_criteria_thresholds(
        reg_criteria.stack,
        lookup,
        crit_map
    )

    # Extract data between lon/lats
    view_of_data = view(res, X(lons[1] .. lons[2]), Y(lat1 .. lat2))
    return rebuild(view_of_data, sparse(convert.(UInt8, view_of_data)))
end

"""
    generate_criteria_mask!(fn::String, rst_stack::RasterStack, lookup::DataFrame, ruleset::Vector{CriteriaBounds{Function}})

Generate mask file for a given region and reef type (slopes or flats) according to thresholds
applied to a set of criteria.

# Notes
- Zeros indicate locations to mask **out**.
- Ones indicate locations to **keep**.

# Arguments
- `fn` : File to write geotiff to
- `reg_criteria` : RegionalCriteria to assess
- `rtype` : reef type to assess (`:slopes` or `:flats`)
- `crit_map` : List of criteria thresholds to apply (see `apply_criteria_thresholds()`)

# Returns
Nothing
"""
function generate_criteria_mask!(
    fn::String,
    rst_stack::RasterStack,
    lookup::DataFrame,
    ruleset::Vector{CriteriaBounds{Function}}
)::Nothing
    # Create the geotiff
    res = spzeros(size(rst_stack))
    tmp_rst = Raster(
        rst_stack[names(rst_stack)[1]];
        data=res,
        missingval=0.0
    )

    res_lookup = trues(nrow(lookup))
    for threshold in ruleset
        res_lookup .= res_lookup .& threshold.rule(lookup[!, threshold.name])
    end

    tmp = lookup[res_lookup, [:lon_idx, :lat_idx]]
    tmp_rst[CartesianIndex.(tmp.lon_idx, tmp.lat_idx)] .= 1.0

    write(
        fn,
        UInt8.(tmp_rst);
        ext=".tiff",
        source="gdal",
        driver="COG",  # GTiff
        options=Dict{String,String}(
            "COMPRESS" => "LZW",
            "SPARSE_OK" => "TRUE",
            "OVERVIEW_COUNT" => "5"
        )
    )

    return nothing
end
