"""Methods for identifying potential deployment locations."""

using FLoops, ThreadsX

"""
    proportion_suitable(subsection::BitMatrix)::Matrix{Int16}

Calculate the the proportion of the subsection that is suitable for deployments.
Subsection is the surrounding a rough hectare area centred on each cell of a raster marked
as being suitable according to user-selected criteria.
"""
function proportion_suitable(x::BitMatrix)::Matrix{Int16}
    x_len, y_len = size(x)
    x′ = zeros(Int16, size(x))

    @floop for row_col in ThreadsX.findall(x)
        (row, col) = Tuple(row_col)
        x_left = max(col - 4, 1)
        x_right = min(col + 4, x_len)

        y_top = max(row - 4, 1)
        y_bottom = min(row + 4, y_len)

        x′[row, col] = Int16(sum(@views x[y_top:y_bottom, x_left:x_right]))
    end

    return x′
end

"""
    filter_distances(
        target_rast::Raster,
        gdf::DataFrame,
        dist_nm
    )::Raster

Exclude pixels in target_rast that are beyond `dist_nm` (nautical miles) from a geometry
in `gdf`. Target_rast and gdf should be in the same CRS (EPSG:7844 / GDA2020 for GBR-reef-guidance-assessment).

# Arguments
- `target_rast` : Raster of suitable pixels (Bool) to filter pixels from.
- `gdf` : DataFrame with `geometry` column that contains vector objects of interest.
- `dist_nm` : Filtering distance from geometry object in nautical miles.

# Returns
- `tmp_areas` : Raster of filtered pixels containing only pixels within target distance
from a geometry centroid.
"""
function filter_distances(target_rast::Raster, gdf::DataFrame, dist; units::String="NM")::Raster
    tmp_areas = copy(target_rast)

    # First dimension is the rows (latitude)
    # Second dimension is the cols (longitude)
    raster_lat = Vector{Float64}(tmp_areas.dims[1].val)
    raster_lon = Vector{Float64}(tmp_areas.dims[2].val)

    @floop for row_col in findall(tmp_areas)
        (lat_ind, lon_ind) = Tuple(row_col)
        point = AG.createpoint()

        lon = raster_lon[lon_ind]
        lat = raster_lat[lat_ind]
        AG.addpoint!(point, lon, lat)

        pixel_dists = AG.distance.([point], port_locs.geometry)
        geom_point = gdf[argmin(pixel_dists), :geometry]
        geom_point = (AG.getx(geom_point, 0), AG.gety(geom_point, 0))

        dist_nearest = Distances.haversine(geom_point, (lon, lat))

        # Convert from meters to nautical miles
        if units == "NM"
            dist_nearest = dist_nearest / 1852
        end

        # Convert from meters to kilometers
        if units == "km"
            dist_nearest = dist_nearest / 1000
        end

        tmp_areas.data[lon_ind, lat_ind] = dist_nearest < dist ? 1 : 0
    end

    return tmp_areas
end

"""
# TODO: Better name, and address duplication in "/assess/{reg}/{rtype}"
"""
function _temp_assess_region(reg_assess_data, reg, qp, rtype, config)
    criteria_names, lbs, ubs = remove_rugosity(reg, parse_criteria_query(qp)...)

    # Otherwise, create the file
    @debug "$(now()) : Assessing criteria"
    assess = reg_assess_data[reg]
    mask_data = make_threshold_mask(
        assess,
        Symbol(rtype),
        CriteriaBounds.(criteria_names, lbs, ubs)
    )

    return mask_data
end

function assess_region(reg_assess_data, reg, qp, rtype, config)
    @debug "Assessing region's suitability score"

    file_id = string(hash(qp))
    assessed_tmp_path = _cache_location(config)
    assessed_path = joinpath(assessed_tmp_path, file_id*"_suitable.tiff")

    if isfile(assessed_path)
        return file(assessed_path)
    end

    # Make mask of suitable locations
    mask_data = _temp_assess_region(reg_assess_data, reg, qp, rtype, config)

    # Assess remaining pixels for their suitability
    @debug "Calculating proportional suitability score"
    suitability_scores = proportion_suitable(mask_data.data)

    @debug "$(now()) : Running on thread $(threadid())"
    @debug "Writing to $(assessed_path)"
    Rasters.write(
        assessed_path,
        rebuild(mask_data, suitability_scores);
        ext=".tiff",
        source="gdal",
        driver="COG",
        options=Dict{String,String}(
            "COMPRESS"=>"DEFLATE",
            "SPARSE_OK"=>"TRUE",
            "OVERVIEW_COUNT"=>"5",
            "BLOCKSIZE"=>"256",
            "NUM_THREADS"=>n_gdal_threads(config)
        ),
        force=true
    )

    return file(assessed_path)


    # Apply rotation-based polygon search
    # assess_reef_site()

    # Filter overlapping polygons.

    # Return geojson of suitable deployment "plots"
    # output_geojson()
end
