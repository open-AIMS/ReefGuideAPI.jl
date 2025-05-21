"""Methods to prep and cache server data"""

"""
    _prep_criteria_ranges(reg_assess_data::OrderedDict)::Nothing

Determine min/max values for each criteria layer and store as entry in data store.
"""
function _prep_criteria_ranges(reg_assess_data::OrderedDict)::Nothing
    regions = get_regions()

    # Collate available criteria across all regions
    criteria_names = []
    for r in regions
        criteria_names = append!(criteria_names, names(reg_assess_data[regions[1]].stack))
    end

    unique!(criteria_names)

    # Ignore the valid data layer (not a "true" criteria layer)
    criteria_names = [c for c in criteria_names if c != :ValidSlopes]

    # Create entries for the min and max
    criteria_ranges = DataFrame([c => Number[NaN, NaN] for c in criteria_names]...)

    # Populate dataframe with data from first region
    for cn in criteria_names
        try
            c_values = reg_assess_data[regions[1]].valid_slopes[:, cn]
            criteria_ranges[:, cn] .= extrema(c_values)
        catch err
            # If error is an ArgumentError, then the criteria name was not found
            # for this specific region and can be safely skipped.
            if !(err isa ArgumentError)
                rethrow(err)
            end
        end
    end

    # Check for the extrema in other regions
    for r in regions[2:end], cn in criteria_names
        reg_criteria = names(reg_assess_data[r].valid_slopes)
        if cn ∉ reg_criteria
            # If the criteria column does not exist in `valid_slopes` (data not
            # available/applicable for region) then it can be safely skipped.
            continue
        end

        c_values = reg_assess_data[r].valid_slopes[:, cn]
        min_max = extrema(c_values)

        # Populate with the min/max
        criteria_ranges[1, cn] .= min(criteria_ranges[1, cn], min_max[1])
        criteria_ranges[2, cn] .= max(criteria_ranges[2, cn], min_max[2])
    end

    reg_assess_data["criteria_ranges"] = criteria_ranges

    return nothing
end

"""
    initialize_regional_data_cache(reef_data_path::String, reg_cache_fn::String)

Create initial regional data store with data from `reef_data_path`, excluding geospatial
data and save to `reg_cache_fn` path.
"""
function initialize_regional_data_cache(reef_data_path::String, reg_cache_fn::String)
    regional_assessment_data = OrderedDict{
        String,Union{RegionalCriteria,DataFrame,Dict}
    }()
    for reg in get_regions()
        @debug "$(now()) : Initializing cache for $reg"
        data_paths = String[]
        data_names = String[]

        slope_table = nothing
        flat_table = nothing

        for (k, dp) in criteria_data_map()
            g = glob("$reg*$dp.tif", reef_data_path)
            if length(g) == 0
                continue
            end

            push!(data_paths, first(g))
            push!(data_names, string(k))
            if occursin("valid", string(dp))
                # Load up Parquet files
                parq_file = replace(first(g), ".tif" => "_lookup.parq")

                if occursin("slope", string(dp))
                    slope_table = GeoParquet.read(parq_file)
                elseif occursin("flat", string(dp))
                    @warn "Skipping data for reef flats as it is currently unused"
                    # flat_table = GeoParquet.read(parq_file)
                else
                    msg = "Unknown lookup found: $(parq_file). Must be 'slope' or 'flat'"
                    throw(ArgumentError(msg))
                end
            end
        end

        # Pre-extract long/lat coordinates
        coords = GI.coordinates.(slope_table.geometry)
        # For each entry, extract
        slope_table[!, :lons] .= first.(coords)
        slope_table[!, :lats] .= last.(coords)

        # coords = GI.coordinates.(flat_table.geometry)
        # flat_table[!, :lons] .= first.(coords)
        # flat_table[!, :lats] .= last.(coords)

        rst_stack = RasterStack(data_paths; name=data_names, lazy=true)
        regional_assessment_data[reg] = RegionalCriteria(
            rst_stack,
            slope_table,
            slope_table[[1], :]  # Dummy entry for flat_table
        )

        @debug "$(now()) : Finished initialization for $reg"
    end

    regional_assessment_data["region_long_names"] = Dict(
        "FarNorthern" => "Far Northern Management Area",
        "Cairns-Cooktown" => "Cairns/Cooktown Management Area",
        "Townsville-Whitsunday" => "Townsville/Whitsunday Management Area",
        "Mackay-Capricorn" => "Mackay/Capricorn Management Area"
    )

    _prep_criteria_ranges(regional_assessment_data)

    # Store cache on disk to avoid excessive cold startup times
    @debug "Saving regional data cache to disk"
    serialize(reg_cache_fn, regional_assessment_data)

    return regional_assessment_data
end

"""
    setup_regional_data(config::Dict)

Load regional data to act as an in-memory cache.

# Arguments
- `config` : Configuration settings, typically loaded from a TOML file.
- `reef_data_path` : Path to pre-prepared reef data

# Returns
OrderedDict of `RegionalCriteria` for each region.
"""
function setup_regional_data(config::Dict)
    reef_data_path = config["prepped_data"]["PREPPED_DATA_DIR"]
    reg_cache_dir = config["server_config"]["REGIONAL_CACHE_DIR"]
    reg_cache_fn = joinpath(reg_cache_dir, "regional_cache.dat")

    if @isdefined(REGIONAL_DATA)
        @debug "Using previously generated regional data store."
    elseif isfile(reg_cache_fn)
        @debug "Loading regional data cache from disk"
        # Updates to packages like DiskArrays can break deserialization
        try
            @eval const REGIONAL_DATA = deserialize($(reg_cache_fn))
        catch err
            @warn "Failed to deserialize $(reg_cache_fn) with error:" err
            rm(reg_cache_fn)
        end
    end

    if !@isdefined(REGIONAL_DATA)
        @debug "Setting up regional data store..."
        regional_assessment_data = initialize_regional_data_cache(
            reef_data_path,
            reg_cache_fn
        )
        # Remember, `@eval` runs in global scope.
        @eval const REGIONAL_DATA = $(regional_assessment_data)
    end

    # If REGIONAL_DATA is defined, but failed to load supporting data that cannot be
    # cached to disk, such as the reef outlines, (e.g., due to incorrect config), then it
    # will cause errors later on.
    # Then there's no way to address this, even between web server sessions, as `const`
    # values cannot be modified.
    # Here, we check for existence and try to load again if needed.
    if !haskey(REGIONAL_DATA, "reef_outlines")
        reef_outline_path = joinpath(reef_data_path, "rrap_canonical_outlines.gpkg")
        REGIONAL_DATA["reef_outlines"] = GDF.read(reef_outline_path)
    end

    return REGIONAL_DATA
end

"""
    warmup_cache(config_path::String)

Invokes warm up of regional data cache to reduce later spin up times.
"""
function warmup_cache(config_path::String)
    config = TOML.parsefile(config_path)

    # Create re-usable empty tile
    no_data_path = cache_filename(Dict("no_data" => "none"), config, "no_data", "png")
    if !isfile(no_data_path)
        save(no_data_path, zeros(RGBA, tile_size(config)))
    end

    return setup_regional_data(config)
end
