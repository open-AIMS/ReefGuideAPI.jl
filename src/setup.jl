using Rasters

const REGIONAL_DATA_CACHE_FILENAME = "regional_cache_v2.dat"

@enum Region begin
    TOWNSVILLE_WHITSUNDAY
    CAIRNS_COOKTOWN
    MACKAY_CAPRICORN
    FAR_NORTHERN
end

struct RegionDetails
    region::Region
    display_name::String
    id::String

    function RegionDetails(;
        region::Region,
        display_name::String,
        id::String
    )
        return new(region, display_name, id)
    end
end

const REGIONS::Vector{RegionDetails} = [
    RegionDetails(;
        region=Region.TOWNSVILLE_WHITSUNDAY,
        display_name="Townsville/Whitsunday Management Area",
        id="Townsville-Whitsunday"
    ),
    RegionDetails(;
        region=Region.CAIRNS_COOKTOWN,
        display_name="Cairns/Cooktown Management Area",
        id="Cairns-Cooktown"
    ),
    RegionDetails(;
        region=Region.MACKAY_CAPRICORN,
        display_name="Mackay/Capricorn Management Area",
        id="Mackay-Capricorn"
    ),
    RegionDetails(;
        region=Region.FAR_NORTHERN,
        display_name="Far Northern Management Area",
        id="Far-Northern"
    )
]

"""
Simple struct to contain min/max values
"""
mutable struct Bounds
    min::Float32
    max::Float32

    function Bounds(; min::Number, max::Number)
        return new(parse(Float32, min), parse(Float32, max))
    end
end

function bounds_from_tuple(min_max::Tuple{Number,Number})::Bounds
    return Bounds(; min=min_max[1], max=min_max[2])
end

struct CriteriaMetadata
    "System ID used to uniquely identify this criteria"
    id::String
    "What is the suffix of this criteria in the data?"
    file_suffix::String
    "Pretty display name, can be changed with a guarantee of no runtime consequences/errors"
    display_label::String

    function CriteriaMetadata(;
        id::String,
        file_suffix::String,
        display_label::String
    )
        return new(
            id,
            file_suffix,
            display_label
        )
    end
end

struct AssessmentCriteria
    depth::CriteriaMetadata
    slope::CriteriaMetadata
    turbidity::CriteriaMetadata
    waves_height::CriteriaMetadata
    waves_period::CriteriaMetadata
    rugosity::CriteriaMetadata

    function AssessmentCriteria(;
        depth::CriteriaMetadata,
        slope::CriteriaMetadata,
        turbidity::CriteriaMetadata,
        waves_height::CriteriaMetadata,
        waves_period::CriteriaMetadata,
        rugosity::CriteriaMetadata
    )
        return new(
            depth, slope, turbidity, waves_height, waves_period, rugosity
        )
    end
end

const ASSESSMENT_CRITERIA::AssessmentCriteria = AssessmentCriteria(;
    depth=CriteriaMetadata(;
        id="Depth",
        file_suffix="_bathy",
        display_label="Depth"
    ),
    slope=CriteriaMetadata(;
        id="Slope",
        file_suffix="_slope",
        display_label="Slope"
    ),
    turbidity=CriteriaMetadata(;
        id="Turbidity",
        file_suffix="_turbid",
        display_label="Turbidity"
    ),
    waves_height=CriteriaMetadata(;
        id="WavesHs",
        file_suffix="_waves_Hs",
        display_label="Wave Height (m)"
    ),
    waves_period=CriteriaMetadata(;
        id="WavesTp",
        file_suffix="_waves_Tp",
        display_label="Wave Period (s)"
    ),
    rugosity=CriteriaMetadata(;
        id="Rugosity",
        file_suffix="_rugosity",
        display_label="Rugosity"
    )
)

struct RegionalCriteriaEntry
    "Other details about the criteria"
    metadata::CriteriaMetadata
    "Bounds (min/max)"
    bounds::CriteriaBounds

    function RegionCriteriaEntry(;
        metadata::CriteriaMetadata,
        bounds::Bounds
    )
        return new(
            metadata,
            bounds
        )
    end
end

struct RegionalCriteria
    depth::RegionalCriteriaEntry
    slope::RegionalCriteriaEntry
    turbidity::RegionalCriteriaEntry
    waves_height::RegionalCriteriaEntry
    waves_period::RegionalCriteriaEntry
    rugosity::RegionalCriteriaEntry

    function RegionalCriteria(;
        depth_bounds::Bounds,
        slope_bounds::Bounds,
        turbidity_bounds::Bounds,
        waves_height_bounds::Bounds,
        waves_period_bounds::Bounds,
        rugosity_bounds::Bounds
    )
        return new(
            RegionalCriteriaEntry(;
                metadata=ASSESSMENT_CRITERIA.depth, bounds=depth_bounds
            ),
            RegionalCriteriaEntry(;
                metadata=ASSESSMENT_CRITERIA.slope, bounds=slope_bounds
            ),
            RegionalCriteriaEntry(;
                metadata=ASSESSMENT_CRITERIA.turbidity, bounds=turbidity_bounds
            ),
            RegionalCriteriaEntry(;
                metadata=ASSESSMENT_CRITERIA.waves_height, bounds=waves_height_bounds
            ),
            RegionalCriteriaEntry(;
                metadata=ASSESSMENT_CRITERIA.waves_period, bounds=waves_period_bounds
            ),
            RegionalCriteriaEntry(;
                metadata=ASSESSMENT_CRITERIA.rugosity, bounds=rugosity_bounds
            )
        )
    end
end

"""
Takes a slope_table data frame and builds out the region specific bounds for
each criteria defined as the extrema of present values over the region
"""
function build_assessment_criteria_from_slope_table(table::DataFrame)::RegionalCriteria
    const depth_bounds = bounds_from_tuple(extrema(table[:, ASSESSMENT_CRITERIA.depth.id]))
    const slope_bounds = bounds_from_tuple(extrema(table[:, ASSESSMENT_CRITERIA.slope.id]))
    const turbidity_bounds = bounds_from_tuple(
        extrema(table[:, ASSESSMENT_CRITERIA.turbidity.id])
    )
    const waves_height_bounds = bounds_from_tuple(
        extrema(table[:, ASSESSMENT_CRITERIA.waves_height.id])
    )
    const waves_period_bounds = bounds_from_tuple(
        extrema(table[:, ASSESSMENT_CRITERIA.waves_period.id])
    )
    const rugosity_bounds = bounds_from_tuple(
        extrema(table[:, ASSESSMENT_CRITERIA.rugosity.id])
    )

    return RegionalCriteria(;
        depth_bounds,
        slope_bounds,
        turbidity_bounds,
        waves_height_bounds,
        waves_period_bounds,
        rugosity_bounds
    )
end

# Helps to find the slopes lookup file - note this is a bit fragile
const SLOPES_LOOKUP_SUFFIX = "_lookup.parq"
function get_slope_parquet_filepath(region::RegionDetails)::String
    return "$(region.id)$(ASSESSMENT_CRITERIA.slope.file_suffix)$(SLOPES_LOOKUP_SUFFIX)"
end

mutable struct RegionalDataEntry
    "The unique ID of the region"
    region_id::String
    "Other metadata about the region"
    region_metadata::RegionDetails
    "The rasters associated with this region"
    raster_stack::Rasters.RasterStack
    "A DataFrame containing coordinates of valid slope reefs"
    slope_table::DataFrame
    "Information about the criteria for this region"
    criteria::RegionalCriteria
    "The canonical reef outlines geo data frame"
    reef_outlines::DataFrame

    function RegionalDataEntry(;
        region_id::String,
        region_metadata::RegionDetails,
        raster_stack::Rasters.RasterStack,
        slope_table::DataFrame,
        criteria::RegionalCriteria,
        reef_outlines::DataFrame
    )
        return new(
            region_id, region_metadata, raster_stack, slope_table, criteria, reef_outlines
        )
    end
end

const RegionalDataMapType = Dict{String,RegionalDataEntry}
struct RegionalData
    "Dictionary mapping the region ID (see region details) to the data entry"
    regions::RegionalDataMapType
    "Canonical reef outlines"
    reef_outlines::DataFrame

    function RegionalData(;
        regions::Dict{String,RegionalDataEntry},
        reef_outlines::DataFrame
    )
        return new(regions, reef_outlines)
    end
end

REGIONAL_DATA::OptionalValue{RegionalData} = nothing

function check_existing_regional_data_from_memory()::OptionalValue{RegionalData}
    if !isnothing(REGIONAL_DATA)
        @debug "Using previously generated regional data store."
        return REGIONAL_DATA
    end
    return nothing
end

function check_existing_regional_data_from_disk(
    cache_directory::String
)::OptionalValue{RegionalData}
    # check if we have an existing cache file
    reg_cache_filename = joinpath(cache_directory, REGIONAL_DATA_CACHE_FILENAME)
    if isfile(reg_cache_filename)
        @debug "Loading regional data cache from disk"
        try
            return deserialize(reg_cache_filename)
        catch err
            @warn "Failed to deserialize $(reg_cache_filename) with error:" err
            # Invalidating cache
            rm(reg_cache_filename)
        end
    end
    # no success
    return nothing
end

function find_data_source_for_criteria(;
    data_source_directory::String, region::RegionDetails, criteria::CriteriaMetadata
)::String
    # ascertain the file pattern 
    matched_files = glob("$(region.id)*$(criteria.file_suffix).tif", data_source_directory)

    # ensure there is exactly one matching file
    if length(matched_files) == 0
        throw(ErrorException("Did not find data for the criteria: $(criteria.id)."))
    end
    if length(matched_files) > 1
        throw(
            ErrorException(
                "Found more than one data source match for criteria: $(criteria.id). This is ambiguous, unsure how to proceed."
            )
        )
    end

    return first(matched_files)
end

"""
For each row of the GeoDataFrame, adds a lons/lats entry which is the coordinate
of the centroid of the geometry of that reef
"""
function add_lat_long_columns_to_dataframe(df::DataFrame)::Nothing
    # This returns a Vector{Tuple{number, number}} i.e. list of lat/lons
    coords = GI.coordinates.(df.geometry)
    # This adds a new column where each value is the first entry from the tuple (i.e. lon)
    df[!, :lons] .= first.(coords)
    # This adds a new column where each value is the second entry from the tuple (i.e. lat)
    df[!, :lats] .= last.(coords)
    return nothing
end

const DEFAULT_CANONICAL_REEFS_FILE_NAME = "rrap_canonical_outlines.gpkg"
function load_canonical_reefs(
    source_dir::String; file_name::String=DEFAULT_CANONICAL_REEFS_FILE_NAME
)::DataFrame
    # Load in the reef canonical outlines geodataframe
    return GDF.read(
        joinpath(source_dir, file_name)
    )
end

function initialise_data(config::Dict)::RegionalData
    regional_data::RegionalDataMapType = Dict()
    data_source_directory = config["prepped_data"]["PREPPED_DATA_DIR"]
    # iterate through regions
    for region_metadata::RegionDetails in REGIONS
        @debug "$(now()) : Initializing cache for $(region_metadata)"

        # Setup data paths and names arrays
        data_paths = String[]
        data_names = String[]

        slope_table = nothing

        # list out the criteria we are interested in pulling out raster and
        # criteria ranges for 
        assessable_criteria = [
            ASSESSMENT_CRITERIA.depth,
            ASSESSMENT_CRITERIA.rugosity,
            ASSESSMENT_CRITERIA.slope,
            ASSESSMENT_CRITERIA.turbidity,
            ASSESSMENT_CRITERIA.waves_height,
            ASSESSMENT_CRITERIA.waves_period
        ]

        # Extract the slope table for this region
        slope_table::DataFrame = GeoParquet.read(
            get_slope_parquet_filepath(region_metadata)
        )

        # Extract the lat longs
        add_lat_long_columns_to_dataframe(slope_table)

        for criteria::CriteriaMetadata in assessable_criteria
            # Get the file name of the matching data layer (.tif) and add it to
            # the raster stack
            push!(
                data_paths,
                find_data_source_for_criteria(;
                    data_source_directory,
                    region,
                    criteria
                )
            )
            # Add name to - this is helpful so you can reference the raster in
            # the raster stack
            push!(data_names, criteria.id)
        end

        # Determine the criteria and extrema 
        criteria::AssessmentCriteria = build_assessment_criteria_from_slope_table(
            slope_table
        )

        # Build the raster stack 
        raster_stack = RasterStack(data_paths; name=data_names, lazy=true)

        # Add entry to the regional data dictionary
        regional_data[region_metadata.id] = RegionalDataEntry(;
            region_id=region_metadata.id,
            region_metadata,
            raster_stack,
            slope_table,
            criteria)
    end

    # Load canonical reefs
    canonical_reefs = load_canonical_reefs(data_source_directory)

    # All done - return populated data
    return RegionalData(; regions=regional_data, reef_outlines=canonical_reefs)
end

function get_empty_tile_path(config::Dict)::String
    cache_location = _cache_location(config)
    return joinpath(cache_location, "no_data_tile.png")
end

function setup_empty_tile_cache(config::Dict)::String
    file_path = get_empty_tile_path(config)
    if !isfile(file_path)
        save(file_path, zeros(RGBA, tile_size(config)))
    end
end

function initialise_data_with_cache(config::Dict; force_cache_invalidation::Bool=false)
    # Use module scope variable here
    global REGIONAL_DATA

    # Where is cache located?
    regional_cache_directory = config["server_config"]["REGIONAL_CACHE_DIR"]

    if !force_cache_invalidation
        # do we have an in-memory cache available? 
        local_data = check_existing_regional_data_from_memory()
        if !isnothing(local_data) && !force_cache_invalidation
            REGIONAL_DATA = local_data
            return nothing
        end

        # do we have disk cache available? 
        disk_data = check_existing_regional_data_from_disk(regional_cache_directory)
        if !isnothing(disk_data) && !force_cache_invalidation
            REGIONAL_DATA = disk_data
            return nothing
        end
    end

    # No cache exists OR we are forcing invalidation
    regional_data = initialise_data(config)

    # Update the global variable
    REGIONAL_DATA = regional_data

    # Also ensure our empty tile cache is ready
    setup_empty_tile_cache(config)

    return nothing
end
