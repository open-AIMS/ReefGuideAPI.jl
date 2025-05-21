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

    function AssessmentCriteria(;
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

# Helps to find the slopes lookup file - note this is a bit fragile
const SLOPES_LOOKUP_SUFFIX = "_lookup.parq"
function get_slope_parquet_filepath(region::RegionDetails)::String
    return "$(region.id)$(ASSESSMENT_CRITERIA.slope.file_suffix)$(SLOPES_LOOKUP_SUFFIX)"
end

function populate_data_driven_criteria(region::RegionDetails)::RegionalCriteria
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

mutable struct RegionalDataEntry
    "The unique ID of the region"
    region_id::String
    "Other metadata about the region"
    region_metadata::RegionDetails
    "The rasters associated with this region"
    raster_stack::Rasters.RasterStack
    "A DataFrame containing coordinates of valid slope reefs"
    slope_coordinates::DataFrame
    "Information about the criteria for this region"
    criteria_ranges::AssessmentCriteria

    function RegionalDataEntry(;
        region_id::String,
        region_metadata::RegionDetails,
        raster_stack::Rasters.RasterStack,
        slope_coordinates::DataFrame,
        criteria_ranges::Vector{RegionalCriteria}
    )
        return new(
            region_id, region_metadata, raster_stack, slope_coordinates, criteria_ranges
        )
    end
end

const RegionalDataType = Dict{String,RegionalDataEntry}
REGIONAL_DATA::OptionalValue{RegionalDataType} = nothing

function check_existing_regional_data_from_memory()::OptionalValue{RegionalDataType}
    if !isnothing(REGIONAL_DATA)
        @debug "Using previously generated regional data store."
        return REGIONAL_DATA
    end
    return nothing
end

function check_existing_regional_data_from_disk(
    cache_directory::String
)::OptionalValue{RegionalDataType}
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
For each row of the GeoDataFrame, adds a lons/lats entry which is a vector of
coordinates
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

function initialise_data(config::Dict)::RegionalDataType
    regional_data::RegionalDataType = Dict()
    data_source_directory = config["prepped_data"]["PREPPED_DATA_DIR"]
    # iterate through regions
    for region_details::RegionDetails in REGIONS
        @debug "$(now()) : Initializing cache for $(region_details)"

        # Setup data paths and names arrays
        data_paths = String[]
        data_names = String[]

        slope_table = nothing

        # TODO others
        assessable_criteria = [
            ASSESSMENT_CRITERIA.depth
        ]

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

        # Extract the slope table for this region
        slope_table::DataFrame = GeoParquet.read(
            get_slope_parquet_filepath(region_details)
        )

        # Extract the lat longs
        add_lat_long_columns_to_dataframe(slope_table)

        # Build the raster stack 
        raster_stack = RasterStack(data_paths; name=data_names, lazy=true)

        regional_data[region_details.id] = RegionalDataEntry(;
            region_id=region_details.id,
            region_metadata=region_details,
            raster_stack=raster_stack,
            slope_coordinates=slope_table,
            criteria_ranges=nothing)
    end
end

"""
"""
function initialise_data_with_cache(config::Dict; force_cache_invalidation::Bool=false)
    # Use module scope variable here
    global REGIONAL_DATA

    # Where is cache located?
    regional_cache_directory = config["server_config"]["REGIONAL_CACHE_DIR"]

    # do we have an in-memory cache available? 
    local_data = check_existing_regional_data_from_memory()
    if !isnothing(local_data)
        REGIONAL_DATA = local_data
        return nothing
    end

    # do we have disk cache available? 
    disk_data = check_existing_regional_data_from_disk(regional_cache_directory)
    if !isnothing(disk_data)
        REGIONAL_DATA = disk_data
        return nothing
    end

    # No cache exists - recompute

end
