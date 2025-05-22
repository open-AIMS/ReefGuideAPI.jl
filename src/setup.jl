using Rasters
using Glob
using GeoParquet
using Serialization
using Oxygen: json, Request
using FileIO

# =============================================================================
# Constants and Configuration
# =============================================================================

const REGIONAL_DATA_CACHE_FILENAME = "regional_cache_v2.dat"
const SLOPES_LOOKUP_SUFFIX = "_lookup.parq"
const DEFAULT_CANONICAL_REEFS_FILE_NAME = "rrap_canonical_outlines.gpkg"

# Global variable to store regional data cache
REGIONAL_DATA::OptionalValue{RegionalData} = nothing

# =============================================================================
# Core Data Structures
# =============================================================================

"""
Metadata container for regional information.

# Fields
- `display_name::String` : Human-readable name for UI display
- `id::String` : Unique system identifier (changing this affects data loading)
"""
struct RegionMetadata
    display_name::String
    id::String

    function RegionMetadata(;
        display_name::String,
        id::String
    )
        return new(display_name, id)
    end
end

"""
Simple container for minimum and maximum boundary values.

# Fields
- `min::Float32` : Minimum value
- `max::Float32` : Maximum value
"""
mutable struct Bounds
    min::Float32
    max::Float32

    function Bounds(; min::Number, max::Number)
        return new(Float32(min), Float32(max))
    end
end

"""
Metadata for assessment criteria including file naming conventions.

# Fields
- `id::String` : Unique system identifier for the criteria
- `file_suffix::String` : File suffix pattern for data files
- `display_label::String` : Human-readable label for UI display
"""
struct CriteriaMetadata
    id::String
    file_suffix::String
    display_label::String

    function CriteriaMetadata(;
        id::String,
        file_suffix::String,
        display_label::String
    )
        return new(id, file_suffix, display_label)
    end
end

"""
Container for all assessment criteria metadata.

# Fields
- `depth::CriteriaMetadata` : Bathymetry/depth criteria
- `slope::CriteriaMetadata` : Slope gradient criteria  
- `turbidity::CriteriaMetadata` : Water turbidity criteria
- `waves_height::CriteriaMetadata` : Wave height criteria
- `waves_period::CriteriaMetadata` : Wave period criteria
- `rugosity::CriteriaMetadata` : Seafloor rugosity criteria
"""
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

"""
Combines criteria metadata with regional boundary values.

# Fields
- `metadata::CriteriaMetadata` : Criteria definition and metadata
- `bounds::Bounds` : Min/max values for this criteria in the region
"""
struct RegionalCriteriaEntry
    metadata::CriteriaMetadata
    bounds::Bounds

    function RegionalCriteriaEntry(;
        metadata::CriteriaMetadata,
        bounds::Bounds
    )
        return new(metadata, bounds)
    end
end

"""
Complete set of regional criteria with computed bounds for each parameter.

# Fields
- `depth::RegionalCriteriaEntry` : Depth criteria and bounds
- `slope::RegionalCriteriaEntry` : Slope criteria and bounds
- `turbidity::RegionalCriteriaEntry` : Turbidity criteria and bounds
- `waves_height::RegionalCriteriaEntry` : Wave height criteria and bounds
- `waves_period::RegionalCriteriaEntry` : Wave period criteria and bounds
- `rugosity::RegionalCriteriaEntry` : Rugosity criteria and bounds
"""
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
Complete data package for a single region including rasters and metadata.

# Fields
- `region_id::String` : Unique identifier for the region
- `region_metadata::RegionMetadata` : Display metadata for the region
- `raster_stack::Rasters.RasterStack` : Geospatial raster data layers
- `slope_table::DataFrame` : Coordinates and values for valid slope reef locations
- `criteria::RegionalCriteria` : Computed criteria bounds for this region
"""
mutable struct RegionalDataEntry
    region_id::String
    region_metadata::RegionMetadata
    raster_stack::Rasters.RasterStack
    slope_table::DataFrame
    criteria::RegionalCriteria

    function RegionalDataEntry(;
        region_id::String,
        region_metadata::RegionMetadata,
        raster_stack::Rasters.RasterStack,
        slope_table::DataFrame,
        criteria::RegionalCriteria
    )
        return new(region_id, region_metadata, raster_stack, slope_table, criteria)
    end
end

"""
Top-level container for all regional data and reef outlines.

# Fields
- `regions::Dict{String,RegionalDataEntry}` : Regional data indexed by region ID
- `reef_outlines::DataFrame` : Canonical reef outline geometries
"""
struct RegionalData
    regions::Dict{String,RegionalDataEntry}
    reef_outlines::DataFrame

    function RegionalData(;
        regions::Dict{String,RegionalDataEntry},
        reef_outlines::DataFrame
    )
        return new(regions, reef_outlines)
    end
end

# =============================================================================
# Configuration Constants
# =============================================================================

# Define all available regions for the assessment system
const REGIONS::Vector{RegionMetadata} = [
    RegionMetadata(;
        display_name="Townsville/Whitsunday Management Area",
        id="Townsville-Whitsunday"
    ),
    RegionMetadata(;
        display_name="Cairns/Cooktown Management Area",
        id="Cairns-Cooktown"
    ),
    RegionMetadata(;
        display_name="Mackay/Capricorn Management Area",
        id="Mackay-Capricorn"
    ),
    RegionMetadata(;
        display_name="Far Northern Management Area",
        id="Far-Northern"
    )
]

# Define all assessment criteria with file naming conventions
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

# Convenience list for iteration over all criteria
const ASSESSMENT_CRITERIA_LIST::Vector{CriteriaMetadata} = [
    ASSESSMENT_CRITERIA.depth,
    ASSESSMENT_CRITERIA.slope,
    ASSESSMENT_CRITERIA.turbidity,
    ASSESSMENT_CRITERIA.waves_height,
    ASSESSMENT_CRITERIA.waves_period,
    ASSESSMENT_CRITERIA.rugosity
]

# Type alias for cleaner code
const RegionalDataMapType = Dict{String,RegionalDataEntry}

# =============================================================================
# Utility Functions
# =============================================================================

"""
Convert a min/max tuple to a Bounds struct.

# Arguments
- `min_max::Tuple{Number,Number}` : Tuple containing (minimum, maximum) values

# Returns
`Bounds` struct with converted float values.
"""
function bounds_from_tuple(min_max::Tuple{Number,Number})::Bounds
    return Bounds(; min=min_max[1], max=min_max[2])
end

"""
Generate the filename for slope lookup data for a given region.

# Arguments  
- `region::RegionMetadata` : Region metadata containing ID

# Returns
String filename in format "{region_id}_slope_lookup.parq"
"""
function get_slope_parquet_filename(region::RegionMetadata)::String
    return "$(region.id)$(ASSESSMENT_CRITERIA.slope.file_suffix)$(SLOPES_LOOKUP_SUFFIX)"
end

"""
Create a dictionary mapping criteria IDs to regional criteria entries.

# Arguments
- `region_data::RegionalDataEntry` : Regional data containing criteria information

# Returns
Dictionary with criteria ID strings as keys and RegionalCriteriaEntry as values.
"""
function build_regional_criteria_dictionary(
    region_data::RegionalDataEntry
)::Dict{String,RegionalCriteriaEntry}
    return Dict(
        ASSESSMENT_CRITERIA.depth.id => region_data.criteria.depth,
        ASSESSMENT_CRITERIA.slope.id => region_data.criteria.slope,
        ASSESSMENT_CRITERIA.turbidity.id => region_data.criteria.turbidity,
        ASSESSMENT_CRITERIA.waves_height.id => region_data.criteria.waves_height,
        ASSESSMENT_CRITERIA.waves_period.id => region_data.criteria.waves_period,
        ASSESSMENT_CRITERIA.rugosity.id => region_data.criteria.rugosity
    )
end

"""
Add longitude and latitude columns to a DataFrame based on geometry centroids.

Modifies the input DataFrame by adding 'lons' and 'lats' columns extracted
from the centroid coordinates of each geometry feature.

# Arguments
- `df::DataFrame` : DataFrame with geometry column containing spatial features
"""
function add_lat_long_columns_to_dataframe(df::DataFrame)::Nothing
    # Extract coordinate tuples from geometry centroids
    coords = GI.coordinates.(df.geometry)
    # Add longitude column (first coordinate)
    df[!, :lons] .= first.(coords)
    # Add latitude column (second coordinate) 
    df[!, :lats] .= last.(coords)
    return nothing
end

# =============================================================================
# Data Loading and Processing Functions
# =============================================================================

"""
Build regional criteria bounds from slope table data.

Computes min/max bounds for each assessment criteria by finding extrema
in the slope table data columns.

# Arguments
- `table::DataFrame` : Slope table containing criteria data columns

# Returns
`RegionalCriteria` struct with computed bounds for all criteria.
"""
function build_assessment_criteria_from_slope_table(table::DataFrame)::RegionalCriteria
    # Compute bounds for each criteria using column extrema
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

"""
Find the data source file for a specific criteria and region.

Uses glob pattern matching to locate the appropriate .tif file based on
region ID and criteria file suffix.

# Arguments
- `data_source_directory::String` : Directory containing data files
- `region::RegionMetadata` : Region metadata for ID matching
- `criteria::CriteriaMetadata` : Criteria metadata for suffix matching

# Returns
String path to the matching data file.

# Throws
- `ErrorException` : If zero or multiple matching files are found
"""
function find_data_source_for_criteria(;
    data_source_directory::String,
    region::RegionMetadata,
    criteria::CriteriaMetadata
)::String
    # Search for files matching the pattern: {region_id}*{criteria_suffix}.tif
    matched_files = glob("$(region.id)*$(criteria.file_suffix).tif", data_source_directory)

    # Validate exactly one match exists
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
Load canonical reef outline geometries from geopackage file.

# Arguments
- `source_dir::String` : Directory containing the reef outlines file
- `file_name::String` : Name of the geopackage file (defaults to canonical name)

# Returns
DataFrame containing reef geometry data.
"""
function load_canonical_reefs(
    source_dir::String;
    file_name::String=DEFAULT_CANONICAL_REEFS_FILE_NAME
)::DataFrame
    # Load reef outlines from geopackage format
    return GDF.read(joinpath(source_dir, file_name))
end

# =============================================================================
# Cache Management Functions
# =============================================================================

"""
Check if regional data is available in memory cache.

# Returns
`RegionalData` if available in memory, `nothing` otherwise.
"""
function check_existing_regional_data_from_memory()::OptionalValue{RegionalData}
    if !isnothing(REGIONAL_DATA)
        @debug "Using previously generated regional data store."
        return REGIONAL_DATA
    end
    return nothing
end

"""
Check if regional data cache exists on disk and attempt to load it.

# Arguments
- `cache_directory::String` : Directory where cache files are stored

# Returns
`RegionalData` if successfully loaded from disk, `nothing` otherwise.
"""
function check_existing_regional_data_from_disk(
    cache_directory::String
)::OptionalValue{RegionalData}
    # Construct cache file path
    reg_cache_filename = joinpath(cache_directory, REGIONAL_DATA_CACHE_FILENAME)

    if isfile(reg_cache_filename)
        @debug "Loading regional data cache from disk"
        try
            return deserialize(reg_cache_filename)
        catch err
            @warn "Failed to deserialize $(reg_cache_filename) with error:" err
            # Remove corrupted cache file
            rm(reg_cache_filename)
        end
    end
    # No cache available or load failed
    return nothing
end

"""
Get the file path for the empty tile cache.

# Arguments
- `config::Dict` : Configuration dictionary containing cache settings

# Returns
String path to the empty tile cache file.
"""
function get_empty_tile_path(config::Dict)::String
    cache_location = _cache_location(config)
    return joinpath(cache_location, "no_data_tile.png")
end

"""
Initialize empty tile cache if it doesn't exist.

Creates a blank PNG tile used for areas with no data coverage.

# Arguments
- `config::Dict` : Configuration dictionary containing tile and cache settings
"""
function setup_empty_tile_cache(config::Dict)::Nothing
    file_path = get_empty_tile_path(config)
    if !isfile(file_path)
        # Create empty RGBA tile with configured dimensions
        save(file_path, zeros(RGBA, tile_size(config)))
    end
    return nothing
end

# =============================================================================
# Main Data Initialization Functions
# =============================================================================

"""
Initialize all regional data from source files.

Loads raster data, slope tables, and computes criteria bounds for all regions.
This is the main data loading function that builds the complete data structure.

# Arguments
- `config::Dict` : Configuration dictionary containing data directory paths

# Returns
`RegionalData` struct containing all loaded and processed regional information.
"""
function initialise_data(config::Dict)::RegionalData
    regional_data::RegionalDataMapType = Dict()
    data_source_directory = config["prepped_data"]["PREPPED_DATA_DIR"]

    # Process each region sequentially
    for region_metadata::RegionMetadata in REGIONS
        @debug "$(now()) : Initializing cache for $(region_metadata)"

        # Initialize data collection arrays
        data_paths = String[]
        data_names = String[]

        # Load slope table containing valid reef coordinates and criteria values
        slope_table::DataFrame = GeoParquet.read(
            joinpath(data_source_directory, get_slope_parquet_filename(region_metadata))
        )

        # Add coordinate columns for spatial referencing
        add_lat_long_columns_to_dataframe(slope_table)

        # Collect raster file paths for all criteria
        for criteria::CriteriaMetadata in ASSESSMENT_CRITERIA_LIST
            # Find the corresponding .tif file for this criteria
            push!(
                data_paths,
                find_data_source_for_criteria(;
                    data_source_directory,
                    region=region_metadata,
                    criteria
                )
            )
            # Use criteria ID as the raster layer name
            push!(data_names, criteria.id)
        end

        # Compute regional criteria bounds from slope table data
        criteria::RegionalCriteria = build_assessment_criteria_from_slope_table(slope_table)

        # Create lazy-loaded raster stack from all criteria files
        raster_stack = RasterStack(data_paths; name=data_names, lazy=true)

        # Store complete regional data entry
        regional_data[region_metadata.id] = RegionalDataEntry(;
            region_id=region_metadata.id,
            region_metadata,
            raster_stack,
            slope_table,
            criteria
        )
    end

    # Load canonical reef outlines that apply to all regions
    canonical_reefs = load_canonical_reefs(data_source_directory)

    # Return complete regional data structure
    return RegionalData(; regions=regional_data, reef_outlines=canonical_reefs)
end

"""
Initialize regional data with caching support.

Attempts to load from memory cache, then disk cache, before falling back
to full data initialization. Handles cache invalidation and saves new data to disk.

# Arguments
- `config::Dict` : Configuration dictionary containing cache and data settings
- `force_cache_invalidation::Bool` : If true, bypass all caches and reload data
"""
function initialise_data_with_cache(config::Dict; force_cache_invalidation::Bool=false)
    # Access global cache variable
    global REGIONAL_DATA

    # Determine cache directory location
    regional_cache_directory = config["server_config"]["REGIONAL_CACHE_DIR"]

    if !force_cache_invalidation
        # Try memory cache first (fastest)
        local_data = check_existing_regional_data_from_memory()
        if !isnothing(local_data)
            REGIONAL_DATA = local_data
            return nothing
        end

        # Try disk cache second (faster than full reload)
        disk_data = check_existing_regional_data_from_disk(regional_cache_directory)
        if !isnothing(disk_data)
            REGIONAL_DATA = disk_data
            return nothing
        end
    end

    # No cache available or forced invalidation - load from source
    regional_data = initialise_data(config)

    # Update global cache
    REGIONAL_DATA = regional_data

    # Initialize empty tile cache for map rendering
    setup_empty_tile_cache(config)

    # Save to disk for future use
    @debug "Saving regional data cache to disk"
    serialize(
        joinpath(regional_cache_directory, REGIONAL_DATA_CACHE_FILENAME),
        regional_data
    )

    return nothing
end

"""
Get regional data with automatic cache management.

Primary interface for accessing regional data. Handles initialization
and caching automatically.

# Arguments
- `config::Dict` : Configuration dictionary

# Returns
`RegionalData` struct containing all regional information.
"""
function get_regional_data(config::Dict)::RegionalData
    # Ensure data is loaded (with caching)
    initialise_data_with_cache(config)
    # Return cached data
    return REGIONAL_DATA
end

# =============================================================================
# Display and Routing Functions
# =============================================================================

"""
Custom display format for RegionalDataEntry showing key statistics.
"""
function Base.show(io::IO, ::MIME"text/plain", z::RegionalDataEntry)
    println("""
    Criteria: $(names(z.raster_stack))
    Number of valid slope locations: $(nrow(z.slope_table))
    """)
    return nothing
end

"""
Setup HTTP routes for criteria information endpoints.

Creates REST endpoints for accessing regional criteria bounds and metadata.

# Arguments
- `config` : Configuration object
- `auth` : Authentication/authorization handler
"""
function setup_criteria_routes(config, auth)
    regional_data::RegionalData = get_regional_data(config)

    # Endpoint: GET /criteria/{region}/ranges
    # Returns JSON with min/max values for all criteria in specified region
    @get auth("/criteria/{region}/ranges") function (_::Request, region::String)
        @debug "Transforming criteria information to JSON for region $(region)"
        output_dict = OrderedDict()

        # Build lookup dictionary for regional criteria
        regional_criteria_lookup = build_regional_criteria_dictionary(
            regional_data.regions[region]
        )

        # Format each criteria with min/max bounds
        for (id::String, criteria::RegionalCriteriaEntry) in regional_criteria_lookup
            output_dict[id] = OrderedDict(
                :min_val => criteria.bounds.min,
                :max_val => criteria.bounds.max
            )
        end

        return json(output_dict)
    end
end

"""
==========
DEPRECATED 
==========
"""
function criteria_data_map()
    # TODO: Load from config?
    return OrderedDict(
        :Depth => "_bathy",
        :Benthic => "_benthic",
        :Geomorphic => "_geomorphic",
        :Slope => "_slope",
        :Turbidity => "_turbid",
        :WavesHs => "_waves_Hs",
        :WavesTp => "_waves_Tp",
        :Rugosity => "_rugosity",
        :ValidSlopes => "_valid_slopes",
        :ValidFlats => "_valid_flats"
    )
end

function search_criteria()::Vector{String}
    return string.(keys(criteria_data_map()))
end

function site_criteria()::Vector{String}
    return ["SuitabilityThreshold", "xdist", "ydist"]
end

function suitability_criteria()::Vector{String}
    return vcat(search_criteria(), ["SuitabilityThreshold"])
end
