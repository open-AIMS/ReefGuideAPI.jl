# =============================================================================
# Constants and Configuration
# =============================================================================

const REGIONAL_DATA_CACHE_FILENAME = "regional_cache_v2.dat"
const SLOPES_LOOKUP_SUFFIX = "_valid_slopes_lookup.parq"
const DEFAULT_CANONICAL_REEFS_FILE_NAME = "rrap_canonical_outlines.gpkg"

# =============================================================================
# Core Data Structures
# =============================================================================

"""
Metadata container for regional information.

# Fields
- `display_name::String` : Human-readable name for UI display
- `id::String` : Unique system identifier (changing this affects data loading)
- `available_criteria::Vector{String}` : the criteria IDs that are relevant to this region
"""
struct RegionMetadata
    display_name::String
    id::String
    available_criteria::Vector{String}

    function RegionMetadata(;
        display_name::String,
        id::String,
        available_criteria::Vector{String}
    )
        @debug "Creating RegionMetadata" display_name id
        return new(display_name, id, available_criteria)
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
Helper function to create Bounds from a tuple.
"""
function bounds_from_tuple(extrema_tuple::Tuple{T,T}) where T<:Number
    return Bounds(; min=extrema_tuple[1], max=extrema_tuple[2])
end

"""
Metadata for assessment criteria including file naming conventions.

# Fields
- `id::String` : Unique system identifier for the criteria
- `file_suffix::String` : File suffix pattern for data files
- `display_label::String` : Human-readable label for UI display
- `subtitle::String` : Human-readable info about this criteria on subtitle of slider
- `units::String` : Human-readable info about relevant units
- `payload_prefix::String` : The prefix for building the job payload
- `default_bounds::OptionalValue{Bounds}` : The default bounds for the parameter sliders
- `min_tooltip::String` : Tooltip text on min slider
- `max_tooltip::String` : Tooltip text on max slider
"""
struct CriteriaMetadata
    id::String
    file_suffix::String
    display_label::String
    subtitle::String
    units::String
    payload_prefix::String
    default_bounds::OptionalValue{Bounds}
    min_tooltip::String
    max_tooltip::String

    function CriteriaMetadata(;
        id::String,
        file_suffix::String,
        display_label::String,
        subtitle::String,
        units::String,
        payload_prefix::String,
        default_bounds::OptionalValue{Bounds}=nothing,
        min_tooltip::String,
        max_tooltip::String
    )
        return new(
            id,
            file_suffix,
            display_label,
            subtitle,
            units,
            payload_prefix,
            default_bounds,
            min_tooltip,
            max_tooltip
        )
    end
end

# NOTE: This is where you add to list of all possible criteria
const ASSESSMENT_CRITERIA::Dict{String,CriteriaMetadata} = Dict(
    "Depth" => CriteriaMetadata(;
        id="Depth",
        file_suffix="_bathy",
        display_label="Depth",
        subtitle="TODO",
        units="TODO",
        payload_prefix="depth_",
        default_bounds=Bounds(; min=-10, max=-2),
        min_tooltip="TODO",
        max_tooltip="TODO"
    ),
    "Slope" => CriteriaMetadata(;
        id="Slope",
        file_suffix="_slope",
        display_label="Slope",
        subtitle="TODO",
        units="TODO",
        payload_prefix="slope_",
        min_tooltip="TODO",
        max_tooltip="TODO"
    ),
    "Turbidity" => CriteriaMetadata(;
        id="Turbidity",
        file_suffix="_turbid",
        display_label="Turbidity",
        subtitle="TODO",
        units="TODO",
        payload_prefix="turbidity_",
        min_tooltip="TODO",
        max_tooltip="TODO"
    ),
    "WavesHs" => CriteriaMetadata(;
        id="WavesHs",
        file_suffix="_waves_Hs",
        display_label="Wave Height (m)",
        subtitle="TODO",
        units="TODO",
        payload_prefix="waves_height_",
        default_bounds=Bounds(; min=0, max=1),
        min_tooltip="TODO",
        max_tooltip="TODO"
    ),
    "WavesTp" => CriteriaMetadata(;
        id="WavesTp",
        file_suffix="_waves_Tp",
        display_label="Wave Period (s)",
        subtitle="TODO",
        units="TODO",
        payload_prefix="waves_period_",
        default_bounds=Bounds(; min=0, max=6),
        min_tooltip="TODO",
        max_tooltip="TODO"
    ),
    "Rugosity" => CriteriaMetadata(;
        id="Rugosity",
        file_suffix="_rugosity",
        display_label="Rugosity",
        subtitle="TODO",
        units="TODO",
        payload_prefix="rugosity_",
        min_tooltip="TODO",
        max_tooltip="TODO"
    )
)

# Create list from the dictionary values
const ASSESSMENT_CRITERIA_LIST = collect(values(ASSESSMENT_CRITERIA))

"""
Combines criteria metadata with regional boundary values.

# Fields
- `metadata::CriteriaMetadata` : Criteria definition and metadata
- `bounds::Bounds` : Min/max values for this criteria in the region
"""
struct BoundedCriteria
    metadata::CriteriaMetadata
    bounds::Bounds

    function BoundedCriteria(;
        metadata::CriteriaMetadata,
        bounds::Bounds
    )
        return new(metadata, bounds)
    end
end

# Maps criteria ID to bounded criteria
const BoundedCriteriaDict = Dict{String,BoundedCriteria}

"""
Complete data package for a single region including rasters and metadata.

Validates during construction that all criteria instantiated have a
corresponding layer in the RasterStack

# Fields
- `region_id::String` : Unique identifier for the region
- `region_metadata::RegionMetadata` : Display metadata for the region
- `raster_stack::Rasters.RasterStack` : Geospatial raster data layers
- `slope_table::DataFrame` : Coordinates and values for valid slope reef
  locations
- `criteria::Dict{String, BoundedCriteria}` : Computed criteria bounds for this region
"""
struct RegionalDataEntry
    region_id::String
    region_metadata::RegionMetadata
    raster_stack::Rasters.RasterStack
    slope_table::DataFrame
    criteria::BoundedCriteriaDict

    function RegionalDataEntry(;
        region_id::String,
        region_metadata::RegionMetadata,
        raster_stack::Rasters.RasterStack,
        slope_table::DataFrame,
        criteria::BoundedCriteriaDict
    )
        # Get available layers and expected criteria from metadata
        raster_layer_names = Set(string.(names(raster_stack)))
        expected_criteria_set = Set(region_metadata.available_criteria)

        # Collect criteria that are actually instantiated (non-nothing)
        instantiated_criteria = String[]
        missing_from_criteria = String[]
        missing_from_rasters = String[]

        # Check each criteria field for instantiation
        for desired_criteria in expected_criteria_set
            criteria_entry = get(criteria, desired_criteria, nothing)
            if !isnothing(criteria_entry)
                layer_id = criteria_entry.metadata.id
                push!(instantiated_criteria, layer_id)

                # Validate this instantiated criteria has a corresponding raster layer
                if layer_id ∉ raster_layer_names
                    push!(missing_from_rasters, layer_id)
                end
            end
        end

        # Cross-validate: ensure all expected criteria from metadata are instantiated
        for expected_criteria_id in expected_criteria_set
            if expected_criteria_id ∉ instantiated_criteria
                push!(missing_from_criteria, expected_criteria_id)
            end
        end

        # Cross-validate: ensure all raster layers correspond to expected criteria
        unexpected_layers = String[]
        for layer_name in raster_layer_names
            if layer_name ∉ expected_criteria_set
                push!(unexpected_layers, layer_name)
            end
        end

        # Report validation errors
        validation_errors = String[]

        if !isempty(missing_from_rasters)
            push!(
                validation_errors,
                "RasterStack missing layers for instantiated criteria: $(join(missing_from_rasters, ", "))"
            )
        end

        if !isempty(missing_from_criteria)
            push!(
                validation_errors,
                "RegionalCriteria missing expected criteria from metadata: $(join(missing_from_criteria, ", "))"
            )
        end

        if !isempty(unexpected_layers)
            push!(
                validation_errors,
                "RasterStack contains unexpected layers not in metadata: $(join(unexpected_layers, ", "))"
            )
        end

        # If any validation errors, log and throw
        if !isempty(validation_errors)
            @error "Validation failed for region $(region_metadata.display_name)" region_id validation_errors available_in_metadata = collect(
                expected_criteria_set
            ) instantiated_criteria available_raster_layers = collect(raster_layer_names)

            error("RegionalDataEntry validation failed: $(join(validation_errors, "; "))")
        end

        @info "Created RegionalDataEntry for $(region_metadata.display_name)" region_id slope_locations = nrow(
            slope_table
        ) raster_layers = length(names(raster_stack)) validated_criteria_layers = length(
            instantiated_criteria
        ) expected_criteria = length(expected_criteria_set)

        return new(region_id, region_metadata, raster_stack, slope_table, criteria)
    end
end

# Type alias
const RegionalDataDict = Dict{String,RegionalDataEntry}

"""
Top-level container for all regional data and reef outlines.

# Fields
- `regions::RegionalDataDict` : Regional data indexed by region ID
- `reef_outlines::DataFrame` : Canonical reef outline geometries
"""
struct RegionalData
    regions::RegionalDataDict
    reef_outlines::DataFrame

    function RegionalData(;
        regions::RegionalDataDict,
        reef_outlines::DataFrame
    )
        total_locations = sum(nrow(entry.slope_table) for entry in values(regions))
        @info "RegionalData initialized" num_regions = length(regions) total_valid_locations =
            total_locations reef_outlines = nrow(reef_outlines)
        return new(regions, reef_outlines)
    end
end

# =============================================================================
# Configuration Constants
# =============================================================================

# Normal list - only Townsville has rugosity!
const BASE_CRITERIA_IDS::Vector{String} = [
    criteria.id for criteria in values(ASSESSMENT_CRITERIA) if criteria.id != "Rugosity"
]
# All criteria
const ALL_CRITERIA_IDS::Vector{String} = [
    criteria.id for criteria in values(ASSESSMENT_CRITERIA)
]

# Define all available regions for the assessment system NOTE: Here is where you
# configure which criteria are available per region
const REGIONS::Vector{RegionMetadata} = [
    RegionMetadata(;
        display_name="Townsville/Whitsunday Management Area",
        id="Townsville-Whitsunday",
        available_criteria=ALL_CRITERIA_IDS
    ),
    RegionMetadata(;
        display_name="Cairns/Cooktown Management Area",
        id="Cairns-Cooktown",
        available_criteria=BASE_CRITERIA_IDS
    ),
    RegionMetadata(;
        display_name="Mackay/Capricorn Management Area",
        id="Mackay-Capricorn",
        available_criteria=BASE_CRITERIA_IDS
    ),
    RegionMetadata(;
        display_name="Far Northern Management Area",
        id="FarNorthern",
        available_criteria=BASE_CRITERIA_IDS
    )
]

# GLOBAL variable to store regional data cache
REGIONAL_DATA::OptionalValue{RegionalData} = nothing

# =============================================================================
# Data Loading and Processing Functions
# =============================================================================

"""
Build regional criteria bounds from slope table data.

Computes min/max bounds for each assessment criteria by finding extrema
in the slope table data columns. Only processes criteria that are available
for the specific region as defined in the region metadata.

# Arguments
- `table::DataFrame` : Slope table containing criteria data columns
- `region_metadata::RegionMetadata` : Region metadata specifying available criteria

# Returns
`BoundedCriteriaDict` with computed bounds for relevant criteria
"""
function derive_criteria_bounds_from_slope_table(
    table::DataFrame,
    region_metadata::RegionMetadata
)::BoundedCriteriaDict
    @debug "Computing assessment criteria bounds from slope table" table_rows = nrow(table) region_id =
        region_metadata.id available_criteria = region_metadata.available_criteria

    # Helper function to compute bounds for a specific criteria
    function compute_criteria_bounds(
        criteria::CriteriaMetadata
    )::Union{BoundedCriteria,Nothing}
        if criteria.id ∈ region_metadata.available_criteria
            if hasproperty(table, Symbol(criteria.id))
                bounds = bounds_from_tuple(extrema(table[:, criteria.id]))
                @debug "Computed $(criteria.display_label) bounds" range = "$(bounds.min):$(bounds.max)"
                return BoundedCriteria(; metadata=criteria, bounds=bounds)
            else
                @error "Region metadata lists $(criteria.display_label) with id $(criteria.id) as available but column missing from slope table" region_id =
                    region_metadata.id column = criteria.id
                throw(
                    ErrorException(
                        "Missing required column '$(criteria.id)' in slope table for region $(region_metadata.id)"
                    )
                )
            end
        end
        return nothing
    end

    # For each criteria, if available in region, try to find bounds
    criteria_dict::BoundedCriteriaDict = Dict()
    for criteria in values(ASSESSMENT_CRITERIA)
        bounded_criteria = compute_criteria_bounds(criteria)
        # Only include criteria relevant to the region
        if !isnothing(bounded_criteria)
            criteria_dict[criteria.id] = bounded_criteria
        end
    end

    @debug "Completed computation of criteria bounds, returning populated BoundedCriteriaDict"
    return criteria_dict
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
    pattern = "$(region.id)*$(criteria.file_suffix).tif"
    @debug "Searching for data file" pattern directory = data_source_directory

    # Search for files matching the pattern: {region_id}*{criteria_suffix}.tif
    matched_files = glob(pattern, data_source_directory)

    # Validate exactly one match exists
    if length(matched_files) == 0
        @error "No data file found for criteria" criteria_id = criteria.id region_id =
            region.id pattern
        throw(ErrorException("Did not find data for the criteria: $(criteria.id)."))
    end
    if length(matched_files) > 1
        @error "Multiple data files found for criteria - ambiguous match" criteria_id =
            criteria.id region_id = region.id matched_files
        throw(
            ErrorException(
                "Found more than one data source match for criteria: $(criteria.id). This is ambiguous, unsure how to proceed."
            )
        )
    end

    file_path = first(matched_files)
    @debug "Found data file" criteria_id = criteria.id file_path
    return file_path
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
    file_path = joinpath(source_dir, file_name)
    @info "Loading canonical reef outlines" file_path
    try
        reef_data = GDF.read(file_path)
        @info "Successfully loaded reef outlines" num_reefs = nrow(reef_data)
        return reef_data
    catch e
        @error "Failed to load canonical reef outlines" file_path error = e
        rethrow(e)
    end
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
        @info "Using cached regional data from memory"
        return REGIONAL_DATA
    end
    @debug "No regional data found in memory cache"
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
        @info "Loading regional data from disk cache" cache_file = reg_cache_filename
        try
            data = deserialize(reg_cache_filename)
            @info "Successfully loaded regional data from disk cache"
            return data
        catch err
            @warn "Failed to deserialize regional data cache - removing corrupted file" cache_file =
                reg_cache_filename error = err
            # Remove corrupted cache file
            rm(reg_cache_filename)
        end
    else
        @debug "No disk cache file found" expected_path = reg_cache_filename
    end
    # No cache available or load failed
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
    @info "Starting regional data initialization from source files"

    regional_data::RegionalDataDict = Dict()
    data_source_directory = config["prepped_data"]["PREPPED_DATA_DIR"]
    @info "Using data source directory" directory = data_source_directory

    # Process each region sequentially
    for region_metadata::RegionMetadata in REGIONS
        @info "Processing region" region = region_metadata.display_name region_id =
            region_metadata.id

        # Initialize data collection arrays
        data_paths = String[]
        data_names = String[]

        # Load slope table containing valid reef coordinates and criteria values
        slope_filename = get_slope_parquet_filename(region_metadata)
        slope_file_path = joinpath(data_source_directory, slope_filename)
        @debug "Loading slope table" file_path = slope_file_path

        try
            slope_table::DataFrame = GeoParquet.read(slope_file_path)
            @info "Loaded slope table" region_id = region_metadata.id num_locations = nrow(
                slope_table
            )

            # Add coordinate columns for spatial referencing
            add_lat_long_columns_to_dataframe(slope_table)

            # Filter criteria list to only those available for this region
            available_criteria::Vector{String} = region_metadata.available_criteria
            region_criteria_list::Vector{CriteriaMetadata} = [
                ASSESSMENT_CRITERIA[id] for id in available_criteria
            ]

            @debug "Filtered criteria for region" region_id = region_metadata.id total_criteria = length(
                ASSESSMENT_CRITERIA
            ) available_criteria = length(region_criteria_list) criteria_ids = [
                c.id for c in region_criteria_list
            ]

            # Collect raster file paths for available criteria only
            for criteria::CriteriaMetadata in region_criteria_list
                @debug "Processing criteria" criteria_id = criteria.id region_id =
                    region_metadata.id

                # Find the corresponding .tif file for this criteria
                data_file_path = find_data_source_for_criteria(;
                    data_source_directory,
                    region=region_metadata,
                    criteria
                )

                push!(data_paths, data_file_path)
                # Use criteria ID as the raster layer name
                push!(data_names, criteria.id)
            end

            @info "Found all criteria data files for region" region_id = region_metadata.id num_criteria = length(
                data_paths
            ) available_criteria = join([c.id for c in region_criteria_list], ", ")

            # Compute regional criteria bounds from slope table data
            bounds::BoundedCriteriaDict = derive_criteria_bounds_from_slope_table(
                slope_table, region_metadata
            )

            # Create lazy-loaded raster stack from all criteria files
            @debug "Creating raster stack" region_id = region_metadata.id num_layers = length(
                data_paths
            )
            raster_stack = RasterStack(data_paths; name=data_names, lazy=true)

            # Store complete regional data entry
            regional_data[region_metadata.id] = RegionalDataEntry(;
                region_id=region_metadata.id,
                region_metadata,
                raster_stack,
                slope_table,
                criteria=bounds
            )

        catch e
            @error "Failed to process region data" region_id = region_metadata.id error = e
            rethrow(e)
        end
    end

    @info "Completed processing all regions" num_regions = length(regional_data)

    # Load canonical reef outlines that apply to all regions
    canonical_reefs = load_canonical_reefs(data_source_directory)

    # Return complete regional data structure
    @info "Regional data initialization completed successfully"
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
    @info "Initializing regional data with caching" force_cache_invalidation

    # Access global cache variable
    global REGIONAL_DATA

    # Determine cache directory location
    regional_cache_directory = config["server_config"]["REGIONAL_CACHE_DIR"]
    @debug "Using cache directory" directory = regional_cache_directory

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
    else
        @info "Cache invalidation forced - reloading from source files"
    end

    # No cache available or forced invalidation - load from source
    @info "Loading regional data from source files (no cache available)"
    regional_data = initialise_data(config)

    # Update global cache
    REGIONAL_DATA = regional_data

    # Save to disk for future use
    @info "Saving regional data to disk cache" cache_directory = regional_cache_directory
    try
        serialize(
            joinpath(regional_cache_directory, REGIONAL_DATA_CACHE_FILENAME),
            regional_data
        )
        @info "Successfully saved regional data cache to disk"
    catch e
        @warn "Failed to save regional data cache to disk" error = e
    end

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
    @debug "Getting regional data with automatic cache management"
    # Ensure data is loaded (with caching)
    initialise_data_with_cache(config)
    # Return cached data
    return REGIONAL_DATA
end

# =============================================================================
# Display Methods
# =============================================================================

"""
Enhanced display format for RegionalDataEntry showing key statistics.
"""
function Base.show(io::IO, ::MIME"text/plain", entry::RegionalDataEntry)
    println(io, "RegionalDataEntry: $(entry.region_metadata.display_name)")
    println(io, "  Region ID: $(entry.region_id)")
    println(io, "  Raster layers: $(join(names(entry.raster_stack), ", "))")
    println(io, "  Valid slope locations: $(nrow(entry.slope_table))")
    println(
        io, "  Available criteria: $(join(entry.region_metadata.available_criteria, ", "))"
    )
    println(io, "  Criteria bounds:")

    # Show each criteria with its bounds from the dictionary
    for (_, bounded_criteria) in entry.criteria
        min_val = round(bounded_criteria.bounds.min; digits=2)
        max_val = round(bounded_criteria.bounds.max; digits=2)
        display_label = bounded_criteria.metadata.display_label
        println(io, "    $(display_label): $(min_val) - $(max_val)")
    end

    # Show criteria that are expected but not available
    expected_criteria = Set(entry.region_metadata.available_criteria)
    available_criteria = Set(keys(entry.criteria))
    missing_criteria = setdiff(expected_criteria, available_criteria)

    for criteria_id in missing_criteria
        if haskey(ASSESSMENT_CRITERIA, criteria_id)
            display_label = ASSESSMENT_CRITERIA[criteria_id].display_label
            println(io, "    $(display_label): Not available")
        end
    end
end

"""
Display format for RegionalData showing system overview.
"""
function Base.show(io::IO, ::MIME"text/plain", data::RegionalData)
    total_locations = sum(nrow(entry.slope_table) for entry in values(data.regions))

    println(io, "RegionalData:")
    println(io, "  Regions: $(length(data.regions))")
    println(io, "  Total valid locations: $(total_locations)")
    println(io, "  Reef outlines: $(nrow(data.reef_outlines))")
    println(io, "")

    println(io, "Regional breakdown:")
    for (_, region_entry) in data.regions
        locations = nrow(region_entry.slope_table)
        println(
            io, "  $(region_entry.region_metadata.display_name): $(locations) locations"
        )
    end
    println(io, "")

    return println(
        io,
        "Assessment criteria: $(join([c.display_label for c in values(ASSESSMENT_CRITERIA)], ", "))"
    )
end
