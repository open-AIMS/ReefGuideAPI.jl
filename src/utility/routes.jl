# =============================================================================
# Routes - utility/misc routes not related to jobs or key tasks
# =============================================================================

"""
Setup HTTP routes for criteria information endpoints.

Creates REST endpoints for accessing regional criteria bounds and metadata.

# Arguments
- `config` : Configuration object
- `auth` : Authentication/authorization handler
"""
function setup_utility_routes(config, auth)
    @info "Setting up criteria routes"
    regional_data::RegionalData = get_regional_data(config)

    # Health check
    @get "/health" function ()
        return json(Dict(:status => "healthy"))
    end

    # Endpoint: GET /criteria/{region}/ranges
    # Returns JSON with min/max values for all criteria in specified region
    @get auth("/criteria/{region}/ranges") function (_::Request, region::String)
        @info "Processing criteria ranges request" region

        if !haskey(regional_data.regions, region)
            @warn "Request for unknown region" region available_regions = keys(
                regional_data.regions
            )
            return json(Dict("error" => "Region not found"))
        end

        @debug "Transforming criteria information to JSON for region $(region)"
        output_dict = OrderedDict()

        # Build lookup dictionary for regional criteria
        regional_criteria_lookup = regional_data.regions[region].criteria

        # Format each criteria with min/max bounds
        for (id::String, criteria::BoundedCriteria) in regional_criteria_lookup
            # build default min/max
            default_bounds::Bounds = something(
                criteria.metadata.default_bounds, criteria.bounds
            )

            output_dict[id] = OrderedDict(
                # Unique ID (and data field name)
                :id => id,

                # min/max bounds
                :min_val => criteria.bounds.min,
                :max_val => criteria.bounds.max,

                # display info
                :display_title => criteria.metadata.display_label,
                :display_subtitle => criteria.metadata.description,
                :units => criteria.metadata.units,

                # default min/max
                :default_min_val => default_bounds.min,
                :default_max_val => default_bounds.max,

                # how to build a job payload (prefix of job i.e. depth_) then build depth_min depth_max
                :payload_property_prefix => criteria.metadata.payload_prefix
            )
        end

        @debug "Returning criteria ranges" region num_criteria = length(output_dict)
        return json(output_dict)
    end

    """Obtain the spatial bounds for a given region of interest"""
    @get auth("/bounds/{region}") function (req::Request, region::String)
        rst_stack = regional_data.regions[region].raster_stack

        return json(Rasters.bounds(rst_stack))
    end

    @info "Criteria routes setup completed"
end
