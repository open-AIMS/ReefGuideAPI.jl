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
function setup_criteria_routes(config, auth)
    @info "Setting up criteria routes"
    regional_data::RegionalData = get_regional_data(config)

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

        @debug "Returning criteria ranges" region num_criteria = length(output_dict)
        return json(output_dict)
    end

    @info "Criteria routes setup completed"
end
