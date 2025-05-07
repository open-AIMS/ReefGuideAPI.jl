
function setup_job_routes(config, auth)
    reg_assess_data = setup_regional_data(config)

    @get auth("/job/details/{job_id}") function (req::Request, job_id::String)
        srv = DiskService(_cache_location(config))
        return json(job_details(srv, job_id))
    end

    @get auth("/job/result/{job_id}") function (req::Request, job_id::String)
        srv = DiskService(_cache_location(config))
        return file(job_result(srv, job_id))
    end

    @get auth("/submit/region-assess/{reg}/{rtype}") function (
        req::Request, reg::String, rtype::String
    )
        # 127.0.0.1:8000/submit/region-assess/Mackay-Capricorn/slopes?Depth=-12.0:-2.0&Slope=0.0:40.0&Rugosity=0.0:6.0&SuitabilityThreshold=95
        qp = queryparams(req)
        srv = DiskService(_cache_location(config))
        job_id = create_job_id(qp) * "$(reg)_suitable"

        details = job_details(srv, job_id)
        job_state = job_status(details)
        if (job_state == "no job") && (job_state != "processing")
            @debug "$(now()) : Submitting $(job_id)"
            assessed_fn = cache_filename(
                extract_criteria(qp, suitability_criteria()),
                config,
                "$(reg)_suitable",
                "tiff"
            )

            details = submit_job(srv, job_id, assessed_fn)

            # Do job asyncronously...
            @async assess_region(
                config,
                qp,
                reg,
                rtype,
                reg_assess_data
            )
        end

        if job_state == "processing"
            @debug details
        end

        @debug "$(now()) : Job submitted, return polling url"
        return json(details.access_url)
    end

    @get auth("/submit/site-assess/{reg}/{rtype}") function (
        req::Request, reg::String, rtype::String
    )
        qp = queryparams(req)
        suitable_sites_fn = cache_filename(
            qp, config, "$(reg)_potential_sites", "geojson"
        )

        srv = DiskService(_cache_location(config))
        job_id = create_job_id(qp) * "$(reg)_potential_sites"

        details = job_details(srv, job_id)
        job_state = job_status(details)
        if (job_state != "no job") && (job_state != "completed")
            @debug "$(now()) : Submitting $(job_id)"
            details = submit_job(srv, job_id, suitable_sites_fn)

            # Do job asyncronously...
            @async begin
                assessed_fn = assess_region(
                    config,
                    qp,
                    reg,
                    rtype,
                    reg_assess_data
                )

                assessed = Raster(assessed_fn; missingval=0)

                # Extract criteria and assessment
                pixel_criteria = extract_criteria(qp, search_criteria())
                deploy_site_criteria = extract_criteria(qp, site_criteria())

                best_sites = filter_sites(
                    assess_sites(
                        reg_assess_data, reg, rtype, pixel_criteria, deploy_site_criteria,
                        assessed
                    )
                )

                # Specifically clear from memory to invoke garbage collector
                assessed = nothing

                if nrow(best_sites) == 0
                    open(suitable_sites_fn, "w") do f
                        JSON.print(f, nothing)
                    end
                else
                    output_geojson(suitable_sites_fn, best_sites)
                end

                details.status = "completed"
                update_job!(srv, job_id)
            end
        end

        @debug "$(now()) : Job submitted, return polling url"
        return json(details.access_url)
    end
end

"""
    setup_region_routes(config, auth)

Set up endpoints for regional assessment.
"""
function setup_region_routes(config, auth)
    reg_assess_data = setup_regional_data(config)

    @get auth("/assess/{reg}/{rtype}") function (req::Request, reg::String, rtype::String)
        qp = queryparams(req)
        mask_path = cache_filename(qp, config, "", "tiff")
        if isfile(mask_path)
            return file(mask_path; headers=COG_HEADERS)
        end

        # Otherwise, create the file
        @debug "$(now()) : Assessing criteria"
        mask_data = mask_region(reg_assess_data, reg, qp, rtype)

        @debug "$(now()) : Running on thread $(threadid())"
        @debug "Writing to $(mask_path)"
        # Writing time: ~10-25 seconds
        _write_cog(mask_path, UInt8.(mask_data), config)

        return file(mask_path; headers=COG_HEADERS)
    end

    @get auth("/suitability/assess/{reg}/{rtype}") function (
        req::Request, reg::String, rtype::String
    )
        # somewhere:8000/suitability/assess/region-name/reeftype?criteria_names=Depth,Slope&lb=-9.0,0.0&ub=-2.0,40
        # 127.0.0.1:8000/suitability/assess/Cairns-Cooktown/slopes?Depth=-4.0:-2.0&Slope=0.0:40.0&Rugosity=0.0:6.0
        # 127.0.0.1:8000/suitability/assess/Cairns-Cooktown/slopes?Depth=-4.0:-2.0&Slope=0.0:40.0&Rugosity=0.0:6.0&SuitabilityThreshold=95
        # 127.0.0.1:8000/suitability/assess/Mackay-Capricorn/slopes?Depth=-12.0:-2.0&Slope=0.0:40.0&Rugosity=0.0:6.0&SuitabilityThreshold=95
        qp = queryparams(req)
        assessed_fn = assess_region(config, qp, reg, rtype, reg_assess_data)

        return file(assessed_fn; headers=COG_HEADERS)
    end

    @get auth("/suitability/site-suitability/{reg}/{rtype}") function (
        req::Request, reg::String, rtype::String
    )
        # 127.0.0.1:8000/suitability/site-suitability/Cairns-Cooktown/slopes?Depth=-4.0:-2.0&Slope=0.0:40.0&Rugosity=0.0:6.0&SuitabilityThreshold=95&xdist=450&ydist=50
        # 127.0.0.1:8000/suitability/site-suitability/Mackay-Capricorn/slopes?Depth=-12.0:-2.0&Slope=0.0:40.0&Rugosity=0.0:6.0&SuitabilityThreshold=95&xdist=100&ydist=100
        qp = queryparams(req)
        suitable_sites_fn = cache_filename(
            qp, config, "$(reg)_potential_sites", "geojson"
        )
        if isfile(suitable_sites_fn)
            return file(suitable_sites_fn)
        end

        # Assess location suitability if needed
        assessed_fn = assess_region(config, qp, reg, rtype, reg_assess_data)
        assessed = Raster(assessed_fn; lazy=true, missingval=0)

        # Extract criteria and assessment
        pixel_criteria = extract_criteria(qp, search_criteria())
        deploy_site_criteria = extract_criteria(qp, site_criteria())

        best_sites = filter_sites(
            assess_sites(
                reg_assess_data, reg, rtype, pixel_criteria, deploy_site_criteria, assessed
            )
        )

        # Specifically clear from memory to invoke garbage collector
        assessed = nothing

        if nrow(best_sites) == 0
            open(suitable_sites_fn, "w") do f
                JSON.print(f, nothing)
            end
        else
            output_geojson(suitable_sites_fn, best_sites)
        end

        return file(suitable_sites_fn)
    end

    """Obtain the spatial bounds for a given region of interest"""
    @get auth("/bounds/{reg}") function (req::Request, reg::String)
        rst_stack = reg_assess_data[reg].stack

        return json(Rasters.bounds(rst_stack))
    end

    # Form for testing/dev
    # https:://somewhere:8000/suitability/assess/region-name/reeftype?criteria_names=Depth,Slope&lb=-9.0,0.0&ub=-2.0,40
    @get "/" function ()
        return html("""
               <form action="/assess/Cairns-Cooktown/slopes" method="post">
                   <label for="criteria_names">Criteria Names:</label><br>
                   <input type="text" id="criteria_names" name="criteria"><br>

                   <label for="lb">Lower Bound:</label><br>
                   <input type="text" id="lb" name="lower_bound"><br><br>

                   <label for="ub">Upper Bound:</label><br>
                   <input type="text" id="ub" name="upper_bound"><br><br>

                   <input type="submit" value="Submit">
               </form>
               """)
    end

    # Parse the form data and return it
    @post auth("/form") function (req)
        data = formdata(req)
        return data
    end
end
