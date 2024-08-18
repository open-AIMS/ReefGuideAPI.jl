"""Geometry-based assessment methods."""


include("geom_ops.jl")


function identify_potential_sites(rst_stack, scan_locs, threshold, geom)
    res = abs(step(dims(rst_stack, X)))
    geom_buff = GO.buffer(geom, res)

    # TODO: Dynamically build this ruleset
    ruleset = Dict(
        :Depth => (data) -> within_thresholds(data, -9.0, -2.0),
        :WavesTp => (data) -> within_thresholds(data, 0.0, 5.9)
    )

    # Search each location to assess
    best_score = zeros(length(scan_locs))
    best_poly = Vector(undef, length(scan_locs))
    best_degree = zeros(Int64, length(scan_locs))
    for (i, (lon_idx, lat_idx)) in enumerate(scan_locs)
        # Move geom to new centroid
        lon = dims(rst_stack, X)[lon_idx]
        lat = dims(rst_stack, Y)[lat_idx]
        mv_geom = move_geom(geom_buff, (lon, lat))

        b_score, b_deg, b_poly = assess_reef_site(rst_stack, mv_geom, ruleset; degree_step=15.0)

        best_score[i] = b_score
        best_degree[i] = b_deg
        best_poly[i] = b_poly
    end

    return DataFrame(score=best_score, orientation=best_degree, poly=best_poly)
end

function assess_reef_site(rst, geom, ruleset; degree_step=15.0)
    n_rotations = length(0.0:degree_step:359.0)
    score = zeros(n_rotations)
    best_poly = Vector(undef, n_rotations)

    for (j, r) in enumerate(0.0:degree_step:359.0)
        rot_geom = rotate_geom(geom, r)
        c_rst = crop(rst; to=rot_geom)
        if !all(size(c_rst) .> (0, 0))
            @warn "No data found!"
            continue
        end

        window = trues(size(c_rst))
        for (n, crit_rule) in ruleset
            window .= window .& crit_rule(c_rst[n])
            if count(window) < ceil(length(window) / 3)
                # Stop checking other rules if below hard threshold
                break
            end
        end

        score[j] = mean(window)
        best_poly[j] = rot_geom
    end

    return score[argmax(score)], argmax(score)-1, best_poly[argmax(score)]
end
