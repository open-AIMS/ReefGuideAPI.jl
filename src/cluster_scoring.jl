"""
Cluster scoring methods for assessing defined physical and ecological constraints.
* Note: All scoring functions have been standardised to output metrics in km/sqkm units
"""

using Statistics

import ArchGDAL as AG
import GeometryOps as GO
import GeoInterface as GI

import GeoDataFrames as GDF

include("site_assessment/geom_ops.jl")

## Functions for cluster scoring
"""
    proportionate_score(x, bound; rev=false)

Calculates a score for `x` value based on whether it meets `bound`. This is designed where
`bound` is a maximum value for suitability (i.e. x is suitable if x < bound). If `x` is
outside of `bound` then a proportionate score is calculated. Keyword `rev` is intended to switch
the function to using a minimum bound for suitability.

# Arguments
- `x` : target value to score
- `bound` : Threshold for suitability (if rev = true or unspecified, then bound is a max value).
- `rev` : Switches scoring from using a maximum bound threshold to using a minumum bound threshold
(i.e. switches to only calculating a proportionate score if `x` < `bound`).
"""
function proportionate_score(x, bound; rev=false)
    if rev
        if x >= bound
            return 0.0
        end

        return bound / x
    else
        if x <= bound
            return 0.0
        end

        return x / bound
    end
end

"""
    score_clusters(
        cluster_dataframe,
        max_side_bound,
        max_area_bound,
        max_distance_bound,
        max_benthic_bound,
        min_benthic_bound,
        min_reef_number;
        output_polygons=false
    )

Score a given assignment of reef clusters according to maximum physical thresholds such as
`max_side_bound`, `max_area_bound`, `max_distance_bound`, `max_benthic_bound`. Assessment of
cluster dimensions and area is based on rectangular shapes.

# Arguments
- `cluster_dataframe` : DataFrame containing reef outlines and cluster assignments for each reef in `clusters` column.
- `raster` : Boolean raster indicating which pixels contain suitable habitat.
- `max_side_bound` : Maximum length of a side of a cluster. (Units=km).
- `max_area_bound` : Maximum are of a single cluster. (Units=squared km).
- `max_distance_bound` : Maximum distance from one reef to it's nearest reef (Units=km).
- `max_benthic_bound` : Maximum area of suitable benthic habitat (Units=squared km).
- `min_benthic_bound` : Minimum area of suitable benthic habitat (Units=squared km).
- `min_reef_number` : Minimum number of reef polygons per cluster.

# Returns
DataFrame with a score for each cluster. Scores == 0.0 mean a cluster satisfies all criteria,
scores > 0.0 mean a cluster does not satisfy all criteria, with higher scores meaning values
are increasing, away from the criteria bounds.
"""
function score_clusters(
    cluster_dataframe,
    polygon_distance_matrix,
    max_side_bound,
    max_area_bound,
    max_distance_bound,
    max_benthic_bound,
    min_benthic_bound,
    min_reef_number;
    output_polygons=false
)
    cluster_scores = DataFrame(
        clusters = Vector{Union{Missing, Int64}}(missing, size(unique(cluster_dataframe.clusters),1)),
        score =  Vector{Union{Missing, Float64}}(missing, size(unique(cluster_dataframe.clusters),1))
    )

    # Also include the option to add cluster rectangular polygons to the output, this is intended for plotting
    if output_polygons
        cluster_scores = DataFrame(
            clusters = Vector{Union{Missing, Int64}}(missing, size(unique(cluster_dataframe.clusters),1)),
            score =  Vector{Union{Missing, Float64}}(missing, size(unique(cluster_dataframe.clusters),1)),
            cluster_poly = Vector{Union{Missing, GeoInterface.Wrappers.Polygon}}(missing, size(unique(cluster_dataframe.clusters),1))
        )
    end

    grouped_gdf = groupby(cluster_dataframe, :clusters)

    for (i, groupdf) in enumerate(grouped_gdf)
        side_length = cluster_side_distance(groupdf)

        total_rectangular_area = cluster_rectangular_area(groupdf)

        benthic_area = sum(groupdf.benthic_area)

        cluster_score = proportionate_score(side_length, max_side_bound) +
            proportionate_score(total_rectangular_area, max_area_bound) +
            proportionate_score(benthic_area, max_benthic_bound) +
            proportionate_score(benthic_area, min_benthic_bound; rev=true) +
            proportionate_score(size(groupdf, 1), min_reef_number; rev=true)

        cluster_scores[i, :clusters] = first(unique(groupdf.clusters))
        cluster_scores[i, :score] = cluster_score

        if output_polygons
            cluster_scores[i, :cluster_poly] = cluster_polygons(groupdf)
        end
    end

    cluster_distances = [
        distance_to_next_reef(
            cluster_i,
            cluster_dataframe.clusters,
            polygon_distance_matrix
        ) for cluster_i in cluster_scores.clusters
    ]
    cluster_scores.score = cluster_scores.score .+ proportionate_score.(cluster_distances, [max_distance_bound])

    external_touching_reefs = [
        touching_polygon_cluster_scoring(
            cluster_i,
            cluster_dataframe.clusters,
            polygon_distance_matrix
        ) for cluster_i in cluster_scores.clusters
    ]
    cluster_scores.score = cluster_scores.score .+ external_touching_reefs

    return cluster_scores
end

"""
    cluster_side_distance(groupdf)

Finds the length of the longest side of a cluster (assessing a cluster via rectangular shape).
All scoring functions have been standardised to return kilometre units.
"""
function cluster_side_distance(groupdf)
    lons = vcat([first.(GI.getpoint(x)) for x in groupdf.geometry]...)
    lon_bounds = (minimum(lons), maximum(lons))
    lats = vcat([last.(GI.getpoint(x)) for x in groupdf.geometry]...)
    lat_bounds = (minimum(lats), maximum(lats))

    lon_dist, lat_dist = x_y_dimensions(lon_bounds, lat_bounds)
    max_side_length = maximum([lon_dist, lat_dist])

    return max_side_length / 1000
end

"""
    cluster_rectangular_area(groupdf)

Calculate the rectangular area of a cluster. Assuming a cluster area occupies space from
minimum(cluster lon, cluster lat) to maximum(cluster lon, cluster lat).
All scoring functions have been standardised to return kilometre/squarekm units.
"""
function cluster_rectangular_area(groupdf)
    lons = vcat([first.(GI.getpoint(x)) for x in groupdf.geometry]...)
    lon_bounds = (minimum(lons), maximum(lons))
    lats = vcat([last.(GI.getpoint(x)) for x in groupdf.geometry]...)
    lat_bounds = (minimum(lats), maximum(lats))

    lon_dist, lat_dist = x_y_dimensions(lon_bounds, lat_bounds)

    cluster_area = (lon_dist / 1000) * (lat_dist / 1000)

    return cluster_area
end

"""
    x_y_dimensions(lon_bounds, lat_bounds)

Returns the length and width of a rectangular cluster in the unit of the geometry (dataframe) CRS.
"""
function x_y_dimensions(lon_bounds, lat_bounds)
    lon_dist = euclidean((first(lon_bounds), lat_bounds[1]), (last(lon_bounds), lat_bounds[1]))
    lat_dist = euclidean((lon_bounds[1], first(lat_bounds)), (lon_bounds[1], last(lat_bounds)))

    return lon_dist, lat_dist
end

"""
    cluster_polygons(groupdf; target_crs=EPSG(9473))

Returns a rectangle box polygon that corresponds to the target cluster.
"""
function cluster_polygons(groupdf; target_crs=EPSG(9473))
    lons = vcat([first.(GI.getpoint(x)) for x in groupdf.geometry]...)
    lon_bounds = (minimum(lons), maximum(lons))
    lats = vcat([last.(GI.getpoint(x)) for x in groupdf.geometry]...)
    lat_bounds = (minimum(lats), maximum(lats))

    polygons = create_poly(create_bbox(lon_bounds, lat_bounds), target_crs)

    return polygons
end

"""
    touching_polygon_cluster_scoring(
        cluster::Int64,
        clusters::Vector{Int64},
        distance_matrix::Matrix{Float64}
    )

Find the number of external polygons that are touching polygons within the target cluster.

# Arguments
- `cluster` : Target cluster to investigate.
- `clusters` : Vector of all clusters assigned to each reef.
- `distance_matrix` : n*n distance matrix of distances between reef polygons.
"""
function touching_polygon_cluster_scoring(
    cluster::Int64,
    clusters::Vector{Int64},
    distance_matrix::Matrix{Float64}
)
    ind = clusters .== cluster
    if any(distance_matrix[ind, Not(ind)] .== 0.0) # Check if any of the reefs in the target cluster are touching any reefs outside of the cluster
        n_touching = sum(distance_matrix[ind, Not(ind)] .== 0.0)

        return n_touching
    else
        return 0.0
    end
end

"""
    distance_to_next_reef(
        cluster::Int64,
        clusters::Vector{Int64},
        distance_matrix::Matrix{Float64}
    )

Find the maximum neighbouring reef distance that occurs within the target cluster.

# Arguments
- `cluster` : Target cluster to investigate.
- `clusters` : Vector containing all clusters assigned to each reef.
- `distance_matrix` : n*n matrix containing distances between each reef polygon.

# Returns
The largest neighbouring-reef-distance in the target `cluster` in km units.
"""
function distance_to_next_reef(
    cluster::Int64,
    clusters::Vector{Int64},
    distance_matrix::Matrix{Float64}
)
    ind = clusters .== cluster
    distances = distance_matrix[ind, ind]
    n_clusters = size(distances, 1)
    if n_clusters > 1
        reef_distances = maximum([minimum(distances[i, Not(i)]) for i in 1:n_clusters])
    else
        reef_distances = 1.0
    end
    return reef_distances
end
