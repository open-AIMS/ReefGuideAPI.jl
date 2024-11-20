"""
Clustering approach that includes functions for scoring and assigning optimal clusters.
"""

import GeometryOps as GO
import GeoInterface as GI
import GeoDataFrames as GDF

using Clustering, Random, BlackBoxOptim
using Distances, Statistics, DataFrames


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
    function score_clusters(
        cluster_dataframe,
        max_side_bound,
        max_area_bound,
        max_distance_bound,
        max_benthic_bound,
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
- `max_benthic_bound` : Maximum area of sutiable benthic habitat (Units=squared km).

# Returns
DataFrame with a score for each cluster. Scores == 0.0 mean a cluster satisfies all criteria,
scores > 0.0 mean a cluster does not satisfy all criteria, with higher scores meaning values
are increasing, away from the criteria bounds.
"""
function score_clusters(
    cluster_dataframe,
    max_side_bound,
    max_area_bound,
    max_distance_bound,
    max_benthic_bound,
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

        if size(groupdf,1) == 1
            highest_reef_distance = max_distance_bound + 1
        else
            highest_reef_distance = maximum(distances_between_reefs(groupdf.geometry))
        end

        benthic_area = sum(groupdf.benthic_area)

        cluster_score = proportionate_score(side_length, max_side_bound) +
            proportionate_score(total_rectangular_area, max_area_bound) +
            proportionate_score(highest_reef_distance, max_distance_bound) +
            proportionate_score(benthic_area, max_benthic_bound) +
            proportionate_score(size(groupdf, 1), min_reef_number; rev=true)

        cluster_scores[i, :clusters] = first(unique(groupdf.clusters))
        cluster_scores[i, :score] = cluster_score

        if output_polygons
            cluster_scores[i, :cluster_poly] = cluster_polygons(groupdf)
        end
    end

    return cluster_scores
end

"""
    cluster_side_distance(groupdf)

Finds the length of the longest side of a cluster (assessing a cluster via rectangular shape).
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
    distances_between_reefs(geometries)

Returns the distance from each reef to its nearest neighbouring reef.
"""
function distances_between_reefs(geometries)
    centroids = GO.centroid.(geometries)
    adjacent_distances = []
    for centroid_j in centroids
        distances_j = euclidean.([centroid_j], centroids)
        push!(adjacent_distances, minimum(distances_j[distances_j .!= 0.0]))
    end

    return vcat(adjacent_distances...) / 1000
end

# Functions to apply clustering and scoring to datasets

"""
    opt_kmeans(
        X,
        reef_df,
        clustering_cols,
        max_side_bound,
        max_area_bound,
        max_distance_bound,
        max_benthic_bound,
        min_reef_number
    )

Objective function that applies `X` number of clusters to the `clustering_cols` columns from
`reef_df`. These clusters are scored based on their side distance, area, distance between reefs,
benthic area and number of reefs. Clusters are also scored based on sillhouette scores (reef cluster similarity).
Lower scores are preferred by BlackBoxOptim.
"""
function opt_kmeans(
    X,
    reef_df,
    clustering_cols,
    max_side_bound,
    max_area_bound,
    max_distance_bound,
    max_benthic_bound,
    min_reef_number
)
    n_clusters = X[1]

        local clusters
    try
        clusters = kmeans(Matrix(reef_df[:, clustering_cols])', floor(Int64, n_clusters), display=:none)
        if !clusters.converged
            # Return worst score if k-means has not converged
            return 1.0
        end
    catch err
        if err isa BoundsError
            return 1.0
        else
            rethrow(err)
        end
    end

    reef_df.clusters = clusters.assignments

    sil_score = -1.0
    try
        sil_score = silhouettes(reef_df.clusters, Matrix(reef_df[:, clustering_cols])'; metric=Euclidean())
    catch err
        if !(err isa ArgumentError)
            rethrow(err)
        else
        # All locations assigned to a single cluster so assign worst score
        sil_score = -1.0
        end
    end

    if size(sil_score, 1) > 1
        sil_score = DataFrame(clusters = reef_df.clusters, sil_score = sil_score)
        sil_score = DataFrames.combine(groupby(sil_score, :clusters), :sil_score => mean)

        cluster_criteria_scores = score_clusters(
            reef_df,
            max_side_bound,
            max_area_bound,
            max_distance_bound,
            max_benthic_bound,
            min_reef_number
        )
        sil_score = normalise(cluster_criteria_scores.score, (0,1)) .+ normalise(-sil_score.sil_score_mean, (0,1))
        #sil_score = normalise(sil_score, (-1, 1))
    end

    # Optimization direction is toward the minimum, so invert score.
    return mean(sil_score)
end

"""
    normalise(x, (a, b))

Normalise the vector `x` so that it's minimum value is `a` and its maximum value is `b`.

# Examples
- `normalise([1, 5, 10, 78] (0,1))` to return a vector with min=0 and max=1.
- `normalise([1, 5, 10, 78] (-1,1))` to return a vector with min=-1 and max=1.
"""
function normalise(x, (a, b))
    x_norm = (b - a) .* ((x .- minimum(x)) ./ (maximum(x) .- minimum(x))) .+ a
    return x_norm
end

"""
    cluster(
        df,
        clustering_cols,
        max_side_bound,
        max_area_bound,
        max_distance_bound,
        max_benthic_bound,
        min_reef_number;
        n_steps=100,
        epsilon=0.4,
        min_clusters=1
    )
"""
function cluster(
    df,
    clustering_cols,
    max_side_bound,
    max_area_bound,
    max_distance_bound,
    max_benthic_bound,
    min_reef_number;
    n_steps=100,
    epsilon=0.4,
    min_clusters=1
)
    Random.seed!(101)

    n_locs = length(unique(df.UNIQUE_ID))
    dist_bnds = [(min_clusters, n_locs)]

    opt_func = x -> opt_kmeans(
        x,
        df,
        clustering_cols,
        max_side_bound,
        max_area_bound,
        max_distance_bound,
        max_benthic_bound,
        min_reef_number
    )

    res = bboptimize(
        opt_func;
        SearchRange=dist_bnds,
        MaxSteps=n_steps,
        ϵ=epsilon,
    )
    best_params = best_candidate(res)

    clusters = kmeans(Matrix(df[:, clustering_cols])', floor(Int64, best_params[1]))

    assignments = clusters.assignments
    sil_score = silhouettes(assignments, Matrix(df[:, clustering_cols])'; metric=Euclidean())

    df.silhouette_score = sil_score
    df.clusters = assignments
    return df
end

# Functions for preprocessing data before clustering is applied.
function reef_proportion_in_criteria(x)
    if (length(collect(x)) > 0)# & (sum(x) > 0)
        return sum(collect(x) .> 0) / length(collect(x))
    else
        return 0.0
    end
end

"""
    reef_benthic_areas(geometries, bool_raster)

Calculates the area of each reef in `geometries` that overlaps suitable habitat (Ones) in `bool_raster`.
"""
function reef_benthic_areas(geometry, bool_raster)
    res = abs(step(bool_raster.dims[1]))
    benthic_areas = zonal(sum, bool_raster; of=geometry)

    return (benthic_areas .* (res*res)) / 1e6
end
