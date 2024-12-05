"""
Clustering approach that includes functions for scoring and assigning optimal clusters.
"""

using Random
using Statistics

import ArchGDAL as AG
import GeometryOps as GO
import GeoInterface as GI

using BlackBoxOptim, Distances, Clustering

using DataFrames
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
        min_benthic_bound,
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
    clustering_matrix,
    max_side_bound,
    max_area_bound,
    max_distance_bound,
    max_benthic_bound,
    min_benthic_bound,
    min_reef_number
)
    n_clusters = X[1]

        local clusters
    try
        clusters = kmeans(clustering_matrix, floor(Int64, n_clusters), display=:none)
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
        sil_score = silhouettes(reef_df.clusters, clustering_matrix)
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
            clustering_matrix,
            max_side_bound,
            max_area_bound,
            max_distance_bound,
            max_benthic_bound,
            min_benthic_bound,
            min_reef_number
        )
        sil_score = normalise(cluster_criteria_scores.score, (0,1)) .+ normalise(-sil_score.sil_score_mean, (0,1))
        #sil_score = normalise(sil_score, (-1, 1))
    end

    # Optimization direction is toward the minimum, so invert score.
    return mean(sil_score)
end

"""
    opt_kmedoids(
        X,
        reef_df,
        clustering_matrix::Matrix,
        max_side_bound,
        max_area_bound,
        max_distance_bound,
        max_benthic_bound,
        min_benthic_bound,
        min_reef_number
    )

Objective function that applies `X` number of clusters to the `clustering_matrix` from
`reef_df`. These clusters are scored based on their side distance, area, distance between reefs,
benthic area and number of reefs. Clusters are also scored based on sillhouette scores (reef cluster similarity).
Lower scores are preferred by BlackBoxOptim.
"""
function opt_kmedoids(
    X,
    reef_df,
    clustering_matrix::Matrix,
    max_side_bound,
    max_area_bound,
    max_distance_bound,
    max_benthic_bound,
    min_benthic_bound,
    min_reef_number
)
    n_clusters = X[1]

        local clusters
    try
        clusters = kmedoids(clustering_matrix, floor(Int64, n_clusters), display=:none)
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
        sil_score = silhouettes(reef_df.clusters, clustering_matrix)
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
            clustering_matrix,
            max_side_bound,
            max_area_bound,
            max_distance_bound,
            max_benthic_bound,
            min_benthic_bound,
            min_reef_number
        )
        sil_score = normalise(cluster_criteria_scores.score, (0,1)) .+ normalise(-sil_score.sil_score_mean, (0,1))
        #sil_score = normalise(sil_score, (-1, 1))
    end

    # Optimization direction is toward the minimum, so invert score.
    return mean(sil_score)
end

"""
    opt_kmedoids(
        X,
        reef_df,
        clustering_matrix::Matrix,
        max_side_bound,
        max_area_bound,
        max_distance_bound,
        max_benthic_bound,
        min_benthic_bound,
        min_reef_number
    )

Objective function that applies `X` number of clusters to the `clustering_matrix` from
`reef_df`. These clusters are scored based on their side distance, area, distance between reefs,
benthic area and number of reefs. Clusters are also scored based on sillhouette scores (reef cluster similarity).
Lower scores are preferred by BlackBoxOptim.
"""
function opt_dbscan(
    X,
    reef_df,
    clustering_matrix::Matrix,
    max_side_bound,
    max_area_bound,
    max_distance_bound,
    max_benthic_bound,
    min_benthic_bound,
    min_reef_number
)
    n_clusters = X[1]

        local clusters
    try
        clusters = dbscan(clustering_matrix, n_clusters; metric=nothing)
        # if !clusters.converged
        #     # Return worst score if k-means has not converged
        #     return 1.0
        # end
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
        sil_score = silhouettes(reef_df.clusters, clustering_matrix)
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
            clustering_matrix,
            max_side_bound,
            max_area_bound,
            max_distance_bound,
            max_benthic_bound,
            min_benthic_bound,
            min_reef_number
        )
        sil_score = normalise(cluster_criteria_scores.score, (0,1)) .+ normalise(-sil_score.sil_score_mean, (0,1))
        #sil_score = normalise(sil_score, (-1, 1))
    end

    # Optimization direction is toward the minimum, so invert score.
    return mean(sil_score)
end

"""
    opt_kmedoids(
        X,
        reef_df,
        clustering_matrix::Matrix,
        max_side_bound,
        max_area_bound,
        max_distance_bound,
        max_benthic_bound,
        min_benthic_bound,
        min_reef_number
    )

Objective function that applies `X` number of clusters to the `clustering_matrix` from
`reef_df`. These clusters are scored based on their side distance, area, distance between reefs,
benthic area and number of reefs. Clusters are also scored based on sillhouette scores (reef cluster similarity).
Lower scores are preferred by BlackBoxOptim.
"""
function opt_hclust(
    X,
    reef_df,
    clustering_matrix::Matrix,
    max_side_bound,
    max_area_bound,
    max_distance_bound,
    max_benthic_bound,
    min_benthic_bound,
    min_reef_number
)
    n_clusters = X[1]

        local clusters
    try
        clusters = cutree(hclust(clustering_matrix); k=floor(Int64, n_clusters))
        # if !clusters.converged
        #     # Return worst score if k-means has not converged
        #     return 1.0
        # end
    catch err
        if err isa BoundsError
            return 1.0
        else
            rethrow(err)
        end
    end

    reef_df.clusters = clusters

    sil_score = -1.0
    try
        sil_score = silhouettes(reef_df.clusters, clustering_matrix)
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
            clustering_matrix,
            max_side_bound,
            max_area_bound,
            max_distance_bound,
            max_benthic_bound,
            min_benthic_bound,
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
        min_benthic_bound,
        min_reef_number;
        n_steps=100,
        epsilon=0.4,
        min_clusters=1
    )

Finds the optimal clustering solution for the `clustering_cols` or `clustering_matrix` within
`n_steps`. Clusters are scored based on the specified criteria in the function. If `clustering_cols`
is used (n*d matrix) then kmeans() is used for clustering. If `clustering_matrix` is used (n*n distance matrix)
then kmedoids() is used for clustering.
"""
function cluster(
    df,
    clustering_matrix::Matrix,
    cluster_function::Symbol,
    max_side_bound,
    max_area_bound,
    max_distance_bound,
    max_benthic_bound,
    min_benthic_bound,
    min_reef_number;
    n_steps=100,
    epsilon=0.4,
    min_clusters=1
)
    Random.seed!(101)

    n_locs = length(unique(df.UNIQUE_ID))

    if cluster_function == :kmeans
        dist_bnds = [(min_clusters, n_locs)]
        opt_func = x -> opt_kmeans(
            x,
            df,
            clustering_matrix,
            max_side_bound,
            max_area_bound,
            max_distance_bound,
            max_benthic_bound,
            min_benthic_bound,
            min_reef_number
        )
    elseif cluster_function == :dbscan
        dist_bnds = [(0.001, 1000)]
        opt_func = x -> opt_dbscan(
            x,
            df,
            clustering_matrix,
            max_side_bound,
            max_area_bound,
            max_distance_bound,
            max_benthic_bound,
            min_benthic_bound,
            min_reef_number
        )
    elseif cluster_function == :kmedoids
        dist_bnds = [(min_clusters, n_locs)]
        opt_func = x -> opt_kmedoids(
            x,
            df,
            clustering_matrix,
            max_side_bound,
            max_area_bound,
            max_distance_bound,
            max_benthic_bound,
            min_benthic_bound,
            min_reef_number
        )
    elseif cluster_function == :hclust
        dist_bnds = [(min_clusters, n_locs)]
        opt_func = x -> opt_hclust(
            x,
            df,
            clustering_matrix,
            max_side_bound,
            max_area_bound,
            max_distance_bound,
            max_benthic_bound,
            min_benthic_bound,
            min_reef_number
        )
    else
        throw("$(cluster_function) clustering algorithm not implemented. Functions include :kmeans, :dbscan, :kmedoids, :hclust.")
    end

    res = bboptimize(
        opt_func;
        SearchRange=dist_bnds,
        MaxSteps=n_steps,
        ϵ=epsilon,
    )
    best_params = best_candidate(res)

    if cluster_function == :kmeans
        assignments = kmeans(clustering_matrix, floor(Int64, best_params[1])).assignments
    elseif cluster_function == :dbscan
        assignments = dbscan(clustering_matrix, best_params; metric=nothing).assignments
    elseif cluster_function == :kmedoids
        assignments = kmedoids(clustering_matrix, floor(Int64, best_params[1])).assignments
    elseif cluster_function == :hclust
        assignments = cutree(hclust(clustering_matrix); k=floor(Int64, best_params[1]))
    end

    sil_score = silhouettes(assignments, clustering_matrix)

    df.silhouette_score = sil_score
    df.clusters = assignments
    return df
end

# Functions for preprocessing data before clustering is applied.
function reef_proportion_in_criteria(x)
    if (length(collect(x)) > 0)# & (sum(x) > 0)
        return sum(collect(x) .> 0) / length(collect(x))
    end

    return 0.0
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

function GeoInterface_to_ArchGDAL_polygon(polygon::GeoInterface.Wrappers.MultiPolygon)
    polygon = GI.getgeom(polygon)
    ArchGDAL_polygon = ArchGDAL.createpolygon()
    for sub_polygon in polygon
        points = collect.(GO.Point.(GI.getpoint(sub_polygon)))
        ring = ArchGDAL.createlinearring(points)
        ArchGDAL.addgeom!(ArchGDAL_polygon, ring)
    end

    return ArchGDAL_polygon
end
function GeoInterface_to_ArchGDAL_polygon(polygon)
    n_holes = GI.nhole(polygon)
    if n_holes > 0 # Check if a polygon is a single polygon containing holes (> 1 ring)
        rings = n_holes + 1
        ArchGDAL_polygon = ArchGDAL.createpolygon()
        for subpolygon in 1:rings
            sub_ring = GI.getring(polygon, subpolygon)
            points = collect.(GO.Point.(GI.getpoint(sub_ring)))
            AG_ring = ArchGDAL.createlinearring(points)
            ArchGDAL.addgeom!(ArchGDAL_polygon, AG_ring)
        end

        return ArchGDAL_polygon
    end

    points = GO.Point.(GI.getpoint(polygon))
    points = collect.(points)

    return ArchGDAL.createpolygon(points)
end

"""
    distance_matrix!(reef_geoms::Vector, dist::Matrix{Float64}, func::Function)::Nothing

Create a distance matrix using the provided function.

# Arguments
- `reef_geoms` : Geometries to calculate distances between
- `dist` : distanec/adjacency matrix
- `func` : Distance calculation function
"""
function distance_matrix!(reef_geoms::Vector, dist::Matrix{Float64}, func::Function)::Nothing
    for ii in axes(dist, 1)
        Threads.@threads for jj in axes(dist, 2)
            if ii == jj || !iszero(dist[ii, jj])
                continue
            end

            @views dist[ii, jj] = func(reef_geoms[ii], reef_geoms[jj])
        end

        dist[:, ii] .= dist[ii, :]
    end

    return nothing
end

"""
    create_distance_matrix(reef_geoms::Vector)::Matrix

Calculate matrix of unique distances between reefs.

# Returns
Distance between reefs in meters
"""
function create_distance_matrix(reef_geoms::Vector)::Matrix{Float64}
    n_reefs = size(reef_geoms, 1)
    dist = zeros(n_reefs, n_reefs)
    distance_matrix!(reef_geoms, dist)

    return dist
end
@inline function distance_matrix!(reef_geoms::Vector{Tuple{Float64, Float64}}, dist::Matrix{Float64})::Nothing
    return distance_matrix!(reef_geoms, dist, haversine)
end
@inline function distance_matrix!(reef_geoms::Vector{AG.IGeometry{AG.wkbPolygon}}, dist::Matrix{Float64})::Nothing
    return distance_matrix!(reef_geoms, dist, AG.distance)
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
    return reef_distances ./ 1000
end
