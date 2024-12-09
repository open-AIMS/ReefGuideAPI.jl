"""
Clustering approach that includes functions for preprocessing data for clustering, and assigning optimal clusters.
* Note: All scoring functions have been standardised to output metrics in km/sqkm units
"""

using Random
using BlackBoxOptim, Distances, Clustering
using DataFrames

include("cluster_scoring.jl")

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
All scoring functions have been standardised to kilometres (sqkm; sqm / 1e6).
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
function GeoInterface_to_ArchGDAL_polygon(polygon::GeoInterface.Wrappers.Polygon)
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
@inline function distance_matrix!(reef_geoms::Vector{AG.IGeometry}, dist::Matrix{Float64})::Nothing
    return distance_matrix!(reef_geoms, dist, AG.distance)
end
@inline function distance_matrix!(reef_geoms::Vector{GIWrap.WrapperGeometry{false, false, T, Nothing} where T}, dist::Matrix{Float64})
    return distance_matrix!(reef_geoms, dist, GI.distance)
end

# Functions to apply clustering and scoring to datasets

"""
    _define_clustering(clustering_type)

Define the clustering function to use in assessment based on available algorithms in Clustering.jl.
"""
function _define_clustering(clustering_type)
    if clustering_type == :kmeans
        return (clust_matrix, cluster_param) -> assignments(kmeans(clust_matrix, floor(Int64, cluster_param)))
    elseif clustering_type == :kmedoids
        return (clust_matrix, cluster_param) -> assignments(kmedoids(clust_matrix, floor(Int64, cluster_param)))
    elseif clustering_type == :dbscan
        return (clust_matrix, cluster_param) -> assignments(dbscan(clust_matrix, cluster_param; metric=nothing))
    elseif clustering_type == :hclust
        return (clust_matrix, cluster_param) -> cutree(hclust(clust_matrix); k=cluster_param)
    end

    throw("Clustering type not implemented. Please implement in `clustering_function()")
end

"""
    optimisation_func(
        X,
        clustering_function,
        reef_df,
        clustering_matrix,
        max_side_bound,
        max_area_bound,
        max_distance_bound,
        max_benthic_bound,
        min_benthic_bound,
        min_reef_number
    )

Optimisation function for application of `clustering_function` method to `clustering_matrix`.
`X` is the parameter for optimisation with BlackBoxOptim.jl.

# Returns
Mean cluster score based on physical and ecological parameters combined with silhouette scores.
"""
function optimisation_func(
    X,
    clustering_function,
    reef_df,
    clustering_matrix,
    max_side_bound,
    max_area_bound,
    max_distance_bound,
    max_benthic_bound,
    min_benthic_bound,
    min_reef_number
)
    cluster_param = X[1]

        local clusters
    try
        clusters = clustering_function(clustering_matrix, cluster_param)
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

    sil_score = 1.0
    try
        sil_score = silhouettes(reef_df.clusters, clustering_matrix)
    catch err
        if !(err isa ArgumentError)
            rethrow(err)
        else
        # All locations assigned to a single cluster so assign worst score
        sil_score = 1.0
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
    cluster_type::Symbol,
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

    clustering_function = _define_clustering(cluster_type)
    if cluster_type == :dbscan
        dist_bnds = [(0.001, 1000)]
    else
        dist_bnds = [(min_clusters, n_locs)]
    end

    opt_func = x -> optimisation_func(
            x,
            clustering_function,
            df,
            clustering_matrix,
            max_side_bound,
            max_area_bound,
            max_distance_bound,
            max_benthic_bound,
            min_benthic_bound,
            min_reef_number
    )
    res = bboptimize(
        opt_func;
        SearchRange=dist_bnds,
        MaxSteps=n_steps,
        ϵ=epsilon,
    )
    best_params = best_candidate(res)

    if cluster_type == :dbscan
        return clustering_function(clustering_matrix, best_params[1])
    end

    return clustering_function(clustering_matrix, floor(Int64, best_params[1]))
end
