using Statistics, StatsBase
using Images
using LinearAlgebra

function cluster_correct(obs, null; zval=1.96, pval=0.05, tail=0, cluster_stat="maxsum", min_size=1, compute_stat=false)

    if compute_stat
        pmap = []
        tmap = []
        tmap_null = []

        n, jj, kk, rr = size(null)
        tmap_null = zeros(jj, kk, rr)
        tmap = zeros(jj, kk)
        
        tmap = mean(obs, dims=1)[1,:,:] ./ std(obs, dims=1)[1,:,:] .* sqrt(n)
        tmap_null = mean(null, dims=1)[1, :, :, :] ./ std(null, dims=1)[1, :, :, :] .* sqrt(n)

        tmap = sign.(tmap) .* (tmap.^2)
        tmap_null = sign.(tmap_null) .* (tmap_null.^2)
        # for j in 1:jj
        #     for k in 1:kk
        #         tmap[j, k] = mean(obs[:, j, k]) / std(obs[:, j, k]) * sqrt(n)
        #         for r in 1:rr
        #             tmap_null[j, k, r] = mean(null[:, j, k, r]) / std(null[:, j, k, r]) * sqrt(n)
        #         end
        #     end
        # end
        permmaps = permutedims(tmap_null, (3, 1, 2)) # put perms on the first dimension
        diffmap = tmap
    else
        diffmap = obs
        permmaps = permutedims(null, (3, 1, 2)) # put perms on the first dimension
    end

    n_permutes = size(permmaps, 1)

    mean_h = mean(permmaps, dims=1)[1, :, :]
    std_h = std(permmaps, dims=1)[1, :, :]


    zmap = (diffmap .- mean_h) ./ std_h

    thresh_zmap = copy(zmap)
    thresh_zmap[abs.(thresh_zmap) .< zval] .= 0

    cluster_correct_zmap = thresh_zmap
    labels = label_components(cluster_correct_zmap .!= 0)
    islands, islandvals = cluster_stats.region_props(cluster_correct_zmap, labels)
    too_small = [key for (key, val) in islands if length(val) < min_size]

    if !isempty(too_small)
        for i in too_small
            cluster_correct_zmap[islands[i]] .= 0
        end
        foreach(key -> delete!(islands, key), too_small)
    end

    if isempty(islands)
        println("No clusters found in data")
        return zeros(size(zmap)), zmap, thresh_zmap
    end

    cluster_vals = [mean(val) for (key, val) in islandvals]
    pos_clusters = length(findall(x->(x>0), cluster_vals))
    neg_clusters = length(findall(x->(x<0), cluster_vals))
    println("Found $(pos_clusters) positive clusters")
    println("Found $(neg_clusters) negative clusters")

    positive_max_cluster_sizes = zeros(n_permutes)
    negative_max_cluster_sizes = zeros(n_permutes)

    for permi in 1:n_permutes
        threshimg = permmaps[permi, :, :]
        threshimg = (threshimg .- mean_h) ./ std_h

        positive_thresh = copy(threshimg)
        # positive_thresh[ collect((positive_thresh .> 0.0) .& (abs.(positive_thresh .- mean_h) .< zval)) ] .= 0.0
        positive_thresh[ collect((positive_thresh .> 0.0) .& (abs.(positive_thresh) .< zval)) ] .= 0.0
        positive_thresh[ collect(positive_thresh .< 0.0) ] .= 0.0

        negative_thresh = copy(threshimg)
        # negative_thresh[ collect( (negative_thresh .< 0.0) .& (abs.(negative_thresh .- mean_h) .< zval))] .= 0.0
        negative_thresh[ collect( (negative_thresh .< 0.0) .& (abs.(negative_thresh) .< zval))] .= 0.0
        negative_thresh[ collect(negative_thresh .> 0.0) ] .= 0.0

        if cluster_stat == "maxsum"
            labels = label_components(positive_thresh .!= 0)
            if sum(labels) > 0
                clusters, cluster_vals = region_props(positive_thresh, labels)
                if !isempty(clusters)
                    sums = [sum(vec) for (key, vec) in cluster_vals]
                    positive_max_cluster_sizes[permi] = maximum(sums)
                end
            end

            labels = label_components(negative_thresh .!= 0)
            if sum(labels) > 0
                clusters, cluster_vals = region_props(negative_thresh, labels)
                if !isempty(clusters)
                    sums = [sum(vec) for (key, vec) in cluster_vals]
                    negative_max_cluster_sizes[permi] = minimum(sums)
                end
            end
        else
            # TODO: implement maxsize
            islands_ = label_components(positive_thresh .!= 0)
            if !isempty(islands_)
                tempclustsizes = map(length, islands_)
                positive_max_cluster_sizes[permi] = maximum(tempclustsizes)
            end

            islands_ = label_components(negative_thresh .!= 0)
            if !isempty(islands_)
                tempclustsizes = map(length, islands_)
                negative_max_cluster_sizes[permi] = maximum(tempclustsizes)
            end
        end
    end

    if tail == 0
        positive_cluster_thresh = percentile(positive_max_cluster_sizes, 100 - (100 * pval / 2))
        negative_cluster_thresh = percentile(negative_max_cluster_sizes, 100 - (100 * pval / 2))
    elseif tail == 1
        positive_cluster_thresh = percentile(positive_max_cluster_sizes, 100 - (100 * pval))
        if cluster_stat == "maxsum"
            negative_cluster_thresh = -Inf
        else
            negative_cluster_thresh = Inf
        end
    else
        positive_cluster_thresh = Inf
        negative_cluster_thresh = percentile(negative_max_cluster_sizes, 100 - (100 * pval))
    end

    for i in eachindex(islands)
        island_size = length(islands[i])
        island_mean_value = mean(zmap[islands[i]])
        island_sum = sum(zmap[islands[i]])

        if cluster_stat == "maxsum"
            if island_mean_value > 0 && (island_sum < positive_cluster_thresh)
                cluster_correct_zmap[islands[i]] .= 0
            end
            if island_mean_value < 0 && (island_sum > negative_cluster_thresh)
                cluster_correct_zmap[islands[i]] .= 0
            end
        else
            if island_mean_value > 0 && (island_size < positive_cluster_thresh)
                cluster_correct_zmap[islands[i]] .= 0
            end
            if island_mean_value < 0 && (island_size < negative_cluster_thresh)
                cluster_correct_zmap[islands[i]] .= 0
            end
        end
    end

    return cluster_correct_zmap, zmap, thresh_zmap
end

function region_props(img::Matrix{Float64}, labels::Matrix{Int64})
    unique_labels = sort(unique(labels))[2:end]  # Skip the background (0)
    regions = Dict()
    values = Dict()
    if (length(unique_labels) > 0)
        for label in unique_labels
            # Find the indices of the current component
            component_indices = findall(x -> x == label, labels)
            
            # Get pixel values of the component from the original image
            pixel_values = [img[i] for i in component_indices]
            
            # Store pixel values in the dictionary
            regions[label] = component_indices
            values[label] = pixel_values
        end
    end
    return regions, values
end