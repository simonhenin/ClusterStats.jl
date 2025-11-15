using Statistics, StatsBase
using Images
using LinearAlgebra

function cluster_correct(obs, null; zval=1.96, pval=0.05, tail=0, cluster_stat="maxsum", min_size=1, compute_stat=false, test_type="t", power=1)
# When compute_stat=true:
#   obs = nsubj x dim1 x dim2 (raw subject data)
#   null = nsubj x dim1 x dim2 x n_perms (raw permuted subject data)
# When compute_stat=false:
#   obs = dim1 x dim2 (already computed statistic map)
#   null = dim1 x dim2 x n_perms (already computed permutation statistic maps)

    if compute_stat
        n, dim1, dim2, n_perms = size(null)
        tmap_null = zeros(dim1, dim2, n_perms)
        tmap = zeros(dim1, dim2)
        if test_type == "t"
            println("Using t-statistic")
            obs_std = std(obs, dims=1)[1,:,:]
            null_std = std(null, dims=1)[1, :, :, :]

            # Check for zero standard deviations
            if any(obs_std .== 0)
                @warn "Zero standard deviation detected in observed data. Setting those t-values to 0."
                obs_std[obs_std .== 0] .= Inf  # Will result in t=0
            end
            if any(null_std .== 0)
                @warn "Zero standard deviation detected in null data. Setting those t-values to 0."
                null_std[null_std .== 0] .= Inf  # Will result in t=0
            end

            tmap = mean(obs, dims=1)[1,:,:] ./ obs_std .* sqrt(n)
            tmap_null = mean(null, dims=1)[1, :, :, :] ./ null_std .* sqrt(n)
        elseif test_type == "z"
            println("Using signed rank test statistic")
            for j in 1:dim1
                for k in 1:dim2
                    tmap[j, k] = signed_rank_z(obs[:, j, k])
                    for r in 1:n_perms
                        tmap_null[j, k, r] = signed_rank_z(null[:, j, k, r])
                    end
                end
            end
        else
            error("test_type must be either 't' or 'z'")
        end
        permmaps = permutedims(tmap_null, (3, 1, 2)) # put perms on the first dimension (now: n_perms x dim1 x dim2)
        diffmap = tmap
    else
        # obs is already a statistic map (dim1, dim2), not raw subject data
        diffmap = obs
        # null is already statistic maps (dim1, dim2, n_perms)
        permmaps = permutedims(null, (3, 1, 2)) # put perms on the first dimension (now: n_perms x dim1 x dim2)
    end

    if power > 1
        diffmap = sign.(diffmap) .* (diffmap .^ power)
        permmaps = sign.(permmaps) .* (permmaps .^ power)
    end


    n_permutes = size(permmaps, 1)

    mean_h = mean(permmaps, dims=1)[1, :, :]
    std_h = std(permmaps, dims=1)[1, :, :]


    zmap = (diffmap .- mean_h) ./ std_h

    thresh_zmap = copy(zmap)
    thresh_zmap[abs.(thresh_zmap) .< zval] .= 0

    cluster_correct_zmap = thresh_zmap
    labels = label_components(cluster_correct_zmap .!= 0)
    islands, islandvals = region_props(cluster_correct_zmap, labels)
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
            # maxsize cluster statistic
            labels = label_components(positive_thresh .!= 0)
            if sum(labels) > 0
                clusters, _ = region_props(positive_thresh, labels)
                if !isempty(clusters)
                    sizes = [length(indices) for (key, indices) in clusters]
                    positive_max_cluster_sizes[permi] = maximum(sizes)
                end
            end

            labels = label_components(negative_thresh .!= 0)
            if sum(labels) > 0
                clusters, _ = region_props(negative_thresh, labels)
                if !isempty(clusters)
                    sizes = [length(indices) for (key, indices) in clusters]
                    negative_max_cluster_sizes[permi] = maximum(sizes)
                end
            end
        end
    end

    if tail == 0
        positive_cluster_thresh = percentile(positive_max_cluster_sizes, 100 - (100 * pval / 2))
        if cluster_stat == "maxsum"
            negative_cluster_thresh = percentile(negative_max_cluster_sizes, (100 * pval / 2))
        else
            negative_cluster_thresh = percentile(negative_max_cluster_sizes, 100 - (100 * pval / 2))
        end
    elseif tail == 1
        positive_cluster_thresh = percentile(positive_max_cluster_sizes, 100 - (100 * pval))
        if cluster_stat == "maxsum"
            negative_cluster_thresh = -Inf
        else
            negative_cluster_thresh = Inf
        end
    else
        positive_cluster_thresh = Inf
        if cluster_stat == "maxsum"
            negative_cluster_thresh = percentile(negative_max_cluster_sizes, (100 * pval))
        else
            negative_cluster_thresh = percentile(negative_max_cluster_sizes, 100-(100 * pval))
        end
    end
    println("Positive cluster threshold: $(positive_cluster_thresh)")
    println("Negative cluster threshold: $(negative_cluster_thresh)")

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
    all_labels = sort(unique(labels))
    # Skip the background (0) if it exists
    unique_labels = filter(x -> x != 0, all_labels)
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


function signed_rank_z(x, y=0)
    diffs = x .- y
    nonzero_diffs = diffs[diffs .!= 0]  # Remove zero differences
    n = length(nonzero_diffs)

    if n == 0
        return NaN  # No valid data for test
    end

    ranks = ordinalrank(abs.(nonzero_diffs))  # Compute ranks
    signed_ranks = sign.(nonzero_diffs) .* ranks  # Apply signs

    W = sum(signed_ranks[signed_ranks .> 0])  # Sum of positive ranks

    # Compute z-score
    mean_W = n * (n + 1) / 4
    std_W = sqrt(n * (n + 1) * (2*n + 1) / 24)
    z_stat = (W - mean_W) / std_W

    return z_stat
end