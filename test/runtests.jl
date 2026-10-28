using ClusterStats
using Test
using Statistics
using StatsBase

@testset "ClusterStats.jl" begin

    @testset "signed_rank_z" begin
        # Test basic functionality with known values
        @testset "Simple positive differences" begin
            x = [1.0, 2.0, 3.0, 4.0, 5.0]
            z = signed_rank_z(x)
            @test !isnan(z)
            @test z > 0  # All positive values should give positive z
        end

        @testset "Simple negative differences" begin
            x = [-1.0, -2.0, -3.0, -4.0, -5.0]
            z = signed_rank_z(x)
            @test !isnan(z)
            @test z < 0  # All negative values should give negative z
        end

        @testset "Zero differences" begin
            x = [0.0, 0.0, 0.0]
            z = signed_rank_z(x)
            @test isnan(z)  # Should return NaN for all zeros
        end

        @testset "Mixed with zeros" begin
            x = [0.0, 1.0, 2.0, 0.0, 3.0]
            z = signed_rank_z(x)
            @test !isnan(z)  # Should ignore zeros and compute on non-zero values
            @test z > 0
        end

        @testset "Symmetric differences" begin
            x = [-2.0, -1.0, 1.0, 2.0]
            z = signed_rank_z(x)
            @test abs(z) < 0.5  # Should be close to 0 for symmetric distribution
        end

        @testset "Known value" begin
            # n=5, W=15, mean=7.5, var=5*6*11/24
            @test signed_rank_z([1.0, 2.0, 3.0, 4.0, 5.0]) ≈ 7.5 / sqrt(13.75)
        end

        @testset "Ties use average ranks (order independent)" begin
            x = [-2.0, -1.0, 1.0, 2.0]
            @test isapprox(signed_rank_z(x), 0.0; atol=1e-12)
            @test signed_rank_z(x) ≈ signed_rank_z(reverse(x))
            y = [1.0, -1.0, 2.0, -2.0, 3.0, 1.0]
            @test signed_rank_z(y) ≈ signed_rank_z(y[[6, 3, 1, 5, 2, 4]])
        end

        @testset "Custom comparison value" begin
            x = [3.0, 4.0, 5.0, 6.0]
            y = 2.0
            z = signed_rank_z(x, y)
            @test !isnan(z)
            @test z > 0  # All values > y should give positive z
        end
    end

    @testset "region_props" begin
        using Images

        @testset "Single region" begin
            img = [1.0 2.0; 3.0 4.0]
            labels = [1 1; 1 1]
            regions, values = region_props(img, labels)

            @test length(regions) == 1
            @test length(values) == 1
            @test length(regions[1]) == 4
            @test values[1] == [1.0, 3.0, 2.0, 4.0]
        end

        @testset "Multiple regions" begin
            img = [1.0 2.0 3.0; 4.0 5.0 6.0; 7.0 8.0 9.0]
            labels = [1 1 2; 1 0 2; 0 0 2]
            regions, values = region_props(img, labels)

            @test length(regions) == 2
            @test length(values) == 2
            @test length(regions[1]) == 3  # First region has 3 pixels
            @test length(regions[2]) == 3  # Second region has 3 pixels
        end

        @testset "No regions (all background)" begin
            img = [1.0 2.0; 3.0 4.0]
            labels = [0 0; 0 0]
            regions, values = region_props(img, labels)

            @test length(regions) == 0
            @test length(values) == 0
        end

        @testset "Region pixel values" begin
            img = [10.0 20.0; 30.0 40.0]
            labels = [1 2; 1 2]
            regions, values = region_props(img, labels)

            @test mean(values[1]) == 20.0  # (10 + 30) / 2
            @test mean(values[2]) == 30.0  # (20 + 40) / 2
        end
    end

    @testset "cluster_correct" begin

        @testset "Basic functionality - no clusters" begin
            # Create data with no significant clusters (already as statistics)
            obs = randn(5, 5)  # statistic map, 5x5 grid
            null = randn(5, 5, 100)  # 100 permutation statistic maps

            result, zmap, thresh_zmap = cluster_correct(obs, null, zval=3.0, pval=0.05)

            @test size(result) == (5, 5)
            @test size(zmap) == (5, 5)
            @test size(thresh_zmap) == (5, 5)
            @test all(isfinite.(result))
        end

        @testset "Compute t-statistics" begin
            obs = randn(20, 4, 4)
            null = randn(20, 4, 4, 50)

            result, zmap, thresh_zmap = cluster_correct(obs, null,
                                                        compute_stat=true,
                                                        test_type="t",
                                                        pval=0.05)

            @test size(result) == (4, 4)
            @test all(isfinite.(result))
        end

        @testset "Different tail options" begin
            obs = randn(3, 3)  # statistic map
            null = randn(3, 3, 30)  # permutation statistic maps

            # Two-tailed (default)
            result_two, _, _ = cluster_correct(obs, null, tail=0, pval=0.1)
            @test size(result_two) == (3, 3)

            # Right-tailed
            result_right, _, _ = cluster_correct(obs, null, tail=1, pval=0.1)
            @test size(result_right) == (3, 3)

            # Left-tailed
            result_left, _, _ = cluster_correct(obs, null, tail=-1, pval=0.1)
            @test size(result_left) == (3, 3)
        end

        @testset "Different cluster statistics" begin
            obs = randn(4, 4)  # statistic map
            null = randn(4, 4, 30)  # permutation statistic maps

            # maxsum
            result_sum, _, _ = cluster_correct(obs, null, cluster_stat="maxsum", pval=0.1)
            @test size(result_sum) == (4, 4)

            # maxsize
            result_size, _, _ = cluster_correct(obs, null, cluster_stat="maxsize", pval=0.1)
            @test size(result_size) == (4, 4)
        end

        @testset "Minimum cluster size filtering" begin
            obs = randn(6, 6)  # statistic map
            null = randn(6, 6, 30)  # permutation statistic maps

            result, _, _ = cluster_correct(obs, null, min_size=3, pval=0.1)
            @test size(result) == (6, 6)
        end

        @testset "Power transformation" begin
            obs = randn(4, 4)  # statistic map
            null = randn(4, 4, 30)  # permutation statistic maps

            result, _, _ = cluster_correct(obs, null, power=2, pval=0.1)
            @test size(result) == (4, 4)
            @test all(isfinite.(result))
        end

        @testset "Signed rank test" begin
            obs = randn(15, 3, 3)
            null = randn(15, 3, 3, 20)

            result, zmap, thresh_zmap = cluster_correct(obs, null,
                                                        compute_stat=true,
                                                        test_type="z",
                                                        pval=0.1)

            @test size(result) == (3, 3)
            @test all(isfinite.(result) .|| isnan.(result))  # Allow NaN for zero variance
        end

        @testset "thresh_zmap is not aliased to the corrected result" begin
            base = [0.01 * i * j for i in 1:6, j in 1:6]
            null = cat([base .+ (k - 50.5) / 30 for k in 1:100]...; dims=3)
            obs = copy(base)
            obs[1:2, 1:2] .+= 10.0   # 4-voxel cluster, survives
            obs[6, 6] += 10.0        # isolated voxel, removed by min_size

            result, zmap, thresh_zmap = cluster_correct(obs, null, zval=1.96, pval=0.05, min_size=2)

            @test all(result[1:2, 1:2] .!= 0)
            @test result[6, 6] == 0
            @test thresh_zmap[6, 6] != 0   # pre-correction map keeps the removed voxel
            @test thresh_zmap != result
        end

        @testset "Non-finite statistics are not labelled as clusters" begin
            # voxel (1,1) is all zeros in obs and null -> NaN signed-rank z
            obs = reshape(sin.(1.0:135.0), 15, 3, 3)
            null = reshape(cos.(1.0:(15 * 3 * 3 * 20)), 15, 3, 3, 20)
            obs[:, 1, 1] .= 0.0
            null[:, 1, 1, :] .= 0.0

            result, zmap, thresh_zmap = cluster_correct(obs, null, compute_stat=true,
                                                        test_type="z", pval=0.1)
            @test !any(isnan, result)
            @test !any(isnan, zmap)
            @test !any(isnan, thresh_zmap)
            @test result[1, 1] == 0

            # NaN in a precomputed statistic map
            stat = [0.1 * (i + j) for i in 1:4, j in 1:4]
            stat[2, 2] = NaN
            nullmaps = cat([[0.1 * (i + j) + (k - 10.5) / 10 for i in 1:4, j in 1:4] for k in 1:20]...; dims=3)
            result, zmap, _ = cluster_correct(stat, nullmaps, pval=0.1)
            @test !any(isnan, result)
            @test !any(isnan, zmap)
        end

        @testset "Invalid options throw" begin
            obs = zeros(3, 3)
            null = reshape(collect(1.0:90.0), 3, 3, 10)
            @test_throws ArgumentError cluster_correct(obs, null, tail=2)
            @test_throws ArgumentError cluster_correct(obs, null, cluster_stat="maxsiz")
        end

        @testset "Strong signal detection" begin
            # Create data with a clear cluster
            obs = zeros(20, 8, 8)
            # Add strong signal in center region
            obs[:, 3:6, 3:6] .+= 5.0 .+ randn(20, 4, 4) * 0.5

            null = randn(20, 8, 8, 50) * 0.5

            result, zmap, thresh_zmap = cluster_correct(obs, null,
                                                        compute_stat=true,
                                                        test_type="t",
                                                        zval=1.96,
                                                        pval=0.05,
                                                        min_size=2)

            # Should detect cluster in the center region
            @test sum(result .!= 0) > 0  # At least some significant voxels
            @test maximum(abs.(result)) > 0
        end
    end
end
