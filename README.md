# ClusterStats

[![Build Status](https://github.com/simonhenin/ClusterStats.jl/workflows/CI/badge.svg)](https://github.com/simonhenin/ClusterStats.jl/actions)

A Julia package for cluster-based permutation testing on 2D data (e.g., time-frequency maps, brain imaging data). This package implements family-wise error rate (FWER) correction using cluster-based permutation statistics.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/simonhenin/ClusterStats.jl")
```

## Features

- **Cluster-based permutation testing** with FWER correction
- **Multiple test statistics**: t-test and Wilcoxon signed-rank test
- **Flexible cluster statistics**: maximum cluster sum or maximum cluster size
- **Two-tailed and one-tailed testing**
- **Minimum cluster size filtering**
- **Power transformation** for enhanced sensitivity

## Basic Usage

### Example 1: Pre-computed Statistics

If you already have statistic maps from your analysis:

```julia
using ClusterStats

# Your observed statistic map (e.g., t-values, z-values)
obs = randn(50, 50)  # 50x50 grid

# Permutation statistic maps (1000 permutations)
null = randn(50, 50, 1000)

# Run cluster correction
result, zmap, thresh_zmap = cluster_correct(obs, null,
                                            zval=1.96,     # z-threshold
                                            pval=0.05,     # p-value threshold
                                            tail=0)        # two-tailed test

# result contains cluster-corrected statistic map
# zmap contains z-scored statistic map
# thresh_zmap contains thresholded map before cluster correction
```

### Example 2: Computing Statistics from Raw Data

If you have raw subject data and want to compute t-statistics:

```julia
using ClusterStats

# Raw data: 20 subjects x 50 time points x 40 frequencies
obs = randn(20, 50, 40)
null = randn(20, 50, 40, 1000)  # 1000 permutations

# Compute t-statistics and perform cluster correction
result, zmap, thresh_zmap = cluster_correct(obs, null,
                                            compute_stat=true,
                                            test_type="t",
                                            zval=1.96,
                                            pval=0.05)
```

### Example 3: Non-parametric Test

For non-normal data, use the Wilcoxon signed-rank test:

```julia
result, zmap, thresh_zmap = cluster_correct(obs, null,
                                            compute_stat=true,
                                            test_type="z",  # Wilcoxon signed-rank
                                            zval=1.96,
                                            pval=0.05)
```

### Example 4: One-tailed Test with Cluster Size

```julia
result, zmap, thresh_zmap = cluster_correct(obs, null,
                                            cluster_stat="maxsize",  # use cluster size
                                            tail=1,                  # right-tailed
                                            min_size=5,             # minimum cluster size
                                            pval=0.05)
```

## Parameters

- `obs`: Observed data (2D for pre-computed stats, 3D for raw data)
- `null`: Null distribution data (3D for pre-computed stats, 4D for raw data)
- `zval`: Z-score threshold for initial thresholding (default: 1.96)
- `pval`: P-value threshold for cluster correction (default: 0.05)
- `tail`: Test direction: 0 (two-tailed), 1 (right-tailed), -1 (left-tailed)
- `cluster_stat`: Cluster statistic: "maxsum" or "maxsize"
- `min_size`: Minimum cluster size in pixels (default: 1)
- `compute_stat`: Whether to compute statistics from raw data (default: false)
- `test_type`: Statistical test: "t" or "z" (default: "t")
- `power`: Power transformation exponent (default: 1)

## Exported Functions

- `cluster_correct`: Main function for cluster-based permutation testing
- `region_props`: Extract properties of labeled regions
- `signed_rank_z`: Compute Wilcoxon signed-rank z-statistic

## License

This package is licensed under the MIT License:

```
MIT License

Copyright (c) 2025 Simon Henin

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

## Disclaimer

**USE AT YOUR OWN RISK**

This software is provided for research and educational purposes. While efforts have been made to ensure correctness, users are responsible for validating results and ensuring the methods are appropriate for their specific use case. The authors assume no liability for any errors, omissions, or consequences resulting from the use of this software.

Always verify statistical results independently and consult with domain experts when applying these methods to critical research or clinical applications.

## Contributing

Contributions, bug reports, and feature requests are welcome! Please open an issue or submit a pull request on [GitHub](https://github.com/simonhenin/ClusterStats.jl).

## References

For more information on cluster-based permutation testing, see:

- Maris, E., & Oostenveld, R. (2007). Nonparametric statistical testing of EEG- and MEG-data. *Journal of Neuroscience Methods*, 164(1), 177-190.