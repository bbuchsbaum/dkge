# Compute cluster-level betas from neuroim2 objects

Compute cluster-level betas from neuroim2 objects

## Usage

``` r
dkge_cluster_betas(bv, x_mat, labels)
```

## Arguments

- bv:

  A \`neuroim2::NeuroVec\` object containing 4D time-series data.

- x_mat:

  Subject design matrix (Txq).

- labels:

  A \`neuroim2::NeuroVol\` object with cluster label assignments.

## Value

Matrix of GLM betas (qxP) where P is the number of clusters.

## Examples

``` r
# \donttest{
if (requireNamespace("neuroim2", quietly = TRUE)) {
  labels <- neuroim2::read_vol(system.file("extdata", "global_mask2.nii.gz", package = "neuroim2"))
  vols <- lapply(1:5, function(i) labels * i)
  bv <- neuroim2::vec_from_vols(vols)
  x_mat <- cbind(intercept = 1, trend = seq_len(length(vols)))
  betas <- dkge_cluster_betas(bv, x_mat, labels)
  dim(betas)
}
#> Registered S3 method overwritten by 'fmridesign':
#>   method               from   
#>   print.sampling_frame fmrihrf
#> [1] 1 2
# }
```
