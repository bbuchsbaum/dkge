# Aggregate voxel time series into cluster means

Aggregate voxel time series into cluster means

## Usage

``` r
dkge_cluster_ts(bv, labels, ids = NULL, chunker = NULL)
```

## Arguments

- bv:

  A \`neuroim2::NeuroVec\` object (4D time-series data, TxXxYxZ).

- labels:

  A \`neuroim2::NeuroVol\` with integer cluster identifiers (0 indicates
  background).

- ids:

  Optional subset of cluster IDs to retain.

- chunker:

  Optional function returning list(mat = TxV_block, vox_idx).

## Value

Matrix of dimension TxP where P is the number of clusters retained.

## Examples

``` r
# \donttest{
if (requireNamespace("neuroim2", quietly = TRUE)) {
  labels <- neuroim2::read_vol(system.file("extdata", "global_mask2.nii.gz", package = "neuroim2"))
  vols <- lapply(1:4, function(i) labels * i)
  bv <- neuroim2::vec_from_vols(vols)
  ts <- dkge_cluster_ts(bv, labels)
  dim(ts)
}
#> [1] 4 1
# }
```
