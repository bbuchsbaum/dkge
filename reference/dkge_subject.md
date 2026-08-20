# Construct a DKGE subject record

Construct a DKGE subject record

## Usage

``` r
dkge_subject(x, ...)
```

## Arguments

- x:

  Source object containing subject-level data: - matrix: qxP beta
  coefficients (effects x clusters/voxels) - NeuroVec: 4D time-series
  data (TxXxYxZ), betas computed via GLM - ClusteredNeuroVec: Cluster
  time-series (TxK), betas computed via GLM

- ...:

  Additional arguments passed to methods. For the matrix method:
  \`design\` (Subject design matrix T_s x q), \`id\` (Optional subject
  identifier), \`omega\` (Optional cluster weights - numeric vector
  length P or PxP matrix), \`observed_rows\` (Optional observed
  effect-row indices in a global effect space; defaults to all local
  rows), \`effect_n\` (per-effect trial counts), \`effect_precision\`
  (direct per-effect precision), \`effect_noise_cov\` (q-by-q covariance
  multiplier such as \`(X'X)^-1\`), \`residual_variance\` (per-column
  residual variances), \`noise_trace\` (optional precomputed spatial
  noise trace), and \`split_betas\` (two q-by-P half estimates). For
  ClusteredNeuroVec method: omega defaults to cluster sizes if not
  provided

## Value

Object of class \`dkge_subject\`

## Examples

``` r
betas <- matrix(rnorm(5 * 200), 5, 200)
design <- matrix(rnorm(150 * 5), 150, 5, dimnames = list(NULL, paste0("eff", 1:5)))
subj <- dkge_subject(betas, design, id = "sub01")
str(subj)
#> List of 15
#>  $ id               : chr "sub01"
#>  $ beta             : num [1:5, 1:200] -0.7104 -0.0651 -1.7595 0.5697 1.6123 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:5] "eff1" "eff2" "eff3" "eff4" ...
#>   .. ..$ : chr [1:200] "cluster_1" "cluster_2" "cluster_3" "cluster_4" ...
#>  $ design           : num [1:150, 1:5] 0.329 1.88 2.043 1.329 -0.22 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : chr [1:5] "eff1" "eff2" "eff3" "eff4" ...
#>  $ omega            : NULL
#>  $ effects          : chr [1:5] "eff1" "eff2" "eff3" "eff4" ...
#>  $ observed_rows    : int [1:5] 1 2 3 4 5
#>  $ effect_n         : NULL
#>  $ effect_precision : NULL
#>  $ effect_noise_cov : NULL
#>  $ residual_variance: NULL
#>  $ residual_df      : NULL
#>  $ noise_trace      : NULL
#>  $ split_betas      : NULL
#>  $ n_clusters       : int 200
#>  $ cluster_ids      : chr [1:200] "cluster_1" "cluster_2" "cluster_3" "cluster_4" ...
#>  - attr(*, "class")= chr "dkge_subject"
```
