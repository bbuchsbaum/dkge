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
  length P or PxP matrix). For ClusteredNeuroVec method: omega defaults
  to cluster sizes if not provided

## Value

Object of class \`dkge_subject\`

## Examples

``` r
betas <- matrix(rnorm(5 * 200), 5, 200)
design <- matrix(rnorm(150 * 5), 150, 5, dimnames = list(NULL, paste0("eff", 1:5)))
subj <- dkge_subject(betas, design, id = "sub01")
str(subj)
#> List of 7
#>  $ id         : chr "sub01"
#>  $ beta       : num [1:5, 1:200] -0.7104 -0.0651 -1.7595 0.5697 1.6123 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:5] "eff1" "eff2" "eff3" "eff4" ...
#>   .. ..$ : chr [1:200] "cluster_1" "cluster_2" "cluster_3" "cluster_4" ...
#>  $ design     : num [1:150, 1:5] 0.329 1.88 2.043 1.329 -0.22 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : chr [1:5] "eff1" "eff2" "eff3" "eff4" ...
#>  $ omega      : NULL
#>  $ effects    : chr [1:5] "eff1" "eff2" "eff3" "eff4" ...
#>  $ n_clusters : int 200
#>  $ cluster_ids: chr [1:200] "cluster_1" "cluster_2" "cluster_3" "cluster_4" ...
#>  - attr(*, "class")= chr "dkge_subject"
```
