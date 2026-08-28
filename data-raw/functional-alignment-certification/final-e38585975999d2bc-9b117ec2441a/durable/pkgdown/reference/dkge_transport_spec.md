# Legacy descriptive transport specification helper

Builds a validated descriptive transport configuration for
[`dkge_pipeline()`](https://bbuchsbaum.github.io/dkge/reference/dkge_pipeline.md).
Pipeline transport cannot enter inference. New functional alignment
workflows should use
[`dkge_prepare_alignment()`](https://bbuchsbaum.github.io/dkge/reference/dkge_prepare_alignment.md)
or
[`dkge_transport_contrasts_to_reference()`](https://bbuchsbaum.github.io/dkge/reference/dkge_transport_contrasts_to_reference.md).

## Usage

``` r
dkge_transport_spec(
  centroids,
  sizes = NULL,
  medoid = 1L,
  method = c("sinkhorn", "ridge", "ols", "sinkhorn_cpp"),
  mapper = NULL,
  epsilon = 0.05,
  max_iter = 5000L,
  tol = 1e-04,
  lambda_emb = 1,
  lambda_spa = 0.5,
  sigma_mm = 15,
  lambda_size = 0,
  value_type = c("intensive", "extensive"),
  warm_start = TRUE,
  ...
)
```

## Arguments

- centroids:

  List of subject-specific centroid matrices (P_s x d).

- sizes:

  Optional list of cluster sizes (one numeric vector per subject).

- medoid:

  Legacy integer reference-subject index (default 1). This field fixes a
  support; it does not perform or certify medoid selection. New
  functional-alignment workflows should use
  [`dkge_prepare_alignment()`](https://bbuchsbaum.github.io/dkge/reference/dkge_prepare_alignment.md).

- method:

  Mapper backend. Default "sinkhorn".

- mapper:

  Optional prefit mapper specification (advanced use).

- epsilon:

  Sinkhorn entropic regularisation parameter.

- max_iter:

  Maximum Sinkhorn iterations.

- tol:

  Convergence tolerance for Sinkhorn scaling.

- lambda_emb:

  Weight on embedding distance in the cost matrix.

- lambda_spa:

  Weight on spatial distance in the cost matrix.

- sigma_mm:

  Spatial scale (in millimetres) used when spatial coordinates are
  available.

- lambda_size:

  Weight on size regularisation between clusters.

- value_type:

  Semantics of transported values: `"intensive"` preserves constant
  fields; `"extensive"` preserves the sum of source totals.

- warm_start:

  Reuse converged Sinkhorn duals for identical problems.

- ...:

  Additional fields stored on the spec (e.g., precomputed loadings or
  betas).

## Value

Object with class `dkge_transport_spec`.

## Examples

``` r
transport <- dkge_transport_spec(centroids = list(matrix(runif(12), 4, 3)))
```
