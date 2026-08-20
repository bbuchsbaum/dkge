# Transport specification helper

Builds a validated transport configuration that can be passed to
\[dkge_pipeline()\] or transport utilities. The helper enforces basic
argument checks and provides sensible defaults for Sinkhorn-based
mapping.

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
  provenance = NULL,
  ...
)
```

## Arguments

- centroids:

  List of subject-specific centroid matrices (P_s x d).

- sizes:

  Optional list of cluster sizes (one numeric vector per subject).

- medoid:

  Integer index of the medoid subject (default 1).

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

  Semantics of transported values: \`"intensive"\` preserves constant
  fields; \`"extensive"\` preserves the sum of source totals.

- warm_start:

  Reuse converged Sinkhorn duals for identical problems.

- provenance:

  Optional \[dkge_transport_provenance()\] declaration. When omitted,
  loading-derived transport is marked descriptive and \[dkge_infer()\]
  will reject it rather than assume inferential validity.

- ...:

  Additional fields stored on the spec (e.g., precomputed loadings or
  betas).

## Value

Object with class \`dkge_transport_spec\`.

## Examples

``` r
transport <- dkge_transport_spec(centroids = list(matrix(runif(12), 4, 3)))
```
