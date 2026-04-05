# Analytic first-order bootstrap in the design space

Uses the stored full eigendecomposition to apply first-order
perturbations for each bootstrap draw. When the perturbation exceeds the
validity region, the method falls back to the exact multiplier bootstrap
for that replicate.

## Usage

``` r
dkge_bootstrap_analytic(
  fit,
  contrasts,
  B = 1000L,
  scheme = c("poisson", "exp", "bayes"),
  ridge = 0,
  align = TRUE,
  allow_reflection = FALSE,
  seed = NULL,
  transport_cache = NULL,
  mapper = "sinkhorn",
  centroids = NULL,
  sizes = NULL,
  medoid = 1L,
  voxel_operator = NULL,
  perturb_tol = 0.2,
  gap_tol = 1e-06,
  ...
)
```

## Arguments

- fit:

  A fitted \`dkge\` object.

- contrasts:

  Contrast specification accepted by \[dkge_contrast()\].

- B:

  Number of bootstrap replicates.

- scheme:

  Multiplier distribution (\`"poisson"\`, \`"exp"\`, or \`"bayes"\`).

- ridge:

  Optional ridge added to the reweighted compressed covariance.

- align:

  Logical; when \`TRUE\` the resampled bases are aligned to the baseline
  basis via K-Procrustes before contrasts are evaluated.

- allow_reflection:

  Passed to \[dkge_procrustes_K()\] when aligning bases.

- seed:

  Optional random seed for reproducibility.

- transport_cache:

  Optional cache from \[dkge_prepare_transport()\].

- mapper:

  Mapper specification used when a cache is not supplied.

- centroids:

  Optional centroids overriding those stored on the fit.

- sizes:

  Optional list of cluster weights passed to
  \[dkge_prepare_transport()\].

- medoid:

  Medoid index used during transport.

- voxel_operator:

  Optional matrix mapping medoid values to voxels.

- perturb_tol:

  Maximum absolute coefficient tolerated in the eigenvector
  perturbation; larger changes trigger a fallback to the full
  eigensolve.

- gap_tol:

  Minimum eigen-gap tolerated (in absolute value) before triggering a
  fallback to the full eigensolve.

- ...:

  Additional arguments forwarded to \[dkge_prepare_transport()\] when
  the transport cache needs to be built.

## Value

Same structure as \[dkge_bootstrap_qspace()\] with additional metadata
on the number of fallbacks used.
