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
  allow_approximate_alignment = FALSE,
  ...
)
```

## Arguments

- fit:

  A fitted `dkge` object.

- contrasts:

  Contrast specification accepted by
  [`dkge_contrast()`](https://bbuchsbaum.github.io/dkge/reference/dkge_contrast.md).

- B:

  Number of bootstrap replicates.

- scheme:

  Multiplier distribution (`"poisson"`, `"exp"`, or `"bayes"`). Poisson
  draws are conditioned on at least one positive subject multiplier, so
  every bootstrap replicate contains a non-empty resampled cohort.

- ridge:

  Optional ridge added to the reweighted compressed covariance.

- align:

  Logical; when `TRUE` the resampled bases are aligned to the baseline
  basis via K-Procrustes before contrasts are evaluated.

- allow_reflection:

  Passed to
  [`dkge_procrustes_K()`](https://bbuchsbaum.github.io/dkge/reference/dkge_procrustes_K.md)
  when aligning bases.

- seed:

  Optional random seed for reproducibility.

- transport_cache:

  Required typed fitted alignment object produced by
  [`dkge_prepare_alignment()`](https://bbuchsbaum.github.io/dkge/reference/dkge_prepare_alignment.md).
  Untyped/legacy caches and automatic full-fit correspondence are
  refused.

- mapper:

  Deprecated; automatic mapper fitting is no longer permitted.

- centroids:

  Deprecated; correspondence must be supplied through `transport_cache`.

- sizes:

  Deprecated compatibility argument; correspondence and masses must
  already be fixed in `transport_cache`.

- medoid:

  Deprecated compatibility argument.

- voxel_operator:

  Optional matrix mapping medoid values to voxels.

- perturb_tol:

  Maximum absolute coefficient tolerated in the eigenvector
  perturbation; larger changes trigger a fallback to the full
  eigensolve.

- gap_tol:

  Minimum eigen-gap tolerated (in absolute value) before triggering a
  fallback to the full eigensolve.

- allow_approximate_alignment:

  Permit a typed alignment labelled approximate. Descriptive/ineligible
  alignment is always refused.

- ...:

  Deprecated compatibility arguments; ignored.

## Value

Same structure as
[`dkge_bootstrap_qspace()`](https://bbuchsbaum.github.io/dkge/reference/dkge_bootstrap_qspace.md)
with additional metadata on the number of fallbacks used.

## Details

As in
[`dkge_bootstrap_qspace()`](https://bbuchsbaum.github.io/dkge/reference/dkge_bootstrap_qspace.md),
stored fit weights affect only the reweighted pooled moment. Subject
maps are aggregated with an equal-subject base estimand using the
bootstrap multipliers alone.
