# Multiplier bootstrap in the design space (q-space)

Reweights subject contributions with i.i.d. multiplier weights,
recomputes the tiny qxq eigendecomposition, and propagates contrasts to
an identified reference (and optionally voxel) support using a validated
typed fitted alignment. The cohort-trained truncated span makes this
route approximate; callers must opt in with
`allow_approximate_alignment = TRUE`. Fit-level MFA/energy weights act
only on the reweighted pooled moment used to estimate each latent basis.
The returned group map has an equal-subject base estimand: bootstrap
multipliers resample subjects, but fit-level moment weights are not
reused as second-level aggregation weights.

## Usage

``` r
dkge_bootstrap_qspace(
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

- allow_approximate_alignment:

  Permit a typed alignment labelled approximate. Descriptive/ineligible
  alignment is always refused.

- ...:

  Deprecated compatibility arguments; ignored.

## Value

List containing per-contrast bootstrap summaries and the transport cache
employed during resampling.
