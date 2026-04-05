# Multiplier bootstrap in the design space (q-space)

Reweights subject contributions with i.i.d. multiplier weights,
recomputes the tiny qxq eigendecomposition, and propagates contrasts to
the medoid (and optionally voxel) space using cached transport
operators.

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

- ...:

  Additional arguments forwarded to \[dkge_prepare_transport()\] when
  the transport cache needs to be built.

## Value

List containing per-contrast bootstrap summaries and the transport cache
employed during resampling.
