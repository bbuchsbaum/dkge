# Subject-level projection bootstrap on a reference support

Resamples subject vectors already carried by a typed aligned-map object
to quantify between-subject variability without recomputing the group
basis or correspondence.

## Usage

``` r
dkge_bootstrap_projected(
  values_medoid,
  contrast = 1L,
  B = 1000L,
  aggregate = c("mean", "median"),
  weights = NULL,
  seed = NULL,
  voxel_operator = NULL,
  return_samples = TRUE,
  allow_approximate_alignment = FALSE
)
```

## Arguments

- values_medoid:

  An operator-bound
  [`dkge_aligned_maps()`](https://bbuchsbaum.github.io/dkge/reference/dkge_aligned_maps.md)
  object returned by a fitted alignment application. Raw lists and
  objects from the public descriptive constructor are refused because
  they carry no certified operator application.

- contrast:

  Contrast name or index to bootstrap from `values_medoid`.

- B:

  Number of bootstrap replicates.

- aggregate:

  Aggregation function applied to the resampled subjects (`"mean"` or
  `"median"`).

- weights:

  Optional assertion of the immutable subject weights stored in
  `values_medoid`. If supplied, it must match that receipt exactly. Only
  used when `aggregate = "mean"`.

- seed:

  Optional random seed for reproducibility.

- voxel_operator:

  Optional matrix that maps medoid vectors to voxel space (columns =
  voxels). When supplied, summaries in voxel space are also returned.

- return_samples:

  Logical; when `TRUE` the matrix of bootstrap samples is returned in
  the output bundle.

- allow_approximate_alignment:

  Permit an explicitly labelled approximate aligned-map object.
  Descriptive/ineligible states are always refused.

## Value

A list containing bootstrap summaries (`mean`, `sd`, `z`, confidence
intervals), and optionally the raw bootstrap draws (reference-support
and voxel space).

## Examples

``` r
# \donttest{
toy <- dkge_sim_toy(
  factors = list(A = list(L = 2), B = list(L = 3)),
  active_terms = c("A", "B"), S = 5, P = 20, snr = 5
)
fit <- dkge_fit(toy$B_list, toy$X_list, toy$K, rank = 2)
# Bootstrap requires transport setup - example shows API
# }
```
