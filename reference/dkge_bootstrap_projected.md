# Subject-level projection bootstrap in medoid space

Resamples transported subject vectors (already aligned in the medoid
parcellation) to quantify between-subject variability without
recomputing the group basis.

## Usage

``` r
dkge_bootstrap_projected(
  values_medoid,
  B = 1000L,
  aggregate = c("mean", "median"),
  weights = NULL,
  seed = NULL,
  voxel_operator = NULL,
  return_samples = TRUE
)
```

## Arguments

- values_medoid:

  List of length \`S\` where each element is a numeric vector defined on
  the medoid parcellation (e.g. LOSO contrast values).

- B:

  Number of bootstrap replicates.

- aggregate:

  Aggregation function applied to the resampled subjects (\`"mean"\` or
  \`"median"\`).

- weights:

  Optional subject weights applied when computing the resampled mean.
  Only used when \`aggregate = "mean"\`.

- seed:

  Optional random seed for reproducibility.

- voxel_operator:

  Optional matrix that maps medoid vectors to voxel space (columns =
  voxels). When supplied, summaries in voxel space are also returned.

- return_samples:

  Logical; when \`TRUE\` the matrix of bootstrap samples is returned in
  the output bundle.

## Value

A list containing bootstrap summaries (\`mean\`, \`sd\`, \`z\`,
confidence intervals), and optionally the raw bootstrap draws (medoid
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
