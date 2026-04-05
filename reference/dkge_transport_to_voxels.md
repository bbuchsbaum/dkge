# Transport DKGE quantities directly to voxel space

Transport DKGE quantities directly to voxel space

## Usage

``` r
dkge_transport_to_voxels(
  fit,
  values,
  voxels,
  coords = NULL,
  mapper = "ridge",
  sizes = NULL,
  ...
)
```

## Arguments

- fit:

  A \`dkge\` object containing subject loadings and centroids.

- values:

  List of subject value vectors (one per subject, length P_s).

- voxels:

  List of voxel feature matrices.

- coords:

  Optional list of voxel coordinate matrices.

- mapper:

  Mapper specification or shorthand (defaults to \`"ridge"\`).

- sizes:

  Optional list of cluster masses (one vector per subject).

- ...:

  Additional mapper parameters.

## Value

List with \`subj_values\` (S x V matrix) and \`value\` (mean across
subjects).

## Examples

``` r
# \donttest{
toy <- dkge_sim_toy(
  factors = list(A = list(L = 2), B = list(L = 3)),
  active_terms = c("A", "B"), S = 3, P = 15, snr = 5
)
fit <- dkge(toy$B_list, toy$X_list, kernel = toy$K, rank = 2)
#> Warning: Argument 'kernel' is deprecated; use 'K' instead.
fit$centroids <- lapply(toy$B_list, function(B) matrix(rnorm(ncol(B) * 3), ncol(B), 3))
values <- lapply(toy$B_list, function(B) rnorm(ncol(B)))
voxels <- lapply(toy$B_list, function(B) matrix(rnorm(10 * fit$rank), 10, fit$rank))
out <- dkge_transport_to_voxels(fit, values = values, voxels = voxels, mapper = "ridge")
dim(out$subj_values)
#> [1]  3 10
# }
```
