# Project subject clusters into the DKGE latent space

For each subject, computes \\Z_s = B_s^\top K U\\ so that every row
represents a subject cluster embedded in the \\r\\-dimensional DKGE
latent space. These projections are commonly used for training
classifiers or computing Haufe-style encoders.

## Usage

``` r
dkge_project_clusters_to_latent(fit)
```

## Arguments

- fit:

  Fitted object of class \`dkge\` containing \`Btil\`, \`K\`, and \`U\`.

## Value

A list of length \`S\` (number of subjects). Element \`s\` is a \`P_s x
r\` matrix holding the latent representation of subject \`s\`'s
clusters.

## Examples

``` r
# \donttest{
toy <- dkge_sim_toy(
  factors = list(A = list(L = 2), B = list(L = 3)),
  active_terms = c("A", "B"), S = 3, P = 20, snr = 5
)
fit <- dkge(toy$B_list, toy$X_list, kernel = toy$K, rank = 2)
#> Warning: Argument 'kernel' is deprecated; use 'K' instead.
Z_list <- dkge_project_clusters_to_latent(fit)
dim(Z_list[[1]])  # P_1 x r matrix
#> [1] 20  2
# }
```
