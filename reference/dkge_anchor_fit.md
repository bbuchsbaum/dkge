# Fit DKGE using feature-anchored subject kernels

Projects item-level kernels onto a shared anchor basis (via
\[dkge_build_anchor_kernels()\]) and reuses \[dkge_fit_from_kernels()\]
to enter the DKGE pipeline. Anchor provenance, coverage diagnostics, and
subject item counts are stored in the resulting \`dkge\` object.

## Usage

``` r
dkge_anchor_fit(
  features_list,
  K_item_list,
  folds = NULL,
  anchors = list(),
  design_kernel = NULL,
  dkge_args = list()
)
```

## Arguments

- features_list:

  List of subject feature matrices (\`n_s x d\` each).

- K_item_list:

  List of subject item kernels (\`n_s x n_s\` each).

- folds:

  Optional fold specification passed to \[dkge_build_anchor_kernels()\].

- anchors:

  Named list overriding anchor-building defaults (\`L\`, \`method\`,
  \`rho\`, \`fill\`, \`seed\`, \`sigma\`, \`center\`, \`whiten\`,
  \`eps\`, \`unit_trace\`, \`item_weights\`).

- design_kernel:

  Optional design kernel forwarded to \[dkge_fit_from_kernels()\].
  Defaults to the identity.

- dkge_args:

  Named list of additional arguments forwarded to
  \[dkge_fit_from_kernels()\].

## Value

A \`dkge\` fit with anchor provenance under
\`fit\$provenance\$anchors\`.

## Examples

``` r
# \donttest{
set.seed(1)
features_list <- list(
  s1 = matrix(rnorm(30 * 5), 30, 5),
  s2 = matrix(rnorm(40 * 5), 40, 5),
  s3 = matrix(rnorm(35 * 5), 35, 5)
)
K_item_list <- lapply(features_list, function(X) {
  Z <- matrix(rnorm(nrow(X) * 4), nrow(X), 4)
  tcrossprod(Z)
})
fit <- dkge_anchor_fit(
  features_list, K_item_list,
  anchors = list(L = 8, method = "dkpp", seed = 1),
  dkge_args = list(rank = 2)
)
#> Warning: Subject 's1': beta matrix has reduced rank (4 < 8 effects).
#> Warning: Subject 's2': beta matrix has reduced rank (4 < 8 effects).
#> Warning: Subject 's3': beta matrix has reduced rank (4 < 8 effects).
dkge_anchor_diagnostics(fit)$summary
#> $method
#> [1] "dkpp"
#> 
#> $sigma
#> [1] 2.955892
#> 
#> $L
#> [1] 8
#> 
#> $mean_item_count
#> [1] 35
#> 
# }
```
