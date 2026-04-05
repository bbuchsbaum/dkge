# Extract anchor diagnostics from a DKGE fit

Extract anchor diagnostics from a DKGE fit

## Usage

``` r
dkge_anchor_diagnostics(fit)
```

## Arguments

- fit:

  Object returned by \[dkge_anchor_fit()\] or
  \[dkge_fit_from_kernels()\].

## Value

List containing per-subject coverage quantiles and per-anchor leverage
estimates.

## Examples

``` r
# \donttest{
set.seed(1)
features_list <- list(
  s1 = matrix(rnorm(20 * 5), 20, 5),
  s2 = matrix(rnorm(25 * 5), 25, 5),
  s3 = matrix(rnorm(22 * 5), 22, 5)
)
K_item_list <- lapply(features_list, function(X) tcrossprod(matrix(rnorm(nrow(X) * 3), nrow(X), 3)))
fit <- dkge_anchor_fit(features_list, K_item_list,
                       anchors = list(L = 6, method = "dkpp", seed = 1),
                       dkge_args = list(rank = 2))
#> Warning: Subject 's1': beta matrix has reduced rank (3 < 6 effects).
#> Warning: Subject 's2': beta matrix has reduced rank (3 < 6 effects).
#> Warning: Subject 's3': beta matrix has reduced rank (3 < 6 effects).
diag <- dkge_anchor_diagnostics(fit)
names(diag)
#> [1] "summary"  "coverage" "leverage"
# }
```
