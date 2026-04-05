# Fit DKGE from an input descriptor

Dispatches to the appropriate preprocessing pipeline (anchors, raw
betas, etc.) before invoking the DKGE core. Currently supports anchor
descriptors created via \[dkge_input_anchor()\].

## Usage

``` r
dkge_fit_from_input(input, ...)

# S3 method for class 'dkge_input_anchor'
dkge_fit_from_input(input, ...)

# Default S3 method
dkge_fit_from_input(input, ...)
```

## Arguments

- input:

  Object inheriting from \`dkge_input\`.

- ...:

  Additional arguments merged into the DKGE fitting call (these are
  interpreted as DKGE core options, e.g. \`w_method\`).

## Value

A fitted \`dkge\` object.

## Examples

``` r
# \donttest{
set.seed(1)
features_list <- list(
  s1 = matrix(rnorm(20 * 4), 20, 4),
  s2 = matrix(rnorm(25 * 4), 25, 4),
  s3 = matrix(rnorm(22 * 4), 22, 4)
)
K_item_list <- lapply(features_list, function(X) tcrossprod(matrix(rnorm(nrow(X) * 3), nrow(X), 3)))
input <- dkge_input_anchor(features_list, K_item_list,
                           anchors = list(L = 6, method = "dkpp", seed = 1),
                           dkge_args = list(rank = 2))
fit <- dkge_fit_from_input(input)
#> Warning: Subject 's1': beta matrix has reduced rank (3 < 6 effects).
#> Warning: Subject 's2': beta matrix has reduced rank (3 < 6 effects).
#> Warning: Subject 's3': beta matrix has reduced rank (3 < 6 effects).
fit$rank
#> [1] 2
# }
```
