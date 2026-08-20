# Compute per-component variance explained

Returns the standard deviation, variance, and cumulative variance
explained by the DKGE components extracted in \[dkge_fit()\].

## Usage

``` r
dkge_variance_explained(fit, relative_to = c("kept", "total"))
```

## Arguments

- fit:

  A \`dkge\` object.

- relative_to:

  Compute variance proportions relative to "kept" (default, only the
  retained components) or "total" (all possible components).

## Value

Data frame with columns \`component\`, \`sdev\`, \`variance\`,
\`prop_var\`, and \`cum_prop_var\`.

## Examples

``` r
toy <- dkge_sim_toy(
  factors = list(A = list(L = 2), B = list(L = 3)),
  active_terms = c("A", "B"), S = 3, P = 15, snr = 5
)
fit <- dkge(toy$B_list, toy$X_list, kernel = toy$K, rank = 2)
#> Warning: Argument 'kernel' is deprecated; use 'K' instead.
dkge_variance_explained(fit)
#>   component      sdev  variance  prop_var cum_prop_var
#> 1         1 10.506306 110.38246 0.5383834    0.5383834
#> 2         2  9.728478  94.64328 0.4616166    1.0000000
```
