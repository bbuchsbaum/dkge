# DKGE scree plot with cumulative curve

DKGE scree plot with cumulative curve

## Usage

``` r
dkge_plot_scree(fit, one_se_pick = NULL)
```

## Arguments

- fit:

  Fitted \`dkge\` object.

- one_se_pick:

  Optional integer component chosen by one-SE rule.

## Value

A ggplot object.

## Examples

``` r
toy <- dkge_sim_toy(
  factors = list(A = list(L = 2), B = list(L = 3)),
  active_terms = c("A", "B"), S = 3, P = 15, snr = 5
)
fit <- dkge(toy$B_list, toy$X_list, kernel = toy$K, rank = 2)
#> Warning: Argument 'kernel' is deprecated; use 'K' instead.
dkge_plot_scree(fit)
```
