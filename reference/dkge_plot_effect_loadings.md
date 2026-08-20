# Effect-space loadings heatmap

Displays the design-kernel-weighted component basis \\K U\\.

## Usage

``` r
dkge_plot_effect_loadings(fit, comps = NULL, zscore = FALSE)
```

## Arguments

- fit:

  Fitted \`dkge\` object.

- comps:

  Components to include (defaults to first min(rank,6)).

- zscore:

  Logical; z-score loadings within each effect.

## Value

A ggplot object.

## Examples

``` r
toy <- dkge_sim_toy(
  factors = list(A = list(L = 2), B = list(L = 3)),
  active_terms = c("A", "B"), S = 3, P = 15, snr = 5
)
fit <- dkge(toy$B_list, toy$X_list, K = toy$K, rank = 2)
dkge_plot_effect_loadings(fit, comps = 1:2)
```
