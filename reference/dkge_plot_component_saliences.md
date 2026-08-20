# Plot raw DKGE component saliences

Plot raw DKGE component saliences

## Usage

``` r
dkge_plot_component_saliences(
  fit,
  comps = NULL,
  scale = c("raw", "unit", "zscore"),
  type = c("heatmap", "profile")
)
```

## Arguments

- fit:

  Fitted \`dkge\` object.

- comps:

  Components to include (defaults to first min(rank, 6)).

- scale:

  Optional within-component display scaling. \`"raw"\` leaves saliences
  on their original scale, \`"unit"\` rescales each component to unit
  K-norm of the underlying latent vector (see Details), and \`"zscore"\`
  z-scores each component across effects.

- type:

  Plot type. \`"heatmap"\` is the most general view; \`"profile"\` draws
  one line per component across ordered effect rows.

## Value

A ggplot object.

## Examples

``` r
toy <- dkge_sim_toy(
  factors = list(A = list(L = 2), B = list(L = 3)),
  active_terms = c("A", "B"), S = 3, P = 15, snr = 5
)
fit <- dkge(toy$B_list, toy$X_list, K = toy$K, rank = 2)
dkge_plot_component_saliences(fit, comps = 1:2)
```
