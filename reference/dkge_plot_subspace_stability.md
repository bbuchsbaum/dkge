# Subspace stability via principal angles

Subspace stability via principal angles

## Usage

``` r
dkge_plot_subspace_stability(bases, K, consensus = NULL, labels = NULL)
```

## Arguments

- bases:

  List of basis matrices.

- K:

  Design kernel.

- consensus:

  Optional consensus basis.

- labels:

  Optional labels.

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
bases <- list(fit$U, fit$U)
dkge_plot_subspace_stability(bases, K = fit$K)
```
