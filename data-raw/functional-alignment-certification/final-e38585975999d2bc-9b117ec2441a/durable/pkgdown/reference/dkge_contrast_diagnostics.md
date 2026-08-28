# Diagnose whether a kernel can represent planned contrasts

Performs a read-only preflight check in the same pooled-design and
kernel geometry used by
[`dkge_contrast()`](https://bbuchsbaum.github.io/dkge/reference/dkge_contrast.md).
This is useful before running cross-fitting or inference, especially
when `K` is singular or strongly structured.

## Usage

``` r
dkge_contrast_diagnostics(fit, contrasts, tol = 1e-08, collinearity_tol = 0.05)
```

## Arguments

- fit:

  A fitted `dkge` object.

- contrasts:

  A contrast vector, matrix, or list accepted by
  [`dkge_contrast()`](https://bbuchsbaum.github.io/dkge/reference/dkge_contrast.md).

- tol:

  Numerical tolerance used only for null-space estimability checks.

- collinearity_tol:

  Relative angular tolerance used to flag practically proportional
  transformed queries. The default `0.05` corresponds to
  `abs(query_correlation) >= 0.95`; use `NULL` to disable collision
  flags.

## Value

A list with `estimability` (support and null fractions per contrast),
`pairs` (pairwise transformed-query correlations and collision flags), a
compact `summary` naming the maximally correlated pair, and scalar
`kernel` rank and spectral-concentration diagnostics.

## Examples

``` r
toy <- dkge_sim_toy(
  factors = list(cond = list(L = 3)), active_terms = "cond",
  S = 3, P = 10, snr = 4
)
fit <- dkge(toy$B_list, toy$X_list, K = toy$K, rank = 2)
c_vec <- c(1, -1, rep(0, nrow(fit$U) - 2))
dkge_contrast_diagnostics(fit, c_vec)$estimability
#>    contrast    status support_fraction null_fraction query_norm
#> 1 contrast1 estimable                1  3.697785e-32  0.4714045
```
