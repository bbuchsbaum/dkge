# Kernel alignment pre-screening

Ranks candidate kernels by their alignment with the pooled design-space
covariance produced by \[dkge_pooled_cov_q()\].

## Usage

``` r
dkge_kernel_prescreen(K_grid, C, normalize_k = TRUE, top_k = 3)
```

## Arguments

- K_grid:

  Named list of qxq kernels.

- C:

  Pooled covariance matrix.

- normalize_k:

  Logical; if \`TRUE\`, kernels are scaled to unit trace before
  alignment.

- top_k:

  Number of kernels to retain.

## Value

Data frame sorted by decreasing alignment; the \`top\` attribute carries
the names of the retained kernels.

## Examples

``` r
toy <- dkge_sim_toy(
  factors = list(cond = list(L = 3)),
  active_terms = "cond", S = 3, P = 10, snr = 5
)
pooled <- dkge_pooled_cov_q(toy$B_list, toy$X_list)
q <- nrow(toy$K)
K_grid <- list(base = toy$K, identity = diag(q))
dkge_kernel_prescreen(K_grid, pooled$C, top_k = 1)
#>     kernel     align
#> 2 identity 0.5095343
#> 1     base 0.3682251
```
