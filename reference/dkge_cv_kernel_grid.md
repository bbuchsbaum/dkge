# LOSO kernel grid search

Evaluates a named list of candidate design kernels using LOSO explained
variance at a fixed rank.

## Usage

``` r
dkge_cv_kernel_grid(
  B_list,
  X_list,
  K_grid,
  rank,
  Omega_list = NULL,
  ridge = 0,
  w_method = "mfa_sigma1",
  w_tau = 0.3
)
```

## Arguments

- B_list:

  List of qxP subject beta matrices.

- X_list:

  List of Txq subject design matrices.

- K_grid:

  Named list of candidate kernels.

- rank:

  Rank used for evaluation.

- Omega_list:

  Optional list of spatial weights.

- ridge:

  Optional ridge parameter passed to \[dkge_fit()\].

- w_method:

  Subject-level weighting scheme passed to \[dkge_fit()\].

- w_tau:

  Shrinkage parameter toward equal weights passed to \[dkge_fit()\].

## Value

List with the one-SE pick, best kernel, summary table, and raw scores.

## Examples

``` r
# \donttest{
toy <- dkge_sim_toy(
  factors = list(cond = list(L = 3)),
  active_terms = "cond", S = 4, P = 15, snr = 5
)
q <- nrow(toy$K)
K_grid <- list(base = toy$K, identity = diag(q))
cv <- dkge_cv_kernel_grid(toy$B_list, toy$X_list, K_grid, rank = 1)
cv$pick
#> [1] "base"
# }
```
