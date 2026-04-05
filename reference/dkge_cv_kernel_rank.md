# Combined kernel and rank selection via pre-screening and LOSO CV

Runs kernel alignment pre-screening followed by LOSO explained-variance
cross-validation, applying the one-standard-error rule to pick a
kernel/rank pair.

## Usage

``` r
dkge_cv_kernel_rank(
  B_list,
  X_list,
  K_grid,
  ranks,
  Omega_list = NULL,
  ridge = 0,
  w_method = "mfa_sigma1",
  w_tau = 0.3,
  top_k = 3
)
```

## Arguments

- B_list:

  List of qxP subject beta matrices.

- X_list:

  List of Txq subject design matrices.

- K_grid:

  Named list of candidate kernels.

- ranks:

  Integer vector of ranks to evaluate.

- Omega_list:

  Optional list of spatial weights.

- ridge:

  Optional ridge parameter passed to \[dkge_fit()\].

- w_method:

  Subject-level weighting scheme passed to \[dkge_fit()\].

- w_tau:

  Shrinkage parameter toward equal weights passed to \[dkge_fit()\].

- top_k:

  Number of kernels to keep after pre-screening.

## Value

List with the selected \`kernel\` and \`rank\`, alignment and CV tables,
and a per-kernel summary of scores at the selected rank.

## Examples

``` r
# \donttest{
toy <- dkge_sim_toy(
  factors = list(cond = list(L = 3)),
  active_terms = "cond", S = 4, P = 15, snr = 5
)
q <- nrow(toy$K)
K_grid <- list(base = toy$K, identity = diag(q))
sel <- dkge_cv_kernel_rank(toy$B_list, toy$X_list, K_grid, ranks = 1:2)
sel$pick
#> $kernel
#> [1] "base"
#> 
#> $rank
#> [1] 2
#> 
# }
```
