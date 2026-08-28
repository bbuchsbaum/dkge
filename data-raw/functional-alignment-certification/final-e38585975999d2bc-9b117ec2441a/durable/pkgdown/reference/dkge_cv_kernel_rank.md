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
  top_k = 3,
  validation_K = NULL,
  kernel_rank_policy = c("full", "allow_singular"),
  kernel_concentration_threshold = 0.75,
  spatial = NULL
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

  Optional ridge parameter passed to
  [`dkge_fit()`](https://bbuchsbaum.github.io/dkge/reference/dkge_fit.md).

- w_method:

  Subject-level weighting scheme passed to
  [`dkge_fit()`](https://bbuchsbaum.github.io/dkge/reference/dkge_fit.md).

- w_tau:

  Shrinkage parameter toward equal weights passed to
  [`dkge_fit()`](https://bbuchsbaum.github.io/dkge/reference/dkge_fit.md).

- top_k:

  Number of kernels to keep after pre-screening.

- validation_K:

  Optional fixed q by q PSD matrix used to score every candidate. `NULL`
  (default) uses the identity. It affects scoring only; candidate `K`
  still determines the learned training subspace.

- kernel_rank_policy:

  Kernel-support policy for selection. The default, `"full"`, requires
  `rank(K) = q`. Use `"allow_singular"` only when every candidate is
  intentionally allowed to define a quotient effect space.

- kernel_concentration_threshold:

  Optional warning threshold for spectral concentration. The default
  `0.75` warns when the selected kernel's participation-ratio effective
  rank divided by the selected latent rank is below `0.75`. This does
  not change the CV score or selection; use `NULL` to disable the
  warning.

- spatial:

  Optional fixed
  [`dkge_spatial_regularizer()`](https://bbuchsbaum.github.io/dkge/reference/dkge_spatial_regularizer.md)
  applied inside every candidate fit and training fold.

## Value

List with the selected `kernel` and `rank`, alignment and CV tables
carrying kernel-support and spectral-concentration diagnostics,
per-kernel selections, exclusions, fixed validation diagnostics,
saturation flag, and selected-kernel concentration summary.

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
#> Warning: All comparable held-out scores are >= 0.999 in the fixed validation geometry; this criterion provides little discrimination among candidates.
sel$pick
#> $kernel
#> [1] "base"
#> 
#> $rank
#> [1] 2
#> 
# }
```
