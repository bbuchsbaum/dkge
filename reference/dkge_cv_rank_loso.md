# LOSO cross-validation for rank selection

Evaluates candidate ranks by recomputing LOSO bases and measuring
explained energy on the held-out subject in one candidate-independent
validation geometry. By default this is ordinary effect space after
pooled-design scaling, not the candidate kernel's own metric.

## Usage

``` r
dkge_cv_rank_loso(
  B_list,
  X_list,
  K,
  ranks,
  Omega_list = NULL,
  ridge = 0,
  w_method = "mfa_sigma1",
  w_tau = 0.3,
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

- K:

  qxq design kernel.

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

List containing the one-SE selection (`pick`), the best rank, aggregated
and per-fold score tables, kernel/validation diagnostics, excluded
ranks, and a saturation flag.

## Examples

``` r
# \donttest{
toy <- dkge_sim_toy(
  factors = list(cond = list(L = 3)),
  active_terms = "cond", S = 4, P = 15, snr = 5
)
cv <- dkge_cv_rank_loso(toy$B_list, toy$X_list, toy$K, ranks = 1:2)
cv$pick
#> [1] 2
# }
```
