# LOSO cross-validation for rank selection

Evaluates candidate ranks by recomputing LOSO bases and measuring
explained variance on the held-out subject in the \\K^{1/2}\\ metric.

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
  w_tau = 0.3
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

  Optional ridge parameter passed to \[dkge_fit()\].

- w_method:

  Subject-level weighting scheme passed to \[dkge_fit()\].

- w_tau:

  Shrinkage parameter toward equal weights passed to \[dkge_fit()\].

## Value

List containing the one-SE selection (\`pick\`), the best rank, and the
aggregated score table.

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
