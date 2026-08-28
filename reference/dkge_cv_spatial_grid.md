# LOSO selection of spatial regularization strength

Selects the Laplacian penalty `lambda` while keeping the spatial graph,
design kernel, latent rank, and validation geometry fixed. Each
candidate is fitted inside each LOSO training fold, but scored against
the **unsmoothed** held-out beta block in a candidate-independent
effect-space geometry. Fixed effect scaling, adaptive location weights,
and `Omega_list` remain part of that held-out geometry; only the
candidate Laplacian solve is omitted. Thus a larger `lambda` cannot
improve its score merely by smoothing or shrinking the same field used
in the denominator.

## Usage

``` r
dkge_cv_spatial_grid(
  B_list,
  X_list,
  K,
  spatial,
  lambdas,
  rank,
  Omega_list = NULL,
  ridge = 0,
  w_method = "mfa_sigma1",
  w_tau = 0.3,
  effect_scaling = c("pooled_design", "none"),
  validation_K = NULL,
  kernel_rank_policy = c("full", "allow_singular"),
  kernel_concentration_threshold = 0.75
)
```

## Arguments

- B_list:

  List of qxP subject beta matrices.

- X_list:

  List of Txq subject design matrices.

- K:

  qxq design kernel.

- spatial:

  A
  [`dkge_spatial_regularizer()`](https://bbuchsbaum.github.io/dkge/reference/dkge_spatial_regularizer.md)
  supplying the fixed graph. Its stored `lambda` is ignored while
  evaluating `lambdas`.

- lambdas:

  Non-negative finite candidate penalties.

- rank:

  One positive latent rank used for every candidate.

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

- effect_scaling:

  Effect-space scaling passed to every candidate
  [`dkge_fit()`](https://bbuchsbaum.github.io/dkge/reference/dkge_fit.md).
  Use the same setting planned for the final refit.

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

## Value

A list with selected `pick`, unconstrained `best`, summary `table`,
per-fold `raw` scores, the selected `spatial` specification, validation
and kernel diagnostics, the one-SE threshold, recorded `fit_settings`,
and a saturation flag.

## Details

The returned `pick` is the largest (smoothest) candidate whose mean
score is within one standard error of the best candidate. Include `0` in
`lambdas` to compare against the unsmoothed model. A saturated score is
reported and warned about in the same way as other DKGE CV helpers.

Spatial CV fails closed when any positive candidate is requested but
every supplied graph is edgeless: all candidate resolvents would be the
identity, so `lambda` is not identifiable from the score. If only some
subject-specific graphs are edgeless, CV warns once with their
identifiers and continues with status `partial`.

## Examples

``` r
# \donttest{
toy <- dkge_sim_toy(
  factors = list(cond = list(L = 3)), active_terms = "cond",
  S = 4, P = 15, snr = 4
)
coords <- cbind(x = seq_len(ncol(toy$B_list[[1]])), y = 0, z = 0)
spatial <- dkge_spatial_regularizer(
  coords, lambda = 1, dthresh = 1.01, nnk = 3,
  weight_mode = "binary"
)
cv <- dkge_cv_spatial_grid(
  toy$B_list, toy$X_list, toy$K, spatial,
  lambdas = c(0, 0.25, 1), rank = 1
)
cv$pick
#> [1] 1
# }
```
