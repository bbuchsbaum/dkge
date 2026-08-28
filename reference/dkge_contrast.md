# Compute DKGE contrasts with cross-fitting

Main entry point for computing design contrasts after DKGE fit, with
support for multiple cross-fitting strategies to ensure unbiased
estimation.

## Usage

``` r
dkge_contrast(
  fit,
  contrasts,
  method = c("loso", "kfold", "analytic"),
  folds = NULL,
  ridge = 0,
  parallel = FALSE,
  verbose = FALSE,
  align = TRUE,
  transport = NULL,
  collinearity_tol = 0.05,
  ...
)
```

## Arguments

- fit:

  A `dkge` object from
  [`dkge_fit()`](https://bbuchsbaum.github.io/dkge/reference/dkge_fit.md)
  or [`dkge()`](https://bbuchsbaum.github.io/dkge/reference/dkge.md)

- contrasts:

  Either a q-length numeric vector for a single contrast, a named list
  of contrasts, or a qxk matrix where columns are contrasts

- method:

  Cross-fitting strategy: "loso" (leave-one-subject-out), "kfold"
  (K-fold cross-validation), or "analytic" (first-order approximation)

- folds:

  For method="kfold", either an integer K for random folds, or a list
  defining custom fold assignments (see
  [`dkge_define_folds()`](https://bbuchsbaum.github.io/dkge/reference/dkge_define_folds.md))

- ridge:

  Optional ridge parameter added when recomputing held-out basis

- parallel:

  Logical; if TRUE uses
  [`future.apply::future_lapply()`](https://future.apply.futureverse.org/reference/future_lapply.html)
  for per-subject work (requires the future.apply package)

- verbose:

  Logical; print progress messages

- align:

  Logical; if TRUE (default) align LOSO/K-fold bases to a common
  reference and compute a consensus basis for reporting.

- transport:

  Optional list describing how to transport subject-level contrasts to a
  shared reference parcellation. Supply either explicit transport
  matrices via `transforms`/`matrices`, or configuration for the
  medoid/atlas transport helpers (e.g., `method`, `centroids`,
  `medoid`). When provided, the resulting transport bundle is stored
  under `metadata$transport` for downstream reuse.

- collinearity_tol:

  Relative angular tolerance used to flag distinct input contrasts whose
  kernel-transformed queries are practically proportional. The default
  `0.05` flags absolute query correlations of at least `0.95`; use
  `NULL` to disable this warning. This is deliberately separate from the
  much smaller null-space tolerance used for estimability.

- ...:

  Additional arguments passed to method-specific functions

## Value

A list with class `dkge_contrasts` containing:

- `values`: Named list of contrast values (one P_s vector per subject
  per contrast)

- `method`: Cross-fitting method used

- `contrasts`: Input contrast specifications

- `metadata`: Method-specific metadata (fold assignments, bases, etc.)

## Details

This function provides a unified interface to three cross-fitting
strategies:

1.  **LOSO** (`method = "loso"`): Recomputes the basis excluding each
    subject, then projects that subject's data. This is the gold
    standard for unbiased estimation but requires S
    eigen-decompositions.

2.  **K-fold** (`method = "kfold"`): Splits data into K folds,
    recomputes basis excluding each fold, projects held-out data. More
    efficient than LOSO while maintaining good bias properties. Supports
    time-based, run-based, or custom fold definitions.

3.  **Analytic** (`method = "analytic"`): Uses first-order eigenvalue
    perturbation theory to approximate the LOSO solution without full
    recomputation. Fast but may be less accurate when subjects have high
    leverage.

All methods work entirely in the qxq design space and respect the
K-metric throughout. A semidefinite kernel defines a quotient effect
space: `dkge_contrast()` errors when a contrast lies wholly in
`null(K)`, warns when only part of a contrast is represented, and
reports transformed-query collisions in `metadata$kernel_query_pairs`.
Multiple contrasts can be evaluated simultaneously for efficiency.

Exact fold replay currently supports fits made with `solver = "pooled"`
and `cpca_part = "none"`. CPCA and joint-diagonalization fits fail
closed for all three methods; DKGE does not replace their fitted
estimator with an ordinary pooled eigensolve while calling the result
cross-fitted.

## See also

[`dkge_loso_contrast()`](https://bbuchsbaum.github.io/dkge/reference/dkge_loso_contrast.md),
[`dkge_define_folds()`](https://bbuchsbaum.github.io/dkge/reference/dkge_define_folds.md),
[`dkge_infer()`](https://bbuchsbaum.github.io/dkge/reference/dkge_infer.md)

## Examples

``` r
# Simulate and fit
toy <- dkge_sim_toy(
  factors = list(A = list(L = 2), B = list(L = 3)),
  active_terms = c("A", "B"), S = 4, P = 20, snr = 5
)
fit <- dkge(toy$B_list, toy$X_list, kernel = toy$K, rank = 2)
#> Warning: Argument 'kernel' is deprecated; use 'K' instead.

# Single contrast with LOSO cross-fitting
c1 <- c(1, rep(0, 4))
result <- dkge_contrast(fit, c1, method = "loso")
result
#> DKGE Contrasts
#> --------------
#> Method: loso
#> Contrasts: 1
#> Subjects: 4
#> Contrast names: contrast1 

# Fast analytic approximation
result_fast <- dkge_contrast(fit, c1, method = "analytic")
```
