# Design Kernels and Model Tuning

This vignette demonstrates how the design kernel controls and guides
DKGE model fits, and provides practical guidance on tuning rank and
penalty parameters before conducting large-scale analyses. Understanding
how to construct appropriate kernels and select optimal parameters is
crucial for achieving robust and interpretable results in cluster-level
fMRI analysis.

## Simulated Working Example

``` r
library(dkge)
S <- 4; q <- 4; P <- 30; T <- 80

betas <- replicate(S, matrix(rnorm(q * P, sd = 0.5), q, P), simplify = FALSE)
true_kernel <- 0.9 * (diag(q) + 0.3)
cholK <- chol(true_kernel)
betas <- lapply(betas, function(B) cholK %*% B + 0.1 * matrix(rnorm(q * P), q, P))

designs <- replicate(S, {
  X <- matrix(rnorm(T * q), T, q)
  qr.Q(qr(X))
}, simplify = FALSE)

subjects <- lapply(seq_len(S), function(s) dkge_subject(betas[[s]], designs[[s]], id = paste0("sub", s)))
bundle <- dkge_data(subjects)
```

The
[`dkge_data()`](https://bbuchsbaum.github.io/dkge/reference/dkge_data.md)
function serves an important preprocessing role by aligning effect
ordering across subjects and caching common metadata for efficient
computation. The design kernel that we specify in subsequent steps
determines how strongly different effects are smoothed or grouped
together during the embedding process.

## Building Kernels with `design_kernel()`

``` r
factors <- list(
  A = list(L = 2, type = "nominal"),
  B = list(L = 2, type = "nominal")
)
terms <- list("A", "B", c("A", "B"))
K_struct <- design_kernel(
  factors,
  terms = terms,
  rho = c("A" = 1, "B" = 1, "A:B" = 0.4),
  basis = "cell",
  normalize = "unit_trace"
)

round(K_struct$K, 2)
#>      [,1] [,2] [,3] [,4]
#> [1,] 0.25 0.10 0.10 0.00
#> [2,] 0.10 0.25 0.00 0.10
#> [3,] 0.10 0.00 0.25 0.10
#> [4,] 0.00 0.10 0.10 0.25
```

The kernel construction process involves several key concepts. Each term
in the design introduces a block-similarity component that captures
relationships between effects. The `rho` parameter weights the
contribution of each term, with a default value of 1, allowing you to
control the relative importance of different components. By combining
multiple terms, you can encode complex relationships such as
interactions or ordered trends within your experimental design. The
[`design_kernel()`](https://bbuchsbaum.github.io/dkge/reference/design_kernel.md)
function returns not only the kernel matrix itself but also associated
metadata that will be used by the cross-validation helper functions
described in the following sections.

## Rank and Kernel Tuning via Cross-Validation

The
[`dkge_cv_kernel_rank()`](https://bbuchsbaum.github.io/dkge/reference/dkge_cv_kernel_rank.md)
function provides a systematic approach for evaluating candidate ranks
and kernels using leave-one-subject-out cross-validation objectives. The
strategy involves building a grid of candidate kernels with different
parameter settings and allowing the cross-validation procedure to
identify the optimal combination based on predictive performance.

``` r
rho_vals <- seq(0.4, 1.0, by = 0.3)
K_grid <- list()
for (ra in rho_vals) {
  for (rb in rho_vals) {
    nm <- sprintf("A%.1f_B%.1f", ra, rb)
    K_grid[[nm]] <- design_kernel(
      factors,
      terms = terms,
      rho = c("A" = ra, "B" = rb, "A:B" = 0.4),
      basis = "cell",
      normalize = "unit_trace"
    )$K
  }
}

cv_rows <- lapply(names(K_grid), function(nm) {
  cv <- dkge_cv_rank_loso(bundle$betas, bundle$designs, K_grid[[nm]], ranks = 1:3)
  cbind(kernel = nm, cv$table)
})
cv_table <- do.call(rbind, cv_rows)
head(cv_table)
#>      kernel param  mean      se
#> 1 A0.4_B0.4     1 0.582 0.02751
#> 2 A0.4_B0.4     2 0.761 0.01993
#> 3 A0.4_B0.4     3 0.935 0.01012
#> 4 A0.4_B0.7     1 0.598 0.03118
#> 5 A0.4_B0.7     2 0.812 0.01427
#> 6 A0.4_B0.7     3 0.948 0.00808
```

The cross-validation procedure reports a score that represents the
average reconstruction error of left-out subjects measured in the
K-metric space. Lower scores indicate better predictive performance, as
they reflect more accurate reconstruction of held-out data. Based on
these results, you should select the rank and kernel combination that
achieves the highest cross-validated score, then refit the model using
these optimal parameters.

``` r
best_idx <- which.max(cv_table$mean)
best <- cv_table[best_idx, ]
best
#>       kernel param  mean      se
#> 27 A1.0_B1.0     3 0.968 0.00501

K_best <- K_grid[[best$kernel]]
fit <- dkge(bundle, K = K_best, rank = best$param)
round(fit$sdev, 3)
#> [1] 9.83 5.36 4.94
```

## Penalising Noisy Effects

In many neuroimaging studies, certain effects are known to exhibit high
variance characteristics, such as motion regressors or other nuisance
variables. To handle these problematic effects, you can add diagonal
ridge terms through the `dkge(..., ridge = value)` parameter or
alternatively pre-scale the kernel matrix to de-emphasize these
components. Higher ridge values have the effect of shrinking small
singular values toward zero, which helps stabilize model fits
particularly when working with limited subject counts or noisy data.

``` r
fit_ridge <- dkge(bundle, K = K_best, rank = best$param, ridge = 0.2)
round(fit_ridge$sdev, 3)
#> [1] 9.84 5.38 4.96
```

## Summary

This vignette has demonstrated several key principles for effective DKGE
modeling. First, use the
[`design_kernel()`](https://bbuchsbaum.github.io/dkge/reference/design_kernel.md)
function to encode your scientific priors about effect relationships. If
you are uncertain about the appropriate kernel structure, start with a
simple diagonal matrix `diag(q)` and then gradually add structured terms
as you develop understanding of your data. Second, always run
[`dkge_cv_kernel_rank()`](https://bbuchsbaum.github.io/dkge/reference/dkge_cv_kernel_rank.md)
with leave-one-subject-out cross-validation (the default) to
systematically select both rank and kernel parameters before proceeding
to computationally expensive inference procedures. Finally, consider
adding ridge regularization when singular values drop sharply or when
working with small subject samples, and regularly check the `fit$sdev`
values to verify numerical stability of your fitted models.
