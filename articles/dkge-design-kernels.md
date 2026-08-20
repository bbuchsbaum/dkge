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
#>       A1:B1 A1:B2 A2:B1 A2:B2
#> A1:B1  0.25  0.10  0.10  0.00
#> A1:B2  0.10  0.25  0.00  0.10
#> A2:B1  0.10  0.00  0.25  0.10
#> A2:B2  0.00  0.10  0.10  0.25
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

## Which rank and kernel predict held-out structure?

[`dkge_cv_kernel_rank()`](https://bbuchsbaum.github.io/dkge/reference/dkge_cv_kernel_rank.md)
first screens named kernels by alignment, then evaluates the retained
kernel-rank pairs with leave-one-subject-out explained energy. It
returns the selected pair and the tables needed to audit that choice.

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

selection <- dkge_cv_kernel_rank(
  bundle$betas, bundle$designs,
  K_grid = K_grid,
  ranks = 1:3,
  top_k = length(K_grid)
)
selection$picks_per_kernel
#>      kernel rank score      se
#> 1 A0.4_B0.4    3 0.935 0.01012
#> 2 A0.7_B0.7    3 0.957 0.00669
#> 3 A0.4_B0.7    3 0.948 0.00808
#> 4 A0.7_B0.4    3 0.948 0.00806
#> 5 A0.7_B1.0    3 0.964 0.00573
#> 6 A1.0_B1.0    3 0.968 0.00501
#> 7 A1.0_B0.7    3 0.964 0.00573
#> 8 A0.4_B1.0    3 0.957 0.00675
#> 9 A1.0_B0.4    3 0.957 0.00672
```

The score is the fraction of a held-out subject’s K-metric energy
captured by the training-fold basis. It is descriptive predictive
evidence, not a significance test. Higher values mean that the fitted
subspace captures more of the held-out energy. The selector applies its
documented one-standard-error rule, preferring a lower rank among
eligible pairs rather than blindly taking the largest observed mean.

``` r

best <- selection$pick
best
#> $kernel
#> [1] "A1.0_B1.0"
#> 
#> $rank
#> [1] 3

K_best <- K_grid[[best$kernel]]
fit <- dkge(bundle, K = K_best, rank = best$rank)
round(fit$sdev, 3)
#> [1] 9.83 5.36 4.94
```

## What does ridge regularization do here?

The `ridge` argument adds an isotropic diagonal term to the fitted
q-space moment before eigendecomposition. For a fixed moment this shifts
all eigenvalues equally and leaves its eigenvectors unchanged; it does
not identify or selectively down-weight noisy effects. Encode
effect-specific structure in the kernel or input weights instead. Ridge
is mainly a numerical guard for near-singular moment calculations.

``` r

fit_ridge <- dkge(bundle, K = K_best, rank = best$rank, ridge = 0.2)
round(fit_ridge$sdev, 3)
#> [1] 9.84 5.38 4.96
```

## Summary

Use
[`design_kernel()`](https://bbuchsbaum.github.io/dkge/reference/design_kernel.md)
to encode effect relationships you can defend. If the structure is
uncertain, include `diag(q)` as a baseline rather than assuming a
structured kernel must win. Compare kernels and ranks with
[`dkge_cv_kernel_rank()`](https://bbuchsbaum.github.io/dkge/reference/dkge_cv_kernel_rank.md),
interpreting its scores as held-out explained energy. Refit the selected
configuration before running inference, and treat ridge as isotropic
numerical regularization rather than an effect-selection device.
