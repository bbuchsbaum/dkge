library(testthat)

# test-kernel-invariants.R
# Property-based tests verifying mathematical invariants of kernel construction
# and root computation across all valid kernel configurations.

# ---------------------------------------------------------------------------
# 1. Symmetry property test across all factor types
# ---------------------------------------------------------------------------

test_that("design_kernel produces symmetric matrices for all factor types", {
  factor_specs <- list(
    nominal = list(type = "nominal", L = 4),
    ordinal = list(type = "ordinal", L = 4, l = 1.0),
    circular = list(type = "circular", L = 4, l = 1.0),
    continuous = list(type = "continuous", values = c(1, 2, 3, 4))
  )

  for (name in names(factor_specs)) {
    factors <- list(A = factor_specs[[name]])
    K <- design_kernel(factors, basis = "cell", normalize = "none")$K

    expect_true(isSymmetric(K, tol = 1e-10),
                info = paste("Factor type:", name))
  }
})

# ---------------------------------------------------------------------------
# 2. PSD property test across all factor types
# ---------------------------------------------------------------------------

test_that("design_kernel produces PSD matrices for all factor types", {
  # Note: circular kernels are only guaranteed PSD for short length-scales (l <= 0.5)
  # relative to the number of levels. For longer length-scales, use ridge/jitter.
  factor_specs <- list(
    nominal = list(type = "nominal", L = 4),
    ordinal = list(type = "ordinal", L = 4, l = 1.0),
    circular = list(type = "circular", L = 4, l = 0.5),  # l=0.5 ensures PSD
    continuous = list(type = "continuous", values = c(1, 2, 3, 4))
  )

  for (name in names(factor_specs)) {
    factors <- list(A = factor_specs[[name]])
    K <- design_kernel(factors, basis = "cell", normalize = "none")$K

    # Verify PSD (all eigenvalues >= -epsilon)
    eig <- eigen(K, symmetric = TRUE)$values
    expect_true(all(eig >= -1e-10),
                info = paste("Type:", name, "min eigenvalue:", min(eig)))
  }
})

# ---------------------------------------------------------------------------
# 3. kernel_roots reconstruction identity test: Khalf %*% Khalf = K
# ---------------------------------------------------------------------------

test_that("kernel_roots reconstruction satisfies Khalf %*% Khalf = K", {
  # Test across multiple kernel sizes and types
  test_cases <- list(
    # Simple diagonal
    list(K = diag(3), name = "identity"),
    # Random PSD
    list(K = {
      withr::local_seed(42)
      L <- matrix(rnorm(5 * 5), 5, 5)
      K <- crossprod(L)
      (K + t(K)) / 2
    }, name = "random_5x5"),
    # From design_kernel
    list(K = design_kernel(list(A = list(L = 4, type = "ordinal", l = 1.0)),
                           basis = "cell", normalize = "none")$K,
         name = "ordinal_kernel")
  )

  for (tc in test_cases) {
    roots <- kernel_roots(tc$K)
    reconstructed <- roots$Khalf %*% roots$Khalf

    expect_equal(reconstructed, tc$K, tolerance = 1e-8,
                 label = paste("Reconstruction for", tc$name))
  }
})

# ---------------------------------------------------------------------------
# 4. kernel_roots inverse reconstruction test: Kihalf %*% K %*% Kihalf = I
# ---------------------------------------------------------------------------

test_that("kernel_roots inverse satisfies Kihalf %*% K %*% Kihalf = I", {
  withr::local_seed(123)
  # Generate random well-conditioned PSD matrix
  L <- matrix(rnorm(4 * 4), 4, 4)
  K <- crossprod(L) + 0.1 * diag(4)  # Add ridge for conditioning
  K <- (K + t(K)) / 2

  roots <- kernel_roots(K)
  identity_approx <- roots$Kihalf %*% K %*% roots$Kihalf

  expect_equal(identity_approx, diag(4), tolerance = 1e-8)
})

# ---------------------------------------------------------------------------
# 5. kernel_roots preserves symmetry
# ---------------------------------------------------------------------------

test_that("kernel_roots output matrices are symmetric", {
  factors <- list(A = list(L = 3, type = "nominal"))
  K <- design_kernel(factors, basis = "cell", normalize = "none")$K

  roots <- kernel_roots(K)

  expect_true(isSymmetric(roots$Khalf, tol = 1e-10))
  expect_true(isSymmetric(roots$Kihalf, tol = 1e-10))
})

# ---------------------------------------------------------------------------
# 6. Multi-factor interaction kernel properties
# ---------------------------------------------------------------------------

test_that("multi-factor kernels maintain symmetry and PSD", {
  factors <- list(
    A = list(L = 3, type = "nominal"),
    B = list(L = 4, type = "ordinal", l = 1.5)
  )
  terms <- list("A", "B", c("A", "B"))

  K_res <- design_kernel(factors, terms = terms, basis = "cell", normalize = "unit_trace")
  K <- K_res$K

  # Symmetry
  expect_true(isSymmetric(K, tol = 1e-10))

  # PSD
  eig <- eigen(K, symmetric = TRUE)$values
  expect_true(all(eig >= -1e-10),
              info = paste("Min eigenvalue:", min(eig)))

  # Trace normalization check
  expect_equal(sum(diag(K)), 1, tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# 7. Randomized property test (multiple seeds)
# ---------------------------------------------------------------------------

test_that("kernel_roots is stable across random PSD matrices", {
  for (seed in c(1, 42, 999)) {
    withr::local_seed(seed)

    # Generate random PSD matrix
    q <- sample(3:10, 1)
    L <- matrix(rnorm(q * q), q, q)
    K <- crossprod(L)
    K <- (K + t(K)) / 2

    roots <- kernel_roots(K)

    # Reconstruction check
    reconstructed <- roots$Khalf %*% roots$Khalf
    expect_equal(reconstructed, K, tolerance = 1e-8,
                 info = paste("Seed:", seed, "q:", q))

    # Rank should equal q for full-rank matrices
    expect_equal(roots$rank, q,
                 info = paste("Seed:", seed, "q:", q))
  }
})

# ---------------------------------------------------------------------------
# 8. Circular distance symmetry test
# ---------------------------------------------------------------------------

test_that("circular kernel distances are symmetric around the period", {
  factors <- list(A = list(L = 6, type = "circular", l = 1.0))
  K <- design_kernel(factors, basis = "cell", normalize = "none")$K

  # For a circular factor with L=6, positions 1 and 5 have the same
  # circular distance from position 3 (both are distance 2).
  # Thus K[3,1] should equal K[3,5]
  expect_equal(K[3, 1], K[3, 5], tolerance = 1e-10)

  # Similarly, K[1,4] = distance 3 = K[1,4] (positions 1 and 4 have circular distance 3)
  # and K[1,3] = distance 2 = K[1,5] (positions 1 and 3, 1 and 5)
  expect_equal(K[1, 3], K[1, 5], tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# 9. Effect-basis kernel inherits symmetry and PSD
# ---------------------------------------------------------------------------

test_that("effect-basis kernel is symmetric and PSD", {
  factors <- list(
    A = list(L = 3, type = "nominal"),
    B = list(L = 4, type = "ordinal", l = 1.0)
  )
  Ctrs <- sum_contrasts(c(A = 3, B = 4))
  K_res <- design_kernel(factors, terms = list("A", "B", c("A", "B")),
                         basis = "effect", contrasts = Ctrs, normalize = "none")
  K <- K_res$K

  # Symmetry
  expect_true(isSymmetric(K, tol = 1e-10))

  # PSD
  eig <- eigen(K, symmetric = TRUE)$values
  expect_true(all(eig >= -1e-10),
              info = paste("Min eigenvalue:", min(eig)))
})

# ---------------------------------------------------------------------------
# 10. Cross-factor combinations preserve mathematical properties
# ---------------------------------------------------------------------------

test_that("three-factor kernels maintain symmetry and PSD", {
  factors <- list(
    A = list(L = 2, type = "nominal"),
    B = list(L = 3, type = "ordinal", l = 1.0),
    C = list(L = 2, type = "circular", l = 0.5)
  )

  K_res <- design_kernel(factors, basis = "cell", normalize = "none")
  K <- K_res$K

  # Expected dimensions: 2 * 3 * 2 = 12

  expect_equal(dim(K), c(12, 12))

  # Symmetry
  expect_true(isSymmetric(K, tol = 1e-10))

  # PSD
  eig <- eigen(K, symmetric = TRUE)$values
  expect_true(all(eig >= -1e-10))
})
