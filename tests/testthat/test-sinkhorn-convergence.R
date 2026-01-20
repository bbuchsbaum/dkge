# test-sinkhorn-convergence.R
# Comprehensive tests for Sinkhorn optimal transport mathematical properties

library(testthat)

# =============================================================================
# Test 1: Doubly-stochastic verification (core requirement)
# =============================================================================

test_that("Sinkhorn plan is doubly-stochastic with uniform weights", {
  set.seed(42)

  # Create random cost matrix
  n <- 5
  m <- 6
  C <- matrix(runif(n * m, 0, 10), n, m)
  mu <- rep(1 / n, n)
  nu <- rep(1 / m, m)

  plan <- dkge:::.dkge_sinkhorn_plan(C, mu, nu, epsilon = 0.1, max_iter = 500, tol = 1e-7)

  # Row sums should equal mu within 1e-6 (same tolerance as Sinkhorn convergence)
  expect_equal(rowSums(plan), mu, tolerance = 1e-6)
  # Column sums should equal nu
  expect_equal(colSums(plan), nu, tolerance = 1e-6)
  # All entries non-negative
  expect_true(all(plan >= 0))
})

test_that("Sinkhorn plan is doubly-stochastic with non-uniform weights", {
  set.seed(42)

  n <- 4
  m <- 7
  C <- matrix(abs(rnorm(n * m, mean = 5, sd = 2)), n, m)
  mu <- runif(n)
  mu <- mu / sum(mu)
  nu <- runif(m)
  nu <- nu / sum(nu)

  plan <- dkge:::.dkge_sinkhorn_plan(C, mu, nu, epsilon = 0.05, max_iter = 500, tol = 1e-7)

  expect_equal(rowSums(plan), mu, tolerance = 1e-6)
  expect_equal(colSums(plan), nu, tolerance = 1e-6)
  expect_true(all(plan >= 0))
})

test_that("Sinkhorn plan is doubly-stochastic with square cost matrix", {
  set.seed(123)

  # Use 5x5 matrix for square matrix test
  n <- 5
  # Use bounded cost matrix to help convergence
  C <- matrix(runif(n * n, 0.1, 2), n, n)
  mu <- rep(1 / n, n)
  nu <- rep(1 / n, n)

  # Use larger epsilon (more entropic smoothing) for stable convergence
  plan <- dkge:::.dkge_sinkhorn_plan(C, mu, nu, epsilon = 0.1, max_iter = 500, tol = 1e-7)

  expect_equal(rowSums(plan), mu, tolerance = 1e-6)
  expect_equal(colSums(plan), nu, tolerance = 1e-6)
  expect_true(all(plan >= 0))
  expect_equal(sum(plan), 1, tolerance = 1e-6)  # Total mass should be 1
})

# =============================================================================
# Test 2: Identity transport for medoid
# =============================================================================

test_that("Medoid subject receives identity transport plan (medoid=1)", {
  set.seed(42)
  Q <- 4

  # Create simple test data with medoid = 1
  base_loadings <- matrix(rnorm(Q * 2), Q, 2)
  A_list <- list(
    base_loadings,
    base_loadings + matrix(rnorm(Q * 2, sd = 0.1), Q, 2),
    base_loadings + matrix(rnorm(Q * 2, sd = 0.1), Q, 2)
  )
  centroids <- list(
    matrix(runif(Q * 3), Q, 3),
    matrix(runif(Q * 3), Q, 3),
    matrix(runif(Q * 3), Q, 3)
  )
  v_list <- list(runif(Q), runif(Q), runif(Q))

  res <- dkge_transport_to_medoid_sinkhorn(
    v_list, A_list, centroids,
    medoid = 1,
    epsilon = 0.05,
    max_iter = 300,
    tol = 1e-6
  )

  # Medoid plan should be identity matrix
  expect_equal(res$plans[[1]], diag(1, Q))
})

test_that("Medoid subject receives identity transport plan (medoid=middle)", {
  set.seed(42)
  Q <- 5
  S <- 5
  medoid_idx <- 3

  A_list <- lapply(1:S, function(s) matrix(rnorm(Q * 2), Q, 2))
  centroids <- lapply(1:S, function(s) matrix(runif(Q * 3), Q, 3))
  v_list <- lapply(1:S, function(s) runif(Q))

  res <- dkge_transport_to_medoid_sinkhorn(
    v_list, A_list, centroids,
    medoid = medoid_idx,
    epsilon = 0.05,
    max_iter = 300,
    tol = 1e-6
  )

  # Medoid plan should be identity
  expect_equal(res$plans[[medoid_idx]], diag(1, Q))
})

test_that("Medoid subject receives identity transport plan (medoid=last)", {
  set.seed(42)
  Q <- 3
  S <- 4
  medoid_idx <- S

  A_list <- lapply(1:S, function(s) matrix(rnorm(Q * 2), Q, 2))
  centroids <- lapply(1:S, function(s) matrix(runif(Q * 3), Q, 3))
  v_list <- lapply(1:S, function(s) runif(Q))

  res <- dkge_transport_to_medoid_sinkhorn(
    v_list, A_list, centroids,
    medoid = medoid_idx,
    epsilon = 0.05,
    max_iter = 300,
    tol = 1e-6
  )

  expect_equal(res$plans[[medoid_idx]], diag(1, Q))
})

# =============================================================================
# Test 3: Identical embeddings produce near-diagonal transport
# =============================================================================

test_that("Identical embeddings produce near-diagonal transport with small epsilon", {
  set.seed(42)
  Q <- 5

  # Create identical loadings for source and target
  identical_loadings <- matrix(c(1, 0, 0, 1, 0.5, 0.5, 0.7, 0.3, 0.2, 0.8),
                                nrow = Q, byrow = TRUE)

  # Same loadings, same centroids
  A_list <- list(
    identical_loadings,
    identical_loadings  # Identical to medoid
  )
  centroids <- list(
    matrix(c(0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1), nrow = Q, byrow = TRUE),
    matrix(c(0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1), nrow = Q, byrow = TRUE)  # Same as medoid
  )
  v_list <- list(1:Q, 1:Q)

  res <- dkge_transport_to_medoid_sinkhorn(
    v_list, A_list, centroids,
    medoid = 1,
    epsilon = 1e-4,  # Small epsilon for near-deterministic
    max_iter = 2000,
    tol = 1e-9
  )

  # Subject 2's plan should be near-diagonal since embeddings are identical
  plan2 <- res$plans[[2]]

  # Diagonal entries should sum to > 0.9 (sparsity threshold from research)
  diag_sum <- sum(diag(plan2)) / sum(plan2)
  expect_gt(diag_sum, 0.9)
})

test_that("Similar embeddings produce transport with high sparsity", {
  set.seed(42)
  Q <- 4

  # Create very similar loadings (small perturbation)
  base_loadings <- matrix(c(1, 0, 0, 1, 1, 1, 0.5, 0.5), nrow = Q, byrow = TRUE)
  similar_loadings <- base_loadings + matrix(rnorm(Q * 2, sd = 0.01), Q, 2)

  # Same centroids for spatial alignment
  shared_centroids <- matrix(c(0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0), nrow = Q, byrow = TRUE)

  A_list <- list(base_loadings, similar_loadings)
  centroids <- list(shared_centroids, shared_centroids)
  v_list <- list(1:Q, 1:Q)

  res <- dkge_transport_to_medoid_sinkhorn(
    v_list, A_list, centroids,
    medoid = 1,
    epsilon = 1e-4,
    max_iter = 2000,
    tol = 1e-9
  )

  # Plan should be near-diagonal since embeddings are very similar
  plan2 <- res$plans[[2]]

  # Diagonal entries should dominate (> 80% of mass)
  diag_mass <- sum(diag(plan2))
  total_mass <- sum(plan2)
  expect_gt(diag_mass / total_mass, 0.8)
})

# =============================================================================
# Test 4: R/C++ equivalence
# =============================================================================

test_that("dkge_transport_to_medoid_sinkhorn and _cpp produce identical results", {
  set.seed(42)
  Q <- 5
  S <- 3

  A_list <- lapply(1:S, function(s) matrix(rnorm(Q * 2), Q, 2))
  centroids <- lapply(1:S, function(s) matrix(runif(Q * 3), Q, 3))
  v_list <- lapply(1:S, function(s) runif(Q))

  res_r <- dkge_transport_to_medoid_sinkhorn(
    v_list, A_list, centroids,
    medoid = 1,
    epsilon = 0.1,
    max_iter = 300,
    tol = 1e-8
  )

  res_cpp <- dkge_transport_to_medoid_sinkhorn_cpp(
    v_list, A_list, centroids,
    medoid = 1,
    epsilon = 0.1,
    max_iter = 300,
    tol = 1e-8,
    return_plans = TRUE
  )

  # Results should match within 1e-7 tolerance (minor floating point differences allowed)
  # Note: Both functions use the same C++ backend - differences arise from median computation
  # and floating-point accumulation order. Very small plan entries (1e-12 to 1e-16) have
  # relative errors that can exceed 1e-8 but absolute errors remain tiny.
  expect_equal(res_cpp$value, res_r$value, tolerance = 1e-7)
  expect_equal(res_cpp$subj_values, res_r$subj_values, tolerance = 1e-7)
  expect_equal(res_cpp$plans, res_r$plans, tolerance = 1e-7)
})

test_that("R and C++ produce identical results with non-uniform weights", {
  set.seed(123)
  Q <- 4
  S <- 4

  A_list <- lapply(1:S, function(s) matrix(rnorm(Q * 3), Q, 3))
  centroids <- lapply(1:S, function(s) matrix(runif(Q * 3), Q, 3))
  v_list <- lapply(1:S, function(s) runif(Q))
  sizes <- lapply(1:S, function(s) {
    sz <- runif(Q)
    sz / sum(sz)
  })

  res_r <- dkge_transport_to_medoid_sinkhorn(
    v_list, A_list, centroids,
    sizes = sizes,
    medoid = 2,
    epsilon = 0.05,
    max_iter = 500,
    tol = 1e-8
  )

  res_cpp <- dkge_transport_to_medoid_sinkhorn_cpp(
    v_list, A_list, centroids,
    sizes = sizes,
    medoid = 2,
    epsilon = 0.05,
    max_iter = 500,
    tol = 1e-8,
    return_plans = TRUE
  )

  expect_equal(res_cpp$value, res_r$value, tolerance = 1e-8)
  expect_equal(res_cpp$subj_values, res_r$subj_values, tolerance = 1e-8)
})

# =============================================================================
# Test 5: Edge cases
# =============================================================================

test_that("Very small epsilon converges (near-deterministic)", {
  set.seed(42)
  n <- 4
  m <- 4
  C <- matrix(c(0, 1, 2, 3,
                1, 0, 1, 2,
                2, 1, 0, 1,
                3, 2, 1, 0), n, m)
  mu <- rep(1 / n, n)
  nu <- rep(1 / m, m)

  # Very small epsilon - should still converge
  plan <- dkge:::.dkge_sinkhorn_plan(C, mu, nu, epsilon = 1e-5, max_iter = 5000, tol = 1e-6)

  # Should still be doubly stochastic
  expect_equal(rowSums(plan), mu, tolerance = 1e-5)  # Looser tolerance for small epsilon
  expect_equal(colSums(plan), nu, tolerance = 1e-5)
  expect_true(all(plan >= 0))

  # With small epsilon and diagonal-favoring cost, should be near-diagonal
  diag_fraction <- sum(diag(plan)) / sum(plan)
  expect_gt(diag_fraction, 0.8)
})

test_that("Non-convergence with very small max_iter returns best-effort result", {
  set.seed(42)
  n <- 6
  m <- 6
  C <- matrix(runif(n * m, 0, 10), n, m)
  mu <- rep(1 / n, n)
  nu <- rep(1 / m, m)

  # Very few iterations - won't fully converge
  plan <- dkge:::.dkge_sinkhorn_plan(C, mu, nu, epsilon = 0.01, max_iter = 5, tol = 1e-9)

  # Should still return a valid plan (no error)
  expect_true(is.matrix(plan))
  expect_equal(dim(plan), c(n, m))

  # All entries should be non-negative
  expect_true(all(plan >= 0))

  # Result should be "best effort" - doubly stochastic but potentially looser tolerance
  # The plan should still approximately sum to 1

  expect_equal(sum(plan), 1, tolerance = 0.1)
})

test_that("Mismatched mass totals produces error", {
  n <- 3
  m <- 4
  C <- matrix(1:12, n, m)
  mu <- c(0.2, 0.3, 0.5)  # sums to 1
  nu <- c(0.1, 0.2, 0.3, 0.5)  # sums to 1.1 - mismatch!

  expect_error(
    dkge:::.dkge_sinkhorn_plan(C, mu, nu, epsilon = 0.1),
    "mu and nu must sum to the same total mass"
  )
})

test_that("Zero or negative weights produce error", {
  n <- 3
  m <- 3
  C <- matrix(1:9, n, m)

  # Zero weight
  mu_zero <- c(0.5, 0, 0.5)
  nu <- c(1 / 3, 1 / 3, 1 / 3)

  expect_error(dkge:::.dkge_sinkhorn_plan(C, mu_zero, nu, epsilon = 0.1))

  # Negative weight
  mu_neg <- c(0.6, -0.1, 0.5)
  expect_error(dkge:::.dkge_sinkhorn_plan(C, mu_neg, nu, epsilon = 0.1))
})

# =============================================================================
# Test 6: Convergence failure detection and behavior
# =============================================================================

test_that("Non-convergence returns valid best-effort result without error", {
  set.seed(42)
  n <- 8
  m <- 8

  # Create a difficult problem: high dynamic range cost, very small epsilon
  # This combination often causes slow convergence
  C <- matrix(runif(n * m, 0, 100), n, m)
  mu <- rep(1 / n, n)
  nu <- rep(1 / m, m)

  # Use very few iterations with tight tolerance - guaranteed non-convergence
  plan <- dkge:::.dkge_sinkhorn_plan(C, mu, nu, epsilon = 0.001, max_iter = 3, tol = 1e-12)

  # The function should NOT error - it returns best-effort result
  expect_true(is.matrix(plan))
  expect_equal(dim(plan), c(n, m))

  # Result should still be a valid transport plan (non-negative)
  expect_true(all(plan >= 0))

  # Mass should be approximately preserved even without full convergence
  # (may have looser tolerance)
  expect_equal(sum(plan), 1, tolerance = 0.5)

  # Document behavior: current implementation returns last plan WITHOUT warning
  # This is accepted behavior per RESEARCH.md findings
})

test_that("Tighter tolerance achieves closer marginal match", {
  set.seed(42)
  n <- 5
  m <- 5
  C <- matrix(runif(n * m, 0, 5), n, m)
  mu <- rep(1 / n, n)
  nu <- rep(1 / n, n)

  # Run with looser tolerance (1e-4)
  plan_loose <- dkge:::.dkge_sinkhorn_plan(C, mu, nu, epsilon = 0.1, max_iter = 500, tol = 1e-4)
  err_loose_row <- max(abs(rowSums(plan_loose) - mu))
  err_loose_col <- max(abs(colSums(plan_loose) - nu))

  # Run with tighter tolerance (1e-8)
  plan_tight <- dkge:::.dkge_sinkhorn_plan(C, mu, nu, epsilon = 0.1, max_iter = 500, tol = 1e-8)
  err_tight_row <- max(abs(rowSums(plan_tight) - mu))
  err_tight_col <- max(abs(colSums(plan_tight) - nu))

  # Tighter tolerance should achieve smaller marginal errors
  expect_lt(err_tight_row, err_loose_row + 1e-5)
  expect_lt(err_tight_col, err_loose_col + 1e-5)

  # Tight tolerance should actually meet target
  expect_lt(err_tight_row, 1e-6)
  expect_lt(err_tight_col, 1e-6)
})

test_that("Looser tolerance converges with fewer iterations", {
  set.seed(42)
  n <- 4
  m <- 4
  C <- matrix(runif(n * m, 0, 3), n, m)
  mu <- rep(1 / n, n)
  nu <- rep(1 / m, m)

  # With loose tolerance (1e-3) and moderate iterations
  plan_loose <- dkge:::.dkge_sinkhorn_plan(C, mu, nu, epsilon = 0.1, max_iter = 50, tol = 1e-3)
  err_loose <- max(abs(rowSums(plan_loose) - mu), abs(colSums(plan_loose) - nu))

  # Loose tolerance should converge within 50 iterations (error < tolerance)
  expect_lt(err_loose, 1e-2)

  # Should be valid doubly-stochastic plan (to tolerance)
  expect_equal(sum(plan_loose), 1, tolerance = 0.01)
})

test_that("Tolerance affects marginal accuracy as expected", {
  set.seed(123)
  n <- 4
  m <- 4
  C <- matrix(runif(n * m, 0.5, 2), n, m)
  mu <- rep(1 / n, n)
  nu <- rep(1 / n, n)

  tolerances <- c(1e-3, 1e-5, 1e-7)
  errors <- numeric(length(tolerances))

  for (i in seq_along(tolerances)) {
    plan <- dkge:::.dkge_sinkhorn_plan(C, mu, nu, epsilon = 0.1, max_iter = 1000, tol = tolerances[i])
    errors[i] <- max(abs(rowSums(plan) - mu), abs(colSums(plan) - nu))
  }

  # Errors should decrease with tighter tolerance
  expect_true(errors[2] <= errors[1] + 1e-10)
  expect_true(errors[3] <= errors[2] + 1e-10)
})

# =============================================================================
# Test 7: Warm-start caching verification
# =============================================================================

test_that("Cache stores and retrieves dual variables", {
  skip_if_not(exists("sinkhorn_plan_cpp", envir = asNamespace("dkge"), inherits = FALSE))

  # Clear cache first
  dkge_clear_sinkhorn_cache()
  env <- dkge:::.dkge_sinkhorn_cache

  # Verify cache is empty
  cached_keys <- setdiff(ls(env, all.names = TRUE), ".order")
  expect_equal(length(cached_keys), 0)

  # Run Sinkhorn - should populate cache
  C <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9), 3, 3)
  mu <- rep(1 / 3, 3)
  nu <- rep(1 / 3, 3)
  plan1 <- dkge:::.dkge_sinkhorn_plan(C, mu, nu, epsilon = 0.1, max_iter = 100)

  # Verify cache now has entry
  cached_keys <- setdiff(ls(env, all.names = TRUE), ".order")
  expect_gt(length(cached_keys), 0)

  # Run same problem again - should use cached warm-start
  plan2 <- dkge:::.dkge_sinkhorn_plan(C, mu, nu, epsilon = 0.1, max_iter = 100)

  # Results should be identical (same inputs, same/warm-started algorithm)
  expect_equal(plan1, plan2, tolerance = 1e-12)

  # Cleanup
  dkge_clear_sinkhorn_cache()
})

test_that("Different inputs create separate cache entries", {
  skip_if_not(exists("sinkhorn_plan_cpp", envir = asNamespace("dkge"), inherits = FALSE))

  dkge_clear_sinkhorn_cache()
  env <- dkge:::.dkge_sinkhorn_cache

  # First problem
  C1 <- matrix(c(1, 2, 2, 1), 2, 2)
  mu1 <- c(0.5, 0.5)
  nu1 <- c(0.5, 0.5)
  dkge:::.dkge_sinkhorn_plan(C1, mu1, nu1, epsilon = 0.1, max_iter = 50)

  entries_after_first <- length(setdiff(ls(env, all.names = TRUE), ".order"))

  # Second problem (different cost matrix)
  C2 <- matrix(c(3, 1, 1, 3), 2, 2)
  dkge:::.dkge_sinkhorn_plan(C2, mu1, nu1, epsilon = 0.1, max_iter = 50)

  entries_after_second <- length(setdiff(ls(env, all.names = TRUE), ".order"))

  # Should have two separate cache entries
  expect_gt(entries_after_second, entries_after_first)

  # Cleanup
  dkge_clear_sinkhorn_cache()
})

test_that("Cache isolation between tests via dkge_clear_sinkhorn_cache", {
  skip_if_not(exists("sinkhorn_plan_cpp", envir = asNamespace("dkge"), inherits = FALSE))

  env <- dkge:::.dkge_sinkhorn_cache

  # Add some entries
  C <- matrix(runif(4), 2, 2)
  mu <- c(0.5, 0.5)
  nu <- c(0.5, 0.5)
  dkge:::.dkge_sinkhorn_plan(C, mu, nu, epsilon = 0.1, max_iter = 50)

  # Clear cache
  dkge_clear_sinkhorn_cache()

  # Verify empty
  cached_keys <- setdiff(ls(env, all.names = TRUE), ".order")
  expect_equal(length(cached_keys), 0)
  expect_equal(length(get(".order", envir = env, inherits = FALSE)), 0)
})
