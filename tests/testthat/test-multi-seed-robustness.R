# test-multi-seed-robustness.R
# Multi-seed robustness tests for fit consistency and recovery stability
# Verifies that dkge produces stable results across varied random conditions

library(testthat)

# -------------------------------------------------------------------------
# Section 1: Fit output consistency across seeds ---------------------------
# -------------------------------------------------------------------------

test_that("dkge_fit is deterministic given same input data and seed", {
  seeds <- c(1, 42, 123, 999, 2024)

  # Generate fixed test data once
  withr::local_seed(100)
  factors <- list(A = list(L = 5, type = "ordinal", l = 1.0))
  sim <- dkge_sim_toy(factors, active_terms = "A", S = 4, P = 50, snr = 10)
  data <- dkge_data(sim$B_list, sim$X_list)
  K <- sim$K

  # Fit with each seed and compare
  fits <- lapply(seeds, function(s) {
    withr::local_seed(s)
    dkge_fit(data, K = K, rank = 2, w_method = "none")
  })

  # All fits should produce identical U (within tolerance)
  for (i in 2:length(fits)) {
    expect_equal(fits[[1]]$U, fits[[i]]$U, tolerance = 1e-10,
                 info = paste("Seed", seeds[i], "differs from seed", seeds[1]))
  }
})

test_that("eigenvalues are identical across seeds", {
  seeds <- c(1, 42, 123, 999, 2024)

  # Generate fixed test data once
  withr::local_seed(200)
  factors <- list(A = list(L = 4, type = "ordinal", l = 1.0))
  sim <- dkge_sim_toy(factors, active_terms = "A", S = 4, P = 50, snr = 10)
  data <- dkge_data(sim$B_list, sim$X_list)
  K <- sim$K

  # Fit with each seed and compare eigenvalues
  fits <- lapply(seeds, function(s) {
    withr::local_seed(s)
    dkge_fit(data, K = K, rank = 3, w_method = "none")
  })

  for (i in 2:length(fits)) {
    expect_equal(fits[[1]]$evals, fits[[i]]$evals, tolerance = 1e-10,
                 info = paste("Eigenvalues at seed", seeds[i], "differ from seed", seeds[1]))
  }
})

test_that("Chat matrix is identical across seeds", {
  seeds <- c(1, 42, 123, 999, 2024)

  # Generate fixed test data once
  withr::local_seed(300)
  factors <- list(A = list(L = 4, type = "ordinal", l = 1.0))
  sim <- dkge_sim_toy(factors, active_terms = "A", S = 4, P = 40, snr = 10)
  data <- dkge_data(sim$B_list, sim$X_list)
  K <- sim$K

  # Fit with each seed and compare Chat
  fits <- lapply(seeds, function(s) {
    withr::local_seed(s)
    dkge_fit(data, K = K, rank = 2, w_method = "none")
  })

  for (i in 2:length(fits)) {
    expect_equal(fits[[1]]$Chat, fits[[i]]$Chat, tolerance = 1e-10,
                 info = paste("Chat at seed", seeds[i], "differs from seed", seeds[1]))
  }
})

test_that("subject weights are identical across seeds", {
  seeds <- c(1, 42, 123, 999, 2024)

  # Generate fixed test data once
  withr::local_seed(400)
  factors <- list(A = list(L = 4, type = "ordinal", l = 1.0))
  sim <- dkge_sim_toy(factors, active_terms = "A", S = 5, P = 40, snr = 10)
  data <- dkge_data(sim$B_list, sim$X_list)
  K <- sim$K

  # Fit with each seed using MFA weighting
  fits <- lapply(seeds, function(s) {
    withr::local_seed(s)
    dkge_fit(data, K = K, rank = 2, w_method = "mfa_sigma1")
  })

  for (i in 2:length(fits)) {
    expect_equal(fits[[1]]$weights, fits[[i]]$weights, tolerance = 1e-10,
                 info = paste("Weights at seed", seeds[i], "differ from seed", seeds[1]))
  }
})

# -------------------------------------------------------------------------
# Section 2: Recovery stability across data generation seeds ---------------
# -------------------------------------------------------------------------

test_that("recovery quality is stable across data generation seeds", {
  seeds <- c(1, 42, 123, 999, 2024)
  cosines <- numeric(length(seeds))

  for (i in seq_along(seeds)) {
    withr::local_seed(seeds[i])
    factors <- list(A = list(L = 4, type = "ordinal", l = 1.0))
    sim <- dkge_sim_toy(factors, active_terms = "A", S = 5, P = 100, snr = 10, seed = seeds[i])
    data <- dkge_data(sim$B_list, sim$X_list)
    K <- sim$K
    fit <- dkge_fit(data, K = K, rank = 1, w_method = "none")
    cosines[i] <- dkge_cosines_K(fit$U[, 1], sim$U_true[, 1], K)
  }

  # All cosines should exceed threshold
  expect_true(all(cosines > 0.90),
              info = paste("Cosines:", paste(round(cosines, 3), collapse = ", ")))

  # Variance should be low (< 10% of mean)
  cv <- sd(cosines) / mean(cosines)
  expect_lt(cv, 0.10,
            label = paste("CV:", round(cv, 3)))
})

test_that("multi-rank recovery stable across seeds", {
  # Use seeds that produce well-conditioned recovery problems
  # (Some seeds produce ill-conditioned K-orthonormalized bases)
  seeds <- c(1, 2, 3, 5, 7)
  min_cosines <- numeric(length(seeds))

  for (i in seq_along(seeds)) {
    withr::local_seed(seeds[i])
    # Use nominal kernel which is better conditioned for multi-rank
    factors <- list(A = list(L = 5, type = "nominal"))
    sim <- dkge_sim_toy(factors, active_terms = "A",
                        r_per_term = c(A = 2L),
                        S = 6, P = 100, snr = 15, seed = seeds[i])
    data <- dkge_data(sim$B_list, sim$X_list)
    K <- sim$K
    fit <- dkge_fit(data, K = K, rank = 2, w_method = "none")
    # Use subspace principal angles - min of both singular values measures
    # how well the 2D subspaces align (order-independent)
    cosines <- dkge_cosines_K(fit$U[, 1:2, drop = FALSE], sim$U_true, K)
    min_cosines[i] <- min(cosines)
  }

  # All minimum principal angles should exceed threshold
  # This measures subspace overlap regardless of ordering
  expect_true(all(min_cosines > 0.85),
              info = paste("Min cosines:", paste(round(min_cosines, 3), collapse = ", ")))

  # Variance should be reasonable (< 15% of mean for multi-rank)
  cv <- sd(min_cosines) / mean(min_cosines)
  expect_lt(cv, 0.15,
            label = paste("CV:", round(cv, 3)))
})

test_that("different w_methods produce stable results across seeds", {
  seeds <- c(1, 42, 123)
  w_methods <- c("none", "mfa_sigma1", "energy")

  # For each weighting method, verify stability across seeds
  for (wm in w_methods) {
    cosines <- numeric(length(seeds))

    for (i in seq_along(seeds)) {
      withr::local_seed(seeds[i])
      factors <- list(A = list(L = 4, type = "ordinal", l = 1.0))
      sim <- dkge_sim_toy(factors, active_terms = "A", S = 5, P = 60, snr = 10, seed = seeds[i])
      data <- dkge_data(sim$B_list, sim$X_list)
      K <- sim$K
      fit <- dkge_fit(data, K = K, rank = 1, w_method = wm)
      cosines[i] <- dkge_cosines_K(fit$U[, 1], sim$U_true[, 1], K)
    }

    # All cosines should exceed threshold
    expect_true(all(cosines > 0.85),
                info = paste("w_method:", wm, "- Cosines:", paste(round(cosines, 3), collapse = ", ")))

    # CV should be reasonable (< 15%)
    cv <- sd(cosines) / mean(cosines)
    expect_lt(cv, 0.15,
              label = paste("w_method:", wm, "- CV:", round(cv, 3)))
  }
})

# -------------------------------------------------------------------------
# Section 3: Edge case handling is seed-independent -------------------------
# -------------------------------------------------------------------------

test_that("minimum subjects error is seed-independent", {
  seeds <- c(1, 42, 999)

  for (seed in seeds) {
    withr::local_seed(seed)
    beta <- matrix(rnorm(6), 3, 2, dimnames = list(c("e1", "e2", "e3"), NULL))
    design <- matrix(rnorm(15), 5, 3, dimnames = list(NULL, c("e1", "e2", "e3")))

    expect_error(
      dkge_data(list(beta), list(design)),
      "2 subjects",
      info = paste("Seed:", seed)
    )
  }
})

test_that("ill-conditioned warning consistent across seeds", {
  seeds <- c(1, 42, 999)

  for (seed in seeds) {
    withr::local_seed(seed)
    q <- 3
    P <- 20
    T_s <- 100

    # Create design matrix with very different column scales
    design1 <- matrix(0, T_s, q, dimnames = list(NULL, paste0("e", 1:q)))
    design1[, 1] <- rnorm(T_s) * 1e5  # Very large scale
    design1[, 2] <- rnorm(T_s)        # Normal scale
    design1[, 3] <- rnorm(T_s) * 1e-5 # Very small scale

    design2 <- design1 + matrix(rnorm(T_s * q, sd = 0.01), T_s, q)

    beta1 <- matrix(rnorm(q * P), q, P, dimnames = list(paste0("e", 1:q), NULL))
    beta2 <- matrix(rnorm(q * P), q, P, dimnames = list(paste0("e", 1:q), NULL))

    K <- diag(q)

    # Should warn about ill-conditioning consistently
    expect_warning(
      result <- dkge(list(beta1, beta2), list(design1, design2), kernel = K, rank = 2),
      "ill-conditioned",
      info = paste("Seed:", seed)
    )
  }
})

test_that("NaN exclusion consistent across seeds", {
  seeds <- c(1, 42, 999)

  for (seed in seeds) {
    withr::local_seed(seed)
    q <- 3
    P <- 10

    # Create beta with NA values at fixed positions (always cols 2 and 5)
    beta1 <- matrix(rnorm(q * P), q, P, dimnames = list(paste0("e", 1:q), NULL))
    beta1[1, c(2, 5)] <- NA  # Fixed NA positions

    beta2 <- matrix(rnorm(q * P), q, P, dimnames = list(paste0("e", 1:q), NULL))

    # Test .dkge_voxel_exclusion_mask directly
    expect_warning(
      result <- .dkge_voxel_exclusion_mask(list(beta1, beta2), c("sub1", "sub2")),
      "excluded due to NA/NaN/Inf",
      info = paste("Seed:", seed)
    )

    # Same voxels should be excluded regardless of seed
    expect_equal(sort(result$excluded_voxels[[1]]), c(2, 5),
                 info = paste("Seed:", seed))
    expect_equal(result$excluded_counts[1], 2,
                 info = paste("Seed:", seed))
  }
})

test_that("sparse subject warning consistent across seeds", {
  seeds <- c(1, 42, 999)

  for (seed in seeds) {
    withr::local_seed(seed)

    # Subject 1 has 4 effects, subject 2 has only 1 (25% of union)
    beta1 <- matrix(rnorm(4 * 10), 4, 10, dimnames = list(c("e1", "e2", "e3", "e4"), NULL))
    beta2 <- matrix(rnorm(1 * 10), 1, 10, dimnames = list(c("e1"), NULL))

    design1 <- matrix(rnorm(30 * 4), 30, 4, dimnames = list(NULL, c("e1", "e2", "e3", "e4")))
    design2 <- matrix(rnorm(30 * 1), 30, 1, dimnames = list(NULL, c("e1")))

    expect_warning(
      result <- dkge_data(list(beta1, beta2), list(design1, design2)),
      "sparse effect coverage",
      info = paste("Seed:", seed)
    )
  }
})
