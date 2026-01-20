library(testthat)

set.seed(1)

test_that("DKGE recovers a single main-effect component and contrasts line up", {
  skip_on_cran()
  factors <- list(A = list(L = 2, type = "nominal"))
  sim <- dkge_sim_toy(factors, active_terms = "A", S = 3, P = 12, snr = 10)

  fit <- dkge_fit(sim$B_list,
                  designs = sim$X_list,
                  K = sim$K,
                  w_method = "none",
                  rank = 1)

  cos1 <- dkge_cosines_K(sim$U_true[, 1, drop = FALSE],
                         fit$U[, 1, drop = FALSE],
                         sim$K)
  expect_gt(cos1[1], 0.98)

  c_true <- as.numeric(sim$U_true[, 1])
  cres <- dkge_contrast(fit, c_true, method = "loso", align = FALSE)
  cors <- vapply(seq_along(sim$B_list), function(s) {
    v_s <- cres$values[[1]][[s]]
    stats::cor(v_s, sim$M_list[[s]][, 1])
  }, numeric(1))
  expect_true(all(abs(cors) > 0.9))
})

# -------------------------------------------------------------------------
# SNR-level recovery tests
# -------------------------------------------------------------------------

test_that("Recovery achieves high cosine (> 0.98) at very high SNR (20)", {
  skip_on_cran()
  withr::local_seed(2001)
  factors <- list(A = list(L = 3, type = "nominal"))
  sim <- dkge_sim_toy(factors, active_terms = "A", S = 5, P = 30, snr = 20, seed = 2001)

  fit <- dkge_fit(sim$B_list,
                  designs = sim$X_list,
                  K = sim$K,
                  w_method = "none",
                  rank = 1)

  cosines <- dkge_cosines_K(sim$U_true[, 1, drop = FALSE],
                            fit$U[, 1, drop = FALSE],
                            sim$K)
  expect_gt(cosines[1], 0.98,
            label = paste("High SNR (20) recovery cosine:", round(cosines[1], 4)))
})

test_that("Recovery achieves good cosine (> 0.90) at medium SNR (8)", {
  skip_on_cran()
  withr::local_seed(2002)
  factors <- list(A = list(L = 3, type = "nominal"))
  sim <- dkge_sim_toy(factors, active_terms = "A", S = 5, P = 30, snr = 8, seed = 2002)

  fit <- dkge_fit(sim$B_list,
                  designs = sim$X_list,
                  K = sim$K,
                  w_method = "none",
                  rank = 1)

  cosines <- dkge_cosines_K(sim$U_true[, 1, drop = FALSE],
                            fit$U[, 1, drop = FALSE],
                            sim$K)
  expect_gt(cosines[1], 0.90,
            label = paste("Medium SNR (8) recovery cosine:", round(cosines[1], 4)))
})

test_that("Recovery achieves acceptable cosine (> 0.70) at low SNR (2)", {
  skip_on_cran()
  withr::local_seed(2003)
  factors <- list(A = list(L = 3, type = "nominal"))
  # Use more subjects and clusters for stability at low SNR
  sim <- dkge_sim_toy(factors, active_terms = "A", S = 8, P = 50, snr = 2, seed = 2003)

  fit <- dkge_fit(sim$B_list,
                  designs = sim$X_list,
                  K = sim$K,
                  w_method = "none",
                  rank = 1)

  cosines <- dkge_cosines_K(sim$U_true[, 1, drop = FALSE],
                            fit$U[, 1, drop = FALSE],
                            sim$K)
  expect_gt(cosines[1], 0.70,
            label = paste("Low SNR (2) recovery cosine:", round(cosines[1], 4)))
})

test_that("Recovery quality degrades gracefully with decreasing SNR", {
  skip_on_cran()
  withr::local_seed(2004)
  factors <- list(A = list(L = 3, type = "nominal"))

  snr_levels <- c(20, 10, 5, 2)
  cosines <- numeric(length(snr_levels))

  for (i in seq_along(snr_levels)) {
    sim <- dkge_sim_toy(factors, active_terms = "A", S = 5, P = 30,
                        snr = snr_levels[i], seed = 2004 + i)
    fit <- dkge_fit(sim$B_list, designs = sim$X_list, K = sim$K,
                    w_method = "none", rank = 1)
    cosines[i] <- dkge_cosines_K(sim$U_true[, 1, drop = FALSE],
                                  fit$U[, 1, drop = FALSE],
                                  sim$K)[1]
  }

  # Verify monotonic decrease (with small tolerance for sampling variance)
  for (i in 2:length(snr_levels)) {
    expect_gte(cosines[i - 1], cosines[i] - 0.05,
               label = paste("Cosine at SNR", snr_levels[i - 1],
                            "should be >= cosine at SNR", snr_levels[i]))
  }

  # Verify bounds at extremes

  expect_gt(cosines[1], 0.95, label = "Highest SNR should have high cosine")
  expect_gt(cosines[length(cosines)], 0.50, label = "Even lowest SNR should recover some signal")
})

# -------------------------------------------------------------------------
# Multi-rank recovery tests
# -------------------------------------------------------------------------

test_that("Multi-rank recovery (r=2) recovers both components", {
  skip_on_cran()
  withr::local_seed(3001)
  # Use a factor with enough levels to have 2 effects
  factors <- list(A = list(L = 4, type = "nominal"))
  # Active two effect columns from factor A
  sim <- dkge_sim_toy(factors, active_terms = "A",
                      r_per_term = c(A = 2L),
                      S = 5, P = 30, snr = 15, seed = 3001)

  fit <- dkge_fit(sim$B_list,
                  designs = sim$X_list,
                  K = sim$K,
                  w_method = "none",
                  rank = 2)

  # Check that both components are recovered
  r_true <- ncol(sim$U_true)
  expect_equal(r_true, 2)

  cosines <- dkge_cosines_K(sim$U_true,
                            fit$U[, 1:r_true, drop = FALSE],
                            sim$K)

  # Both cosines should exceed threshold
  expect_true(all(cosines > 0.90),
              label = paste("Multi-rank (r=2) cosines:", paste(round(cosines, 4), collapse = ", ")))
})

test_that("Multi-rank recovery (r=3) recovers all components", {
  skip_on_cran()
  withr::local_seed(3002)
  # Use a factor with enough levels to have 3 effects
  factors <- list(A = list(L = 5, type = "nominal"))
  sim <- dkge_sim_toy(factors, active_terms = "A",
                      r_per_term = c(A = 3L),
                      S = 6, P = 40, snr = 15, seed = 3002)

  fit <- dkge_fit(sim$B_list,
                  designs = sim$X_list,
                  K = sim$K,
                  w_method = "none",
                  rank = 3)

  r_true <- ncol(sim$U_true)
  expect_equal(r_true, 3)

  cosines <- dkge_cosines_K(sim$U_true,
                            fit$U[, 1:r_true, drop = FALSE],
                            sim$K)

  # All three cosines should exceed threshold
  expect_true(all(cosines > 0.85),
              label = paste("Multi-rank (r=3) cosines:", paste(round(cosines, 4), collapse = ", ")))
})

test_that("Multi-factor recovery works with two nominal factors", {
  skip_on_cran()
  withr::local_seed(3003)
  factors <- list(
    A = list(L = 3, type = "nominal"),
    B = list(L = 2, type = "nominal")
  )
  # Activate both A and B main effects (1 column each)
  sim <- dkge_sim_toy(factors, active_terms = c("A", "B"),
                      r_per_term = c(A = 1L, B = 1L),
                      S = 5, P = 30, snr = 12, seed = 3003)

  fit <- dkge_fit(sim$B_list,
                  designs = sim$X_list,
                  K = sim$K,
                  w_method = "none",
                  rank = 2)

  r_true <- ncol(sim$U_true)
  expect_equal(r_true, 2)

  cosines <- dkge_cosines_K(sim$U_true,
                            fit$U[, 1:r_true, drop = FALSE],
                            sim$K)

  # Both factor components should be recovered
  expect_true(all(cosines > 0.85),
              label = paste("Multi-factor cosines:", paste(round(cosines, 4), collapse = ", ")))
})

test_that("Multi-factor with interaction recovers all components", {
  skip_on_cran()
  withr::local_seed(3004)
  factors <- list(
    A = list(L = 3, type = "nominal"),
    B = list(L = 2, type = "nominal")
  )
  # Activate A, B, and A:B interaction
  sim <- dkge_sim_toy(factors, active_terms = c("A", "B", "A:B"),
                      r_per_term = c(A = 1L, B = 1L, "A:B" = 1L),
                      S = 6, P = 40, snr = 15, seed = 3004)

  fit <- dkge_fit(sim$B_list,
                  designs = sim$X_list,
                  K = sim$K,
                  w_method = "none",
                  rank = 3)

  r_true <- ncol(sim$U_true)
  expect_equal(r_true, 3)

  cosines <- dkge_cosines_K(sim$U_true,
                            fit$U[, 1:r_true, drop = FALSE],
                            sim$K)

  # All components including interaction should be recovered
  expect_true(all(cosines > 0.80),
              label = paste("Multi-factor+interaction cosines:", paste(round(cosines, 4), collapse = ", ")))
})
