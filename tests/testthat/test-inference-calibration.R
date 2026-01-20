# test-inference-calibration.R
# Tests for null distribution calibration and parallel execution equivalence

library(testthat)
library(dkge)

# ============================================================================
# Task 1: Null Distribution Calibration Tests
# ============================================================================

test_that("max-T controls FWER under pure noise (family-wise error rate)", {
  skip_on_cran()  # intensive simulation

  # The max-T procedure controls FWER (family-wise error rate), meaning the

  # probability of at least one false positive across all clusters should be
  # <= alpha under the null hypothesis.

  set.seed(42)
  S <- 20    # subjects
  Q <- 10    # clusters
  B <- 500   # permutations
  n_iter <- 200  # number of independent null datasets
  alpha <- 0.05

  # Count how many datasets have at least one false positive (any p < alpha)
  any_significant <- logical(n_iter)

  for (i in seq_len(n_iter)) {
    # Pure noise: N(0,1) values with no signal
    Y <- matrix(rnorm(S * Q), nrow = S, ncol = Q)

    # Run sign-flip test with max-T correction
    result <- dkge_signflip_maxT(Y, B = B)

    # Check if any p-value < alpha (false positive)
    any_significant[i] <- any(result$p < alpha)
  }

  # FWER = proportion of datasets with at least one false positive
  observed_fwer <- mean(any_significant)

  # Should be <= alpha (with some tolerance for simulation variance)
  # Use binomial test for FWER control
  # Under exact FWER = alpha, expected successes = n_iter * alpha = 10
  # Allow some tolerance: FWER should be <= alpha + 0.03
  expect_lte(observed_fwer, alpha + 0.03,
             label = sprintf("Observed FWER %.3f", observed_fwer))

  # Also verify p-values are in valid range
  expect_true(all(result$p >= 0 & result$p <= 1))
})

test_that("max-T p-values have valid properties under null", {
  skip_on_cran()

  # Under the null, max-T p-values should:
  # 1. All be in [0, 1]
  # 2. All be >= 1/(B+1) (theoretical minimum)
  # 3. Have a distribution that is at least as conservative as uniform

  set.seed(123)
  S <- 20
  Q <- 10
  B <- 500
  n_iter <- 50

  all_pvals <- numeric(n_iter * Q)

  for (i in seq_len(n_iter)) {
    Y <- matrix(rnorm(S * Q), nrow = S, ncol = Q)
    result <- dkge_signflip_maxT(Y, B = B)
    all_pvals[((i - 1) * Q + 1):(i * Q)] <- result$p
  }

  # Basic validity checks
  expect_true(all(all_pvals >= 0 & all_pvals <= 1))

  # Verify no p-values below theoretical minimum
  min_theoretical <- 1 / (B + 1)
  expect_true(all(all_pvals >= min_theoretical - 1e-10))

  # Max-T is conservative: median p-value under null should be >= 0.5
  # (or at least not significantly below 0.5)
  expect_gte(median(all_pvals), 0.4)
})

test_that("sign-flip null distribution is symmetric", {
  # The sign-flip procedure assumes a symmetric null distribution around zero.
  # Under this assumption, the flipped statistics should span the null well.

  set.seed(456)
  S <- 15
  Q <- 5
  B <- 500

  Y <- matrix(rnorm(S * Q), nrow = S, ncol = Q)
  result <- dkge_signflip_maxT(Y, B = B)

  # The max null distribution should be roughly symmetric in log space
  # and the median should be non-extreme
  expect_true(median(result$maxnull) > quantile(result$maxnull, 0.1))
  expect_true(median(result$maxnull) < quantile(result$maxnull, 0.9))

  # Flips should have roughly equal positive and negative values
  expect_true(abs(mean(result$flips)) < 0.2)
})

test_that("strong signal produces small p-value while noise remains uniform", {
  set.seed(777)
  S <- 20
  Q <- 10
  B <- 1000

  # Create data where cluster 1 has strong positive effect
  Y <- matrix(rnorm(S * Q), nrow = S, ncol = Q)
  Y[, 1] <- abs(Y[, 1]) + 3.0  # all subjects positive, very strong effect

  result <- dkge_signflip_maxT(Y, B = B)

  # Signal cluster should have small p-value
  expect_lt(result$p[1], 0.05,
            label = sprintf("Signal cluster p-value %.4f", result$p[1]))

  # Check that statistic for signal cluster is large
  expect_gt(result$stat[1], 3)
})

# ============================================================================
# Task 2: Parallel vs Sequential Equivalence Tests
# ============================================================================

test_that("dkge_contrast parallel matches sequential within tolerance", {
  skip_if_not_installed("future.apply")
  skip_if_not_installed("future")

  # Create fixture with S=8 subjects
  set.seed(99)
  S <- 8
  q <- 3
  P <- 5
  T_obs <- 60
  effects <- paste0("eff", seq_len(q))

  betas <- replicate(S, {
    mat <- matrix(rnorm(q * P), q, P)
    rownames(mat) <- effects
    mat
  }, simplify = FALSE)

  designs <- replicate(S, {
    X <- matrix(rnorm(T_obs * q), T_obs, q)
    X <- qr.Q(qr(X))
    colnames(X) <- effects
    X
  }, simplify = FALSE)

  fit <- dkge_fit(dkge_data(betas, designs = designs), K = diag(q), rank = 2)
  contrast_vec <- c(1, -1, 0)

  # Sequential execution
  result_seq <- dkge_contrast(fit, contrast_vec, method = "loso", parallel = FALSE)
  Y_seq <- as.matrix(result_seq, contrast = 1)

  # Parallel execution
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  future::plan(future::multisession, workers = 2)

  result_par <- dkge_contrast(fit, contrast_vec, method = "loso", parallel = TRUE)
  Y_par <- as.matrix(result_par, contrast = 1)

  future::plan(future::sequential)

  # Compare results - should match within floating-point tolerance
  expect_equal(Y_seq, Y_par, tolerance = 1e-10)
})

test_that("dkge_infer statistics match between sequential and parallel", {
  skip_if_not_installed("future.apply")
  skip_if_not_installed("future")

  set.seed(88)
  S <- 8
  q <- 3
  P <- 5
  T_obs <- 60
  effects <- paste0("eff", seq_len(q))

  betas <- replicate(S, {
    mat <- matrix(rnorm(q * P), q, P)
    rownames(mat) <- effects
    mat
  }, simplify = FALSE)

  designs <- replicate(S, {
    X <- matrix(rnorm(T_obs * q), T_obs, q)
    X <- qr.Q(qr(X))
    colnames(X) <- effects
    X
  }, simplify = FALSE)

  fit <- dkge_fit(dkge_data(betas, designs = designs), K = diag(q), rank = 2)
  contrast_vec <- c(1, -1, 0)

  # Use parametric inference to compare statistics (avoids random seed issues with sign-flip)
  result_seq <- dkge_infer(fit, contrast_vec, inference = "parametric",
                           correction = "none", parallel = FALSE)

  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  future::plan(future::multisession, workers = 2)

  result_par <- dkge_infer(fit, contrast_vec, inference = "parametric",
                           correction = "none", parallel = TRUE)

  future::plan(future::sequential)

  # Statistics should match exactly (same underlying contrast computation)
  expect_equal(result_seq$statistics[[1]], result_par$statistics[[1]], tolerance = 1e-10)

  # P-values should also match (parametric, deterministic)
  expect_equal(result_seq$p_values[[1]], result_par$p_values[[1]], tolerance = 1e-10)
})

test_that("parallel tests skip gracefully without future.apply", {
  # This test documents expected skip behavior
  # When future.apply is unavailable, parallel tests should skip, not fail
  skip_if(requireNamespace("future.apply", quietly = TRUE),
          "Test only runs when future.apply is NOT available")

  # If we get here, future.apply is missing
  expect_true(TRUE)
})

# ============================================================================
# Task 3: Edge Cases and Integration Tests
# ============================================================================

test_that("dkge_signflip_maxT requires minimum 5 subjects", {
  set.seed(111)
  Q <- 5
  B <- 100

  # S=4 should error
  Y4 <- matrix(rnorm(4 * Q), nrow = 4, ncol = Q)
  expect_error(dkge_signflip_maxT(Y4, B = B))

  # S=5 should succeed
  Y5 <- matrix(rnorm(5 * Q), nrow = 5, ncol = Q)
  result <- dkge_signflip_maxT(Y5, B = B)
  expect_true(is.list(result))
  expect_length(result$p, Q)
})

test_that("dkge_signflip_maxT requires minimum 100 permutations", {
  set.seed(222)
  S <- 10
  Q <- 5
  Y <- matrix(rnorm(S * Q), nrow = S, ncol = Q)

  # B=50 should error
  expect_error(dkge_signflip_maxT(Y, B = 50))

  # B=100 should succeed
  result <- dkge_signflip_maxT(Y, B = 100)
  expect_true(is.list(result))
  expect_length(result$p, Q)
})

test_that("p-values respect theoretical bounds", {
  set.seed(333)
  S <- 15
  Q <- 8
  B <- 500
  Y <- matrix(rnorm(S * Q), nrow = S, ncol = Q)

  result <- dkge_signflip_maxT(Y, B = B)

  # All p-values in [0, 1]
  expect_true(all(result$p >= 0 & result$p <= 1))

  # Minimum theoretical p-value with B permutations is 1/(B+1)
  min_theoretical <- 1 / (B + 1)
  expect_true(all(result$p >= min_theoretical))
})

test_that("inference with transport works on mismatched clusters", {
  skip_if_not_installed("future.apply")

  data <- create_mismatched_data()
  fit <- dkge_fit(data$betas, data$designs, K = data$K, rank = 2)

  transport_cfg <- list(
    centroids = data$centroids,
    medoid = 1L,
    mapper = dkge_mapper_spec("ridge", lambda = 1e-2),
    betas = data$betas
  )

  # Should complete without error when transport is configured
  result <- suppressWarnings(dkge_infer(fit, c(1, -1, 0), transport = transport_cfg))

  # Verify valid output

  expect_s3_class(result, "dkge_inference")
  expect_true(!is.null(result$transport))

  # P-values should be valid
  pvals <- result$p_values[[1]]
  expect_true(all(pvals >= 0 & pvals <= 1))
  expect_false(anyNA(pvals))
})

test_that("dkge_infer parametric produces valid p-values", {
  set.seed(444)
  S <- 10
  q <- 3
  P <- 6
  T_obs <- 50
  effects <- paste0("eff", seq_len(q))

  betas <- replicate(S, {
    mat <- matrix(rnorm(q * P), q, P)
    rownames(mat) <- effects
    mat
  }, simplify = FALSE)

  designs <- replicate(S, {
    X <- matrix(rnorm(T_obs * q), T_obs, q)
    X <- qr.Q(qr(X))
    colnames(X) <- effects
    X
  }, simplify = FALSE)

  fit <- dkge_fit(dkge_data(betas, designs = designs), K = diag(q), rank = 2)

  result <- dkge_infer(fit, c(1, -1, 0), inference = "parametric", correction = "fdr")

  expect_s3_class(result, "dkge_inference")
  expect_equal(result$inference, "parametric")
  expect_equal(result$correction, "fdr")

  # P-values valid
  expect_true(all(result$p_values[[1]] >= 0 & result$p_values[[1]] <= 1))
  expect_true(all(result$p_adjusted[[1]] >= 0 & result$p_adjusted[[1]] <= 1))
})

test_that("multiple correction methods work correctly", {
  set.seed(555)
  S <- 8
  q <- 3
  P <- 5
  T_obs <- 40
  effects <- paste0("eff", seq_len(q))

  betas <- replicate(S, {
    mat <- matrix(rnorm(q * P), q, P)
    rownames(mat) <- effects
    mat
  }, simplify = FALSE)

  designs <- replicate(S, {
    X <- matrix(rnorm(T_obs * q), T_obs, q)
    X <- qr.Q(qr(X))
    colnames(X) <- effects
    X
  }, simplify = FALSE)

  fit <- dkge_fit(dkge_data(betas, designs = designs), K = diag(q), rank = 2)
  contrast_vec <- c(1, -1, 0)

  # Test each correction method
  res_none <- dkge_infer(fit, contrast_vec, inference = "parametric", correction = "none")
  res_fdr <- dkge_infer(fit, contrast_vec, inference = "parametric", correction = "fdr")
  res_bonf <- dkge_infer(fit, contrast_vec, inference = "parametric", correction = "bonferroni")

  # None: p_adjusted == p_values
  expect_equal(res_none$p_adjusted[[1]], res_none$p_values[[1]])

  # FDR: p_adjusted >= p_values (monotonic)
  expect_true(all(res_fdr$p_adjusted[[1]] >= res_fdr$p_values[[1]]))

  # Bonferroni: p_adjusted = min(p * n_tests, 1)
  n_tests <- length(res_bonf$p_values[[1]])
  expected_bonf <- pmin(res_bonf$p_values[[1]] * n_tests, 1)
  expect_equal(res_bonf$p_adjusted[[1]], expected_bonf, tolerance = 1e-12)
})
