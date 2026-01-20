# test-numerical-robustness.R
# Numerical edge case tests for degenerate inputs and graceful degradation

library(testthat)

# -------------------------------------------------------------------------
# Section 1: Rank-deficient inputs ----------------------------------------
# -------------------------------------------------------------------------

test_that("rank-deficient design matrix produces warning", {
  withr::local_seed(42)
  # Create design with collinear columns
  design <- matrix(c(1, 2, 3, 2, 4, 6), 3, 2)  # col2 = 2 * col1
  colnames(design) <- c("e1", "e2")
  beta <- matrix(rnorm(2 * 10), 2, 10, dimnames = list(c("e1", "e2"), NULL))

  expect_warning(
    dkge_subject(beta, design),
    "rank-deficient"
  )
})

test_that("rank-deficient beta matrix produces warning", {
  withr::local_seed(42)
  design <- matrix(rnorm(20 * 2), 20, 2, dimnames = list(NULL, c("e1", "e2")))
  # Create beta with linearly dependent rows (row2 = 2 * row1)
  beta <- matrix(c(1:5, 2:6), 2, 5, dimnames = list(c("e1", "e2"), NULL), byrow = TRUE)
  beta[2, ] <- 2 * beta[1, ]

  expect_warning(
    dkge_subject(beta, design),
    "reduced rank"
  )
})

test_that("dkge_fit completes with rank-deficient input after warning", {
  withr::local_seed(42)
  # Create well-conditioned case but with slightly rank-reduced structure
  q <- 3
  P <- 20
  T_s <- 50

  # Normal subject
  beta1 <- matrix(rnorm(q * P), q, P, dimnames = list(paste0("e", 1:q), NULL))
  design1 <- matrix(rnorm(T_s * q), T_s, q, dimnames = list(NULL, paste0("e", 1:q)))

  # Subject with collinear design columns
  design2 <- cbind(design1[, 1:2], design1[, 2])  # Third col = second col
  colnames(design2) <- paste0("e", 1:3)
  beta2 <- matrix(rnorm(q * P), q, P, dimnames = list(paste0("e", 1:q), NULL))

  K <- diag(q)  # Identity kernel

  # Should warn but complete
  expect_warning(
    result <- dkge(list(beta1, beta2), list(design1, design2), kernel = K, rank = 2),
    "rank-deficient"
  )

  expect_s3_class(result, "dkge")
  expect_true(!is.null(result$U))
})

test_that("effective_rank stored in fit object when reduced", {
  withr::local_seed(42)
  # Create a case where requested rank > effective rank
  q <- 5
  P <- 10
  T_s <- 30

  # Create betas with limited variation (only 2 effective dimensions)
  base1 <- rnorm(P)
  base2 <- rnorm(P)
  beta1 <- matrix(0, q, P, dimnames = list(paste0("e", 1:q), NULL))
  beta1[1, ] <- base1
  beta1[2, ] <- base2
  # Other rows are noise
  beta1[3:q, ] <- matrix(rnorm((q-2) * P, sd = 0.001), q-2, P)

  beta2 <- beta1 + matrix(rnorm(q * P, sd = 0.01), q, P)

  design1 <- matrix(rnorm(T_s * q), T_s, q, dimnames = list(NULL, paste0("e", 1:q)))
  design2 <- matrix(rnorm(T_s * q), T_s, q, dimnames = list(NULL, paste0("e", 1:q)))

  K <- diag(q)

  # Request high rank on low-dimensional data
  result <- suppressWarnings(
    dkge(list(beta1, beta2), list(design1, design2), kernel = K, rank = q)
  )

  expect_true(!is.null(result$effective_rank))
  expect_true(result$effective_rank > 0)
})

test_that("rank below minimum viable threshold produces appropriate response", {
  withr::local_seed(42)
  q <- 3
  P <- 5
  T_s <- 20

  # Create degenerate case with only 1 effective dimension
  base <- rnorm(P)
  beta1 <- matrix(0, q, P, dimnames = list(paste0("e", 1:q), NULL))
  beta1[1, ] <- base
  beta1[2:q, ] <- matrix(rnorm((q-1) * P, sd = 1e-10), q-1, P)

  beta2 <- beta1 + matrix(rnorm(q * P, sd = 1e-10), q, P)

  design1 <- matrix(rnorm(T_s * q), T_s, q, dimnames = list(NULL, paste0("e", 1:q)))
  design2 <- matrix(rnorm(T_s * q), T_s, q, dimnames = list(NULL, paste0("e", 1:q)))

  K <- diag(q)

  # Should still complete - verifying effective_rank < requested_rank
  result <- suppressWarnings(
    dkge(list(beta1, beta2), list(design1, design2), kernel = K, rank = q)
  )

  expect_s3_class(result, "dkge")
  expect_true(result$rank <= result$effective_rank)
})

# -------------------------------------------------------------------------
# Section 2: Ill-conditioned matrices -------------------------------------
# -------------------------------------------------------------------------

test_that("ill-conditioned pooled Gram produces warning", {
  withr::local_seed(42)
  q <- 3
  P <- 20
  T_s <- 100

  # Create design matrix with very different column scales
  design1 <- matrix(0, T_s, q, dimnames = list(NULL, paste0("e", 1:q)))
  design1[, 1] <- rnorm(T_s) * 1e5  # Very large scale
  design1[, 2] <- rnorm(T_s)         # Normal scale
  design1[, 3] <- rnorm(T_s) * 1e-5  # Very small scale

  design2 <- design1 + matrix(rnorm(T_s * q, sd = 0.01), T_s, q)

  beta1 <- matrix(rnorm(q * P), q, P, dimnames = list(paste0("e", 1:q), NULL))
  beta2 <- matrix(rnorm(q * P), q, P, dimnames = list(paste0("e", 1:q), NULL))

  K <- diag(q)

  # Should warn about ill-conditioning
  expect_warning(
    result <- dkge(list(beta1, beta2), list(design1, design2), kernel = K, rank = 2),
    "ill-conditioned"
  )
})

test_that("fit succeeds despite ill-conditioning warning", {
  withr::local_seed(42)
  q <- 3
  P <- 20
  T_s <- 100

  # Create moderately ill-conditioned design
  design1 <- matrix(0, T_s, q, dimnames = list(NULL, paste0("e", 1:q)))
  design1[, 1] <- rnorm(T_s) * 1e4
  design1[, 2] <- rnorm(T_s)
  design1[, 3] <- rnorm(T_s) * 1e-4

  design2 <- design1 + matrix(rnorm(T_s * q, sd = 0.01), T_s, q)

  beta1 <- matrix(rnorm(q * P), q, P, dimnames = list(paste0("e", 1:q), NULL))
  beta2 <- matrix(rnorm(q * P), q, P, dimnames = list(paste0("e", 1:q), NULL))

  K <- diag(q)

  # Should complete despite warning
  result <- suppressWarnings(
    dkge(list(beta1, beta2), list(design1, design2), kernel = K, rank = 2)
  )

  expect_s3_class(result, "dkge")
  expect_true(all(is.finite(result$U)))
})

test_that("well-conditioned matrices produce no ill-conditioning warning", {
  withr::local_seed(42)
  q <- 3
  P <- 20
  T_s <- 50

  # Create well-conditioned designs
  design1 <- matrix(rnorm(T_s * q), T_s, q, dimnames = list(NULL, paste0("e", 1:q)))
  design2 <- matrix(rnorm(T_s * q), T_s, q, dimnames = list(NULL, paste0("e", 1:q)))

  beta1 <- matrix(rnorm(q * P), q, P, dimnames = list(paste0("e", 1:q), NULL))
  beta2 <- matrix(rnorm(q * P), q, P, dimnames = list(paste0("e", 1:q), NULL))

  K <- diag(q)

  # Should NOT warn about ill-conditioning
  # Capture all warnings and verify none match "ill-conditioned"
  warnings_caught <- character(0)
  result <- withCallingHandlers(
    dkge(list(beta1, beta2), list(design1, design2), kernel = K, rank = 2),
    warning = function(w) {
      warnings_caught <<- c(warnings_caught, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  ill_cond_warnings <- grep("ill-conditioned", warnings_caught, value = TRUE)
  expect_length(ill_cond_warnings, 0)
})

# -------------------------------------------------------------------------
# Section 3: NaN/Inf handling ---------------------------------------------
# -------------------------------------------------------------------------

test_that("beta with NA values triggers exclusion warning", {
  withr::local_seed(42)
  q <- 3
  P <- 10

  # Create beta with NA values
  beta1 <- matrix(rnorm(q * P), q, P, dimnames = list(paste0("e", 1:q), NULL))
  beta1[1, 3] <- NA
  beta1[2, 5] <- NA

  beta2 <- matrix(rnorm(q * P), q, P, dimnames = list(paste0("e", 1:q), NULL))

  # Test .dkge_voxel_exclusion_mask directly
  expect_warning(
    result <- .dkge_voxel_exclusion_mask(list(beta1, beta2), c("sub1", "sub2")),
    "excluded due to NA/NaN/Inf"
  )

  expect_equal(result$excluded_counts[1], 2)  # Two voxels with NA in subject 1
  expect_equal(result$excluded_counts[2], 0)  # No NA in subject 2
  expect_equal(result$total_excluded, 2)
})

test_that("beta with Inf values triggers exclusion warning", {
  withr::local_seed(42)
  q <- 3
  P <- 10

  # Create beta with Inf values
  beta1 <- matrix(rnorm(q * P), q, P, dimnames = list(paste0("e", 1:q), NULL))
  beta1[1, 2] <- Inf
  beta1[2, 7] <- -Inf

  beta2 <- matrix(rnorm(q * P), q, P, dimnames = list(paste0("e", 1:q), NULL))

  # Test .dkge_voxel_exclusion_mask directly
  expect_warning(
    result <- .dkge_voxel_exclusion_mask(list(beta1, beta2), c("sub1", "sub2")),
    "excluded due to NA/NaN/Inf"
  )

  expect_equal(result$excluded_counts[1], 2)  # Two voxels with Inf in subject 1
  expect_equal(result$total_excluded, 2)
})

test_that("exclusion metadata correctly identifies affected voxels", {
  withr::local_seed(42)
  q <- 2
  P <- 8

  # Create beta with specific NA positions
  beta1 <- matrix(1:(q * P), q, P, dimnames = list(c("e1", "e2"), NULL))
  beta1[1, c(2, 5)] <- NA  # Columns 2 and 5 should be excluded

  beta2 <- matrix(1:(q * P), q, P, dimnames = list(c("e1", "e2"), NULL))

  result <- suppressWarnings(
    .dkge_voxel_exclusion_mask(list(beta1, beta2))
  )

  expect_equal(sort(result$excluded_voxels[[1]]), c(2, 5))
  expect_equal(result$excluded_voxels[[2]], integer(0))
})

test_that("all-NaN subject triggers exclusion warning with correct count", {
  withr::local_seed(42)
  q <- 2
  P <- 5

  # Create all-NaN beta
  beta_all_na <- matrix(NA_real_, q, P, dimnames = list(c("e1", "e2"), NULL))
  beta_normal <- matrix(rnorm(q * P), q, P, dimnames = list(c("e1", "e2"), NULL))

  expect_warning(
    result <- .dkge_voxel_exclusion_mask(list(beta_all_na, beta_normal), c("na_sub", "normal_sub")),
    "excluded due to NA/NaN/Inf"
  )

  # All P voxels in subject 1 should be excluded
  expect_equal(result$excluded_counts[1], P)
  expect_equal(result$excluded_counts[2], 0)
})

# -------------------------------------------------------------------------
# Section 4: Partial effect overlap ---------------------------------------
# -------------------------------------------------------------------------

test_that("sparse subject (>50% missing) produces warning", {
  withr::local_seed(42)

  # Subject 1 has 4 effects, subject 2 has only 1 (25% of union)
  beta1 <- matrix(rnorm(4 * 10), 4, 10, dimnames = list(c("e1", "e2", "e3", "e4"), NULL))
  beta2 <- matrix(rnorm(1 * 10), 1, 10, dimnames = list(c("e1"), NULL))

  design1 <- matrix(rnorm(30 * 4), 30, 4, dimnames = list(NULL, c("e1", "e2", "e3", "e4")))
  design2 <- matrix(rnorm(30 * 1), 30, 1, dimnames = list(NULL, c("e1")))

  expect_warning(
    result <- dkge_data(list(beta1, beta2), list(design1, design2)),
    "sparse effect coverage"
  )
})

test_that("fit succeeds with partial overlap", {
  withr::local_seed(42)

  # Subject 1 has effects e1, e2, e3
  # Subject 2 has effects e2, e3, e4
  # Overlap: e2, e3
  beta1 <- matrix(rnorm(3 * 10), 3, 10, dimnames = list(c("e1", "e2", "e3"), NULL))
  beta2 <- matrix(rnorm(3 * 10), 3, 10, dimnames = list(c("e2", "e3", "e4"), NULL))

  design1 <- matrix(rnorm(30 * 3), 30, 3, dimnames = list(NULL, c("e1", "e2", "e3")))
  design2 <- matrix(rnorm(30 * 3), 30, 3, dimnames = list(NULL, c("e2", "e3", "e4")))

  K <- diag(4)  # Identity kernel for 4 effects in union

  result <- suppressWarnings(
    dkge(list(beta1, beta2), list(design1, design2), kernel = K, rank = 2)
  )

  expect_s3_class(result, "dkge")
  expect_equal(length(result$effects), 4)  # Union has 4 effects
})

test_that("disjoint effects work with zero cross-pairs", {
  withr::local_seed(42)

  # Subject 1 has effects e1, e2 (no overlap with subject 2)
  # Subject 2 has effects e3, e4
  beta1 <- matrix(rnorm(2 * 10), 2, 10, dimnames = list(c("e1", "e2"), NULL))
  beta2 <- matrix(rnorm(2 * 10), 2, 10, dimnames = list(c("e3", "e4"), NULL))

  design1 <- matrix(rnorm(30 * 2), 30, 2, dimnames = list(NULL, c("e1", "e2")))
  design2 <- matrix(rnorm(30 * 2), 30, 2, dimnames = list(NULL, c("e3", "e4")))

  # Should produce sparse subject warnings (each has 50% of effects)
  data <- suppressWarnings(
    dkge_data(list(beta1, beta2), list(design1, design2))
  )

  expect_s3_class(data, "dkge_data")
  expect_equal(length(data$effects), 4)

  # Check provenance for zero cross-pairs
  expect_equal(data$provenance$pair_counts["e1", "e3"], 0L)
  expect_equal(data$provenance$pair_counts["e2", "e4"], 0L)
})

# -------------------------------------------------------------------------
# Section 5: Minimum subjects ---------------------------------------------
# -------------------------------------------------------------------------

test_that("single subject errors with informative message", {
  withr::local_seed(42)
  beta <- matrix(rnorm(3 * 10), 3, 10, dimnames = list(c("e1", "e2", "e3"), NULL))
  design <- matrix(rnorm(20 * 3), 20, 3, dimnames = list(NULL, c("e1", "e2", "e3")))

  expect_error(
    dkge_data(list(beta), list(design)),
    "At least 2 subjects required"
  )
})

test_that("two subjects is minimum viable", {
  withr::local_seed(42)
  beta1 <- matrix(rnorm(3 * 10), 3, 10, dimnames = list(c("e1", "e2", "e3"), NULL))
  beta2 <- matrix(rnorm(3 * 10), 3, 10, dimnames = list(c("e1", "e2", "e3"), NULL))
  design <- matrix(rnorm(20 * 3), 20, 3, dimnames = list(NULL, c("e1", "e2", "e3")))

  result <- dkge_data(list(beta1, beta2), list(design, design))

  expect_s3_class(result, "dkge_data")
  expect_equal(result$n_subjects, 2)
})

# -------------------------------------------------------------------------
# Additional robustness tests ---------------------------------------------
# -------------------------------------------------------------------------

test_that("check_rank helper returns correct metadata", {
  # Well-formed matrix
  design_good <- matrix(rnorm(20 * 3), 20, 3)
  result_good <- .dkge_check_rank(design_good)
  expect_equal(result_good$design_rank, 3)
  expect_null(result_good$beta_rank)

  # With beta
  beta_good <- matrix(rnorm(3 * 10), 3, 10)
  result_both <- .dkge_check_rank(design_good, beta_good)
  expect_equal(result_both$design_rank, 3)
  expect_equal(result_both$beta_rank, 3)
})

test_that("check_condition returns condition number", {
  withr::local_seed(42)
  M_good <- crossprod(matrix(rnorm(100 * 3), 100, 3))
  cond <- .dkge_check_condition(M_good, threshold = 1e8, name = "test")
  expect_true(is.numeric(cond))
  expect_true(cond > 0)
})
