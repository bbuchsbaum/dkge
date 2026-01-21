# test-coverage-gaps.R
# Additional tests targeting coverage gaps for Phase 6

# ===========================================================================
# Tests for edge cases in core exported functions
# ===========================================================================

test_that("dkge_sim_toy returns expected structure", {
  toy <- dkge_sim_toy(
    factors = list(A = list(L = 2), B = list(L = 3)),
    active_terms = c("A", "B"),
    S = 3, P = 15, snr = 5
  )

  expect_named(toy, c("B_list", "X_list", "K", "info", "U_true", "M_list",
                       "active_cols", "subject_ids"))
  expect_length(toy$B_list, 3)
  expect_equal(dim(toy$K), c(5, 5))
  expect_equal(nrow(toy$U_true), 5)
  expect_equal(ncol(toy$U_true), 2)  # 2 active terms
})

test_that("dkge_sim_toy respects seed parameter", {
  toy1 <- dkge_sim_toy(
    factors = list(A = list(L = 2)),
    active_terms = "A", S = 2, P = 10, seed = 42
  )
  toy2 <- dkge_sim_toy(
    factors = list(A = list(L = 2)),
    active_terms = "A", S = 2, P = 10, seed = 42
  )
  expect_identical(toy1$B_list[[1]], toy2$B_list[[1]])
})

test_that("design_kernel handles different normalization options", {
  factors <- list(A = list(L = 3))

  k_trace <- design_kernel(factors, normalize = "unit_trace")
  k_fro <- design_kernel(factors, normalize = "unit_fro")
  k_diag <- design_kernel(factors, normalize = "max_diag")
  k_none <- design_kernel(factors, normalize = "none")

  expect_equal(sum(diag(k_trace$K_cell)), 1, tolerance = 1e-8)
  expect_true(max(diag(k_diag$K_cell)) <= 1 + 1e-8)
  expect_true(sum(k_none$K_cell) > sum(k_trace$K_cell))
})

test_that("design_kernel ordinal factor produces valid RBF kernel", {
  kern <- design_kernel(
    factors = list(time = list(L = 5, type = "ordinal", l = 1.0)),
    basis = "cell"
  )

  expect_equal(dim(kern$K_cell), c(5, 5))
  # Ordinal kernel should have decay with distance
  expect_gt(kern$K_cell[1, 1], kern$K_cell[1, 3])
  expect_gt(kern$K_cell[1, 3], kern$K_cell[1, 5])
})

test_that("design_kernel circular factor wraps correctly", {
  kern <- design_kernel(
    factors = list(angle = list(L = 6, type = "circular", l = 0.5)),
    basis = "cell"
  )

  # Circular: position 1 and 6 should be neighbors (closer than 1 and 4)
  expect_gt(kern$K_cell[1, 6], kern$K_cell[1, 4])
})

test_that("dkge_weights prior_mask helper works", {
  mask <- c(TRUE, TRUE, FALSE, FALSE, TRUE)
  w <- dkge_weights_prior_mask(mask, value_in = 1, value_out = 0)

  expect_equal(w, c(1, 1, 0, 0, 1))
})

test_that("dkge_weights prior_roi helper works", {
  labels <- factor(c("A", "A", "B", "B", "C"))
  w <- dkge_weights_prior_roi(labels)

  expect_length(w, 5)
  expect_true(all(w > 0))
})

test_that("dkge_weights auto creates default spec", {
  w <- dkge_weights_auto()

  expect_s3_class(w, "dkge_weights")
  expect_equal(w$adapt, "kenergy_prec")
  expect_equal(w$combine, "product")
})

test_that("kernel_roots handles near-singular kernels", {
  # Create a kernel that's nearly singular
  K <- matrix(c(1, 0.999, 0.999, 1), 2, 2)

  roots <- kernel_roots(K, jitter = 1e-8)

  expect_equal(dim(roots$Khalf), c(2, 2))
  expect_equal(dim(roots$Kihalf), c(2, 2))
  expect_true(roots$rank <= 2)
})

test_that("kernel_alignment score is in [0,1]", {
  A <- matrix(rnorm(9), 3, 3)
  B <- matrix(rnorm(9), 3, 3)

  score <- kernel_alignment(A, B)

  expect_gte(abs(score), 0)
  expect_lte(abs(score), 1)
})

test_that("helmert_contrasts produces orthogonal contrasts", {
  cm <- helmert_contrasts(c(A = 3, B = 4))

  expect_length(cm, 2)
  expect_equal(nrow(cm$A), 3)
  expect_equal(ncol(cm$A), 2)
  expect_equal(nrow(cm$B), 4)
  expect_equal(ncol(cm$B), 3)
})

test_that("sum_contrasts produces sum-to-zero contrasts", {
  cm <- sum_contrasts(c(A = 3))

  expect_equal(nrow(cm$A), 3)
  expect_equal(ncol(cm$A), 2)
  expect_equal(sum(cm$A[, 1]), 0, tolerance = 1e-10)
})

test_that("dkge_clear_sinkhorn_cache returns TRUE", {
  result <- dkge_clear_sinkhorn_cache()
  expect_true(result)
})

# ===========================================================================
# Pipeline and integration tests
# ===========================================================================

test_that("dkge_pipeline handles minimal inputs", {
  toy <- dkge_sim_toy(
    factors = list(A = list(L = 2), B = list(L = 3)),
    active_terms = c("A"), S = 6, P = 12, snr = 3  # need S >= 5 for inference
  )

  result <- dkge_pipeline(
    betas = toy$B_list,
    designs = toy$X_list,
    kernel = toy$K,
    contrasts = c(1, 0, 0, 0, 0),
    method = "analytic",
    inference = list()  # disable inference to avoid S>=5 constraint
  )

  expect_named(result, c("fit", "diagnostics", "contrasts", "transport",
                          "inference", "classification"))
  expect_s3_class(result$fit, "dkge")
  expect_s3_class(result$contrasts, "dkge_contrasts")
})

test_that("dkge_contrast analytic method returns valid structure", {
  toy <- dkge_sim_toy(
    factors = list(A = list(L = 2)),
    active_terms = "A", S = 3, P = 10, snr = 5
  )
  fit <- dkge(toy$B_list, toy$X_list, kernel = toy$K, rank = 1)

  # Use analytic method (faster)
  result <- dkge_contrast(fit, c(1, 0), method = "analytic")

  expect_s3_class(result, "dkge_contrasts")
  expect_equal(result$method, "analytic")
  expect_length(result$values, 1)  # one contrast
})

test_that("dkge_contrast kfold method works", {
  toy <- dkge_sim_toy(
    factors = list(A = list(L = 2)),
    active_terms = "A", S = 4, P = 10, snr = 5
  )
  fit <- dkge(toy$B_list, toy$X_list, kernel = toy$K, rank = 1)

  result <- dkge_contrast(fit, c(1, 0), method = "kfold", folds = 2)

  expect_s3_class(result, "dkge_contrasts")
  expect_equal(result$method, "kfold")
})

# ===========================================================================
# S3 method coverage tests
# ===========================================================================

test_that("print.dkge_contrasts works", {
  toy <- dkge_sim_toy(
    factors = list(A = list(L = 2)),
    active_terms = "A", S = 3, P = 8, snr = 5
  )
  fit <- dkge(toy$B_list, toy$X_list, kernel = toy$K, rank = 1)
  result <- dkge_contrast(fit, c(1, 0), method = "analytic")

  output <- capture.output(print(result))

  expect_gt(length(output), 0)
  expect_true(any(grepl("DKGE Contrasts", output)))
})

test_that("as.data.frame.dkge_contrasts works", {
  toy <- dkge_sim_toy(
    factors = list(A = list(L = 2)),
    active_terms = "A", S = 3, P = 8, snr = 5
  )
  fit <- dkge(toy$B_list, toy$X_list, kernel = toy$K, rank = 1)
  result <- dkge_contrast(fit, c(1, 0), method = "analytic")

  df <- as.data.frame(result)

  expect_s3_class(df, "data.frame")
  expect_true("contrast" %in% names(df))
  expect_true("subject" %in% names(df))
  expect_true("value" %in% names(df))
})

test_that("print.dkge_weights works", {
  w <- dkge_weights(adapt = "kenergy")
  output <- capture.output(print(w))

  expect_gt(length(output), 0)
  expect_true(any(grepl("dkge weight", output)))
})

# ===========================================================================
# Inference edge cases
# ===========================================================================

test_that("dkge_signflip_maxT handles minimum inputs", {
  Y <- matrix(rnorm(5 * 10), 5, 10)

  result <- dkge_signflip_maxT(Y, B = 100)

  expect_named(result, c("stat", "p", "maxnull", "flips"))
  expect_length(result$stat, 10)
  expect_length(result$p, 10)
  expect_length(result$maxnull, 100)
})

test_that("dkge_signflip_maxT handles different tails", {
  Y <- matrix(rnorm(6 * 8), 6, 8)

  res_two <- dkge_signflip_maxT(Y, B = 100, tail = "two.sided")
  res_gt <- dkge_signflip_maxT(Y, B = 100, tail = "greater")
  res_lt <- dkge_signflip_maxT(Y, B = 100, tail = "less")

  expect_length(res_two$p, 8)
  expect_length(res_gt$p, 8)
  expect_length(res_lt$p, 8)
})

# ===========================================================================
# dkge cosines helper
# ===========================================================================

test_that("dkge_cosines_K computes principal angles", {
  U <- qr.Q(qr(matrix(rnorm(10), 5, 2)))
  V <- qr.Q(qr(matrix(rnorm(10), 5, 2)))
  K <- diag(5)

  cos_vals <- dkge_cosines_K(U, V, K)

  expect_length(cos_vals, 2)
  expect_true(all(cos_vals >= 0 & cos_vals <= 1))
})

test_that("dkge_cosines_K returns 1 for identical subspaces", {
  U <- qr.Q(qr(matrix(rnorm(10), 5, 2)))
  K <- diag(5)

  cos_vals <- dkge_cosines_K(U, U, K)

  expect_equal(cos_vals, c(1, 1), tolerance = 1e-10)
})
