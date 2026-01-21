# test-s3-contracts.R
# S3 method contract verification for predict and as.data.frame methods

library(testthat)

# -----------------------------------------------------------------------------
# predict.dkge contracts
# -----------------------------------------------------------------------------

test_that("predict.dkge returns matrix with expected structure", {
  fixture <- make_small_fit()
  new_betas <- replicate(2, matrix(rnorm(fixture$q * 3), fixture$q, 3), simplify = FALSE)
  contrasts <- diag(fixture$q)

  result <- predict(fixture$fit, newdata = list(betas = new_betas, contrasts = contrasts))

  # Returns list with values component

  expect_true(is.list(result))
  expect_true("values" %in% names(result))

  # values is list of subjects
  expect_equal(length(result$values), length(new_betas))
})

test_that("predict.dkge returns matrix per subject with contrasts as columns", {
  fixture <- make_small_fit()
  new_betas <- replicate(2, matrix(rnorm(fixture$q * 3), fixture$q, 3), simplify = FALSE)
  n_contrasts <- 2
  contrasts <- diag(fixture$q)[, 1:n_contrasts, drop = FALSE]

  result <- predict(fixture$fit, B_list = new_betas, contrasts = contrasts)

  # Each subject's values is a matrix with n_clusters rows and n_contrasts columns
  expect_true(is.matrix(result$values[[1]]))
  expect_equal(ncol(result$values[[1]]), n_contrasts)
})

test_that("predict.dkge_model dispatches correctly", {
  fixture <- make_small_fit()
  model <- dkge_freeze(fixture$fit)
  new_betas <- fixture$betas
  contrasts <- diag(fixture$q)[, 1:2, drop = FALSE]

  # Verify predict generic dispatches to correct method
  expect_s3_class(model, "dkge_model")

  result <- predict(model, newdata = list(betas = new_betas, contrasts = contrasts))

  expect_true(is.list(result))
  expect_true("values" %in% names(result))
})

test_that("predict.dkge handles return_loadings flag", {
  fixture <- make_small_fit()
  new_betas <- replicate(2, matrix(rnorm(fixture$q * 3), fixture$q, 3), simplify = FALSE)
  contrasts <- diag(fixture$q)

  with_loadings <- predict(fixture$fit, B_list = new_betas, contrasts = contrasts,
                           return_loadings = TRUE)
  without_loadings <- predict(fixture$fit, B_list = new_betas, contrasts = contrasts,
                              return_loadings = FALSE)

  expect_true("A_list" %in% names(with_loadings))
  expect_false("A_list" %in% names(without_loadings))
})

# -----------------------------------------------------------------------------
# as.data.frame.dkge_contrasts contracts
# -----------------------------------------------------------------------------

test_that("as.data.frame.dkge_contrasts returns data.frame class", {
  fixture <- make_small_fit()
  contrasts <- diag(fixture$q)[, 1:2, drop = FALSE]
  contrast_obj <- dkge_contrast(fixture$fit, contrasts, method = "analytic",
                                parallel = FALSE, align = FALSE)

  df <- as.data.frame(contrast_obj)

  expect_s3_class(df, "data.frame")
})

test_that("as.data.frame.dkge_contrasts has expected columns", {
  fixture <- make_small_fit()
  contrasts <- diag(fixture$q)[, 1:2, drop = FALSE]
  contrast_obj <- dkge_contrast(fixture$fit, contrasts, method = "analytic",
                                parallel = FALSE, align = FALSE)

  df <- as.data.frame(contrast_obj)

  expected_columns <- c("contrast", "subject", "component", "value", "method")
  expect_true(all(expected_columns %in% names(df)))
})

test_that("as.data.frame.dkge_contrasts row count matches structure", {
  fixture <- make_small_fit()
  n_contrasts <- 2
  contrasts <- diag(fixture$q)[, 1:n_contrasts, drop = FALSE]
  contrast_obj <- dkge_contrast(fixture$fit, contrasts, method = "analytic",
                                parallel = FALSE, align = FALSE)

  df <- as.data.frame(contrast_obj)

  n_subjects <- length(contrast_obj$values[[1]])
  # Component count comes from actual values structure, not fit rank
  n_components <- length(contrast_obj$values[[1]][[1]])
  expected_rows <- n_contrasts * n_subjects * n_components

  expect_equal(nrow(df), expected_rows)
})

# -----------------------------------------------------------------------------
# as.data.frame.dkge_inference contracts
# -----------------------------------------------------------------------------

test_that("as.data.frame.dkge_inference returns data.frame class", {
  fixture <- make_small_fit()
  contrasts <- diag(fixture$q)[, 1:2, drop = FALSE]
  infer_obj <- dkge_infer(fixture$fit, contrasts,
                          method = "analytic",
                          inference = "parametric",
                          correction = "none",
                          n_perm = 10)

  df <- as.data.frame(infer_obj)

  expect_s3_class(df, "data.frame")
})

test_that("as.data.frame.dkge_inference has expected columns", {
  fixture <- make_small_fit()
  contrasts <- diag(fixture$q)[, 1:2, drop = FALSE]
  infer_obj <- dkge_infer(fixture$fit, contrasts,
                          method = "analytic",
                          inference = "parametric",
                          correction = "none",
                          n_perm = 10)

  df <- as.data.frame(infer_obj)

  expected_columns <- c("contrast", "component", "statistic", "p_value",
                        "p_adjusted", "significant", "alpha", "method",
                        "inference", "correction")
  expect_true(all(expected_columns %in% names(df)))
})

test_that("as.data.frame.dkge_inference has no NA statistics", {
  fixture <- make_small_fit()
  contrasts <- diag(fixture$q)[, 1:2, drop = FALSE]
  infer_obj <- dkge_infer(fixture$fit, contrasts,
                          method = "analytic",
                          inference = "parametric",
                          correction = "none",
                          n_perm = 10)

  df <- as.data.frame(infer_obj)

  expect_false(any(is.na(df$statistic)))
})

# -----------------------------------------------------------------------------
# as.data.frame.dkge_classification contracts
# -----------------------------------------------------------------------------

test_that("as.data.frame.dkge_classification returns data.frame class", {
  # Create minimal fixture
  obj <- structure(
    list(
      results = list(target1 = list(
        metrics = c(accuracy = 0.8, logloss = 0.3),
        p_values = c(accuracy = 0.05, logloss = 0.1),
        diagnostics = list(folds = list())
      )),
      method = "lda",
      metric = c("accuracy", "logloss"),
      n_perm = 10
    ),
    class = "dkge_classification"
  )

  df <- as.data.frame(obj, what = "summary")

  expect_s3_class(df, "data.frame")
})

test_that("as.data.frame.dkge_classification summary has expected columns", {
  obj <- structure(
    list(
      results = list(target1 = list(
        metrics = c(accuracy = 0.8),
        p_values = c(accuracy = 0.05),
        diagnostics = list(folds = list())
      )),
      method = "lda",
      metric = "accuracy",
      n_perm = 10
    ),
    class = "dkge_classification"
  )

  df <- as.data.frame(obj, what = "summary")

  expect_true("target" %in% names(df))
})
