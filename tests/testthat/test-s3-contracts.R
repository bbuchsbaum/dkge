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

# -----------------------------------------------------------------------------
# S3 method dispatch registration tests
# -----------------------------------------------------------------------------

test_that("print methods are registered in S3 method table", {
  # Get all registered print methods for dkge classes
  print_methods <- methods(print)
  dkge_print_methods <- print_methods[grepl("^print\\.dkge_", print_methods)]

  # Verify key print methods exist
  expected_methods <- c(
    "print.dkge_classification",
    "print.dkge_classification_spec",
    "print.dkge_contrast_validated",
    "print.dkge_contrasts",
    "print.dkge_folds",
    "print.dkge_inference",
    "print.dkge_inference_spec",
    "print.dkge_regress",
    "print.dkge_transport_spec",
    "print.dkge_weights"
  )

  for (method in expected_methods) {
    expect_true(method %in% dkge_print_methods,
                info = sprintf("Method %s should be registered", method))
  }
})

test_that("predict methods are registered in S3 method table", {
  predict_methods <- methods(predict)
  dkge_predict_methods <- predict_methods[grepl("^predict\\.dkge", predict_methods)]

  expected_methods <- c("predict.dkge", "predict.dkge_model")

  for (method in expected_methods) {
    expect_true(method %in% dkge_predict_methods,
                info = sprintf("Method %s should be registered", method))
  }
})

test_that("as.data.frame methods are registered in S3 method table", {
  adf_methods <- methods(as.data.frame)
  dkge_adf_methods <- adf_methods[grepl("^as\\.data\\.frame\\.dkge", adf_methods)]

  expected_methods <- c(
    "as.data.frame.dkge_classification",
    "as.data.frame.dkge_contrasts",
    "as.data.frame.dkge_inference"
  )

  for (method in expected_methods) {
    expect_true(method %in% dkge_adf_methods,
                info = sprintf("Method %s should be registered", method))
  }
})

test_that("generic dispatch for print works correctly", {
  # Create objects and verify dispatch by output type

  # dkge_contrasts
  fixture <- make_small_fit()
  contrasts <- diag(fixture$q)[, 1:2, drop = FALSE]
  contrast_obj <- dkge_contrast(fixture$fit, contrasts, method = "analytic",
                                parallel = FALSE, align = FALSE)
  output <- capture.output(result <- print(contrast_obj))
  expect_identical(result, contrast_obj)

  # dkge_weights
  weights_obj <- dkge_weights(adapt = "none")
  output <- capture.output(result <- print(weights_obj))
  expect_identical(result, weights_obj)

  # dkge_inference_spec
  spec_obj <- dkge_inference_spec(B = 100)
  output <- capture.output(result <- print(spec_obj))
  expect_identical(result, spec_obj)
})

test_that("generic dispatch for predict works correctly", {
  fixture <- make_small_fit()
  new_betas <- fixture$betas
  contrasts <- diag(fixture$q)

  # predict.dkge dispatch
  expect_s3_class(fixture$fit, "dkge")
  result <- predict(fixture$fit, B_list = new_betas, contrasts = contrasts)
  expect_true(is.list(result))
  expect_true("values" %in% names(result))

  # predict.dkge_model dispatch
  model <- dkge_freeze(fixture$fit)
  expect_s3_class(model, "dkge_model")
  result <- predict(model, newdata = list(betas = new_betas, contrasts = contrasts))
  expect_true(is.list(result))
})

test_that("generic dispatch for as.data.frame works correctly", {
  fixture <- make_small_fit()
  contrasts <- diag(fixture$q)[, 1:2, drop = FALSE]

  # as.data.frame.dkge_contrasts dispatch
  contrast_obj <- dkge_contrast(fixture$fit, contrasts, method = "analytic",
                                parallel = FALSE, align = FALSE)
  expect_s3_class(contrast_obj, "dkge_contrasts")
  df <- as.data.frame(contrast_obj)
  expect_s3_class(df, "data.frame")

  # as.data.frame.dkge_inference dispatch
  infer_obj <- dkge_infer(fixture$fit, contrasts, method = "analytic",
                          inference = "parametric", correction = "none", n_perm = 10)
  expect_s3_class(infer_obj, "dkge_inference")
  df <- as.data.frame(infer_obj)
  expect_s3_class(df, "data.frame")
})

test_that("class-specific dispatch works for dkge_model", {
  # dkge_model is a separate class (frozen fit) with its own predict method
  fixture <- make_small_fit()
  model <- dkge_freeze(fixture$fit)

  # Model has dkge_model class
  expect_true("dkge_model" %in% class(model))

  # Verify dispatch to predict.dkge_model works correctly
  new_betas <- fixture$betas
  contrasts <- diag(fixture$q)
  result <- predict(model, newdata = list(betas = new_betas, contrasts = contrasts))
  expect_true(is.list(result))
  expect_true("values" %in% names(result))
})
