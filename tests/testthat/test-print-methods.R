# test-print-methods.R
# Tests for S3 print method contracts - all print methods return invisibly

library(testthat)

# -----------------------------------------------------------------------------
# Helper to create minimal fixture objects for each class
# -----------------------------------------------------------------------------

make_classification_fixture <- function() {
  fixture <- make_small_fit()
  # Create a minimal dkge_classification object manually
  structure(
    list(
      results = list(target1 = list(
        metrics = c(accuracy = 0.8),
        p_values = c(accuracy = 0.05)
      )),
      method = "lda",
      metric = "accuracy",
      n_perm = 10
    ),
    class = "dkge_classification"
  )
}

make_classification_spec_fixture <- function() {
  dkge_classification_spec(targets = ~ group, method = "lda")
}

make_contrast_validated_fixture <- function() {
  # Create minimal structure mimicking dkge_contrast_validated
  structure(
    list(
      observed = list(values = list(c1 = list(s1 = 1:3))),
      completed = list(values = list(c1 = list(s1 = 1:3))),
      summary = data.frame(
        contrast = "c1",
        component = 1,
        observed_mean = 0.5,
        completed_mean = 0.6,
        stringsAsFactors = FALSE
      ),
      provenance = list()
    ),
    class = "dkge_contrast_validated"
  )
}

make_contrasts_fixture <- function() {
  fixture <- make_small_fit()
  contrasts <- diag(fixture$q)[, 1:2, drop = FALSE]
  dkge_contrast(fixture$fit, contrasts, method = "analytic",
                parallel = FALSE, align = FALSE)
}

make_folds_fixture <- function() {
  fixture <- make_small_fit()
  # Create minimal dkge_folds manually (dkge_define_folds needs a fit with subjects)
  structure(
    list(
      type = "subject",
      k = 2,
      assignments = list(1L, 2:3),
      metadata = list(n_subjects = 3, fold_sizes = c(1, 2))
    ),
    class = "dkge_folds"
  )
}

make_inference_fixture <- function() {
  fixture <- make_small_fit()
  contrasts <- diag(fixture$q)[, 1:2, drop = FALSE]
  dkge_infer(fixture$fit, contrasts,
             method = "analytic",
             inference = "parametric",
             correction = "none",
             n_perm = 10)
}

make_inference_spec_fixture <- function() {
  dkge_inference_spec(B = 500, tail = "two.sided")
}

make_regress_fixture <- function() {
  # Create minimal dkge_regress object manually
  structure(
    list(
      metrics = list(
        rmse_micro = 0.5,
        r2_micro = 0.8,
        r2_macro = 0.7,
        cosine_mean = 0.9
      ),
      pred = matrix(rnorm(6), 3, 2),
      truth = matrix(rnorm(6), 3, 2)
    ),
    class = "dkge_regress"
  )
}

make_transport_spec_fixture <- function() {
  dkge_transport_spec(
    centroids = list(matrix(runif(12), 4, 3), matrix(runif(12), 4, 3)),
    medoid = 1L
  )
}

make_weights_fixture <- function() {
  dkge_weights(adapt = "none")
}

# -----------------------------------------------------------------------------
# Print method contract tests
# -----------------------------------------------------------------------------

test_that("print.dkge_classification returns invisibly", {
  obj <- make_classification_fixture()
  output <- capture.output(result <- print(obj))
  expect_s3_class(result, "dkge_classification")
  expect_identical(result, obj)
  expect_true(length(output) > 0)
})

test_that("print.dkge_classification_spec returns invisibly", {
  obj <- make_classification_spec_fixture()
  output <- capture.output(result <- print(obj))
  expect_s3_class(result, "dkge_classification_spec")
  expect_identical(result, obj)
  expect_true(length(output) > 0)
})

test_that("print.dkge_contrast_validated returns invisibly", {
  obj <- make_contrast_validated_fixture()
  output <- capture.output(result <- print(obj))
  expect_s3_class(result, "dkge_contrast_validated")
  expect_identical(result, obj)
  expect_true(length(output) > 0)
})

test_that("print.dkge_contrasts returns invisibly", {
  obj <- make_contrasts_fixture()
  output <- capture.output(result <- print(obj))
  expect_s3_class(result, "dkge_contrasts")
  expect_identical(result, obj)
  expect_true(length(output) > 0)
})

test_that("print.dkge_folds returns invisibly", {
  obj <- make_folds_fixture()
  output <- capture.output(result <- print(obj))
  expect_s3_class(result, "dkge_folds")
  expect_identical(result, obj)
  expect_true(length(output) > 0)
})

test_that("print.dkge_inference returns invisibly", {
  obj <- make_inference_fixture()
  output <- capture.output(result <- print(obj))
  expect_s3_class(result, "dkge_inference")
  expect_identical(result, obj)
  expect_true(length(output) > 0)
})

test_that("print.dkge_inference_spec returns invisibly", {
  obj <- make_inference_spec_fixture()
  output <- capture.output(result <- print(obj))
  expect_s3_class(result, "dkge_inference_spec")
  expect_identical(result, obj)
  expect_true(length(output) > 0)
})

test_that("print.dkge_regress returns invisibly", {
  obj <- make_regress_fixture()
  output <- capture.output(result <- print(obj))
  expect_s3_class(result, "dkge_regress")
  expect_identical(result, obj)
  expect_true(length(output) > 0)
})

test_that("print.dkge_transport_spec returns invisibly", {
  obj <- make_transport_spec_fixture()
  output <- capture.output(result <- print(obj))
  expect_s3_class(result, "dkge_transport_spec")
  expect_identical(result, obj)
  expect_true(length(output) > 0)
})

test_that("print.dkge_weights returns invisibly", {
  obj <- make_weights_fixture()
  output <- capture.output(result <- print(obj))
  expect_s3_class(result, "dkge_weights")
  expect_identical(result, obj)
  expect_true(length(output) > 0)
})

# -----------------------------------------------------------------------------
# Verify print output contains expected content
# -----------------------------------------------------------------------------

test_that("print.dkge_contrasts includes method and counts", {
  obj <- make_contrasts_fixture()
  output <- capture.output(print(obj))
  output_text <- paste(output, collapse = "\n")
  expect_match(output_text, "Method:", fixed = TRUE)
  expect_match(output_text, "Contrasts:", fixed = TRUE)
  expect_match(output_text, "Subjects:", fixed = TRUE)
})

test_that("print.dkge_inference includes method and correction", {
  obj <- make_inference_fixture()
  output <- capture.output(print(obj))
  output_text <- paste(output, collapse = "\n")
  expect_match(output_text, "Cross-fitting:", fixed = TRUE)
  expect_match(output_text, "Inference:", fixed = TRUE)
  expect_match(output_text, "Correction:", fixed = TRUE)
})

test_that("print.dkge_folds includes type and fold count", {
  obj <- make_folds_fixture()
  output <- capture.output(print(obj))
  output_text <- paste(output, collapse = "\n")
  expect_match(output_text, "Type:", fixed = TRUE)
  expect_match(output_text, "Folds:", fixed = TRUE)
})

test_that("print.dkge_weights includes key fields", {
  obj <- make_weights_fixture()
  output <- capture.output(print(obj))
  output_text <- paste(output, collapse = "\n")
  expect_match(output_text, "prior", fixed = TRUE)
  expect_match(output_text, "adapt", fixed = TRUE)
  expect_match(output_text, "combine", fixed = TRUE)
})
