# test-integration-pipeline.R
# Comprehensive integration tests for the dkge_pipeline() end-to-end workflow

library(testthat)

# ---- Local helper functions ----

#' Create multi-subject synthetic data for integration tests
#'
#' @param S Number of subjects
#' @param q Number of design effects
#' @param P Number of clusters/voxels per subject
#' @param T Number of timepoints per subject design matrix
#' @param seed Random seed for reproducibility
make_integration_data <- function(S = 6, q = 5, P = 20, T = 24, seed = 123) {
  set.seed(seed)
  betas <- replicate(S, matrix(rnorm(q * P), q, P), simplify = FALSE)
  designs <- replicate(S, {
    X <- matrix(rnorm(T * q), T, q)
    colnames(X) <- paste0("eff", seq_len(q))
    qr.Q(qr(X))
  }, simplify = FALSE)
  K <- diag(q)
  list(betas = betas, designs = designs, K = K, S = S, q = q, P = P)
}

#' Create centroids for transport tests with matched dimensions
make_integration_centroids <- function(S, P, n_features = 3, seed = 456) {
  set.seed(seed)
  replicate(S, matrix(rnorm(P * n_features), P, n_features), simplify = FALSE)
}

# ---- Task 1: Comprehensive pipeline integration tests ----

test_that("dkge_pipeline completes full workflow with valid multi-subject data", {
  # Integration test: Full pipeline from scratch
  dat <- make_integration_data(S = 6, q = 5, P = 20)
  cvec <- c(1, -1, 0, 0, 0)

  res <- dkge_pipeline(
    betas = dat$betas,
    designs = dat$designs,
    kernel = dat$K,
    contrasts = cvec,
    rank = 3
  )

  # Verify all expected components present

  expect_true(all(c("fit", "diagnostics", "contrasts") %in% names(res)))
  expect_s3_class(res$fit, "dkge")
  expect_type(res$diagnostics, "list")
  expect_s3_class(res$contrasts, "dkge_contrasts")

  # Verify fit metadata
  expect_equal(res$diagnostics$rank, 3)
  expect_equal(length(res$fit$Btil), 6)
})

test_that("dkge_pipeline chains fit -> LOSO contrast -> transport -> inference correctly", {
  dat <- make_integration_data(S = 6, q = 5, P = 20)
  centroids <- make_integration_centroids(S = 6, P = 20)

  # Pre-fit so we can test pipeline with existing fit
  fit <- dkge_fit(dat$betas, dat$designs, K = dat$K, rank = 3)

  transport_cfg <- list(
    centroids = centroids,
    medoid = 1L
  )

  inference_cfg <- list(B = 100)

  cvec <- c(1, -1, 0, 0, 0)

  res <- suppressWarnings(dkge_pipeline(
    fit = fit,
    contrasts = cvec,
    transport = transport_cfg,
    inference = inference_cfg
  ))

  # Verify full chain
  expect_s3_class(res$fit, "dkge")
  expect_s3_class(res$contrasts, "dkge_contrasts")
  expect_false(is.null(res$transport))
  expect_false(is.null(res$inference))

  # Verify transport produced results for each subject
  expect_equal(length(res$transport), 1)  # one contrast
  expect_false(is.null(res$transport[[1]]$subj_values))

  # Verify inference produced results
  expect_equal(length(res$inference), 1)  # one contrast
  expect_false(is.null(res$inference[[1]]$stat))
  expect_false(is.null(res$inference[[1]]$p))
})

test_that("dkge_pipeline handles both pre-computed fit and fit-from-scratch modes", {
  dat <- make_integration_data(S = 5, q = 4, P = 15)
  cvec <- c(1, -1, 0, 0)

  # Mode 1: Fit from scratch (betas, designs, kernel provided)
  res_scratch <- dkge_pipeline(
    betas = dat$betas,
    designs = dat$designs,
    kernel = dat$K,
    contrasts = cvec,
    rank = 2
  )

  expect_s3_class(res_scratch$fit, "dkge")
  expect_equal(res_scratch$diagnostics$rank, 2)

  # Mode 2: Pre-computed fit
  fit <- dkge_fit(dat$betas, dat$designs, K = dat$K, rank = 2)
  res_precomputed <- dkge_pipeline(
    fit = fit,
    contrasts = cvec
  )

  expect_s3_class(res_precomputed$fit, "dkge")
  expect_equal(res_precomputed$diagnostics$rank, 2)

  # Both modes produce valid contrast results
  expect_s3_class(res_scratch$contrasts, "dkge_contrasts")
  expect_s3_class(res_precomputed$contrasts, "dkge_contrasts")
})

test_that("LOSO method produces cross-fitted results with U_minus bases", {
  dat <- make_integration_data(S = 5, q = 4, P = 15)
  cvec <- c(1, -1, 0, 0)

  res <- dkge_pipeline(
    betas = dat$betas,
    designs = dat$designs,
    kernel = dat$K,
    contrasts = cvec,
    method = "loso",
    rank = 2
  )

  # LOSO should populate bases in metadata
  expect_false(is.null(res$contrasts$metadata$bases))
  bases <- res$contrasts$metadata$bases

  # Should have S fold bases (one per subject in LOSO)
  expect_equal(length(bases), 5)

  # Each basis should be qxr
  for (basis in bases) {
    expect_equal(nrow(basis), 4)  # q
    expect_equal(ncol(basis), 2)  # rank
  }
})

test_that("dkge_pipeline errors on missing required arguments", {
  dat <- make_integration_data(S = 5, q = 4, P = 15)

  cvec <- c(1, -1, 0, 0)

  # Error: No fit AND no betas/designs/kernel (stopifnot produces specific message)
  expect_error(
    dkge_pipeline(contrasts = cvec),
    regexp = "is\\.null\\(betas\\)|betas"
  )

  # Error: Invalid contrast dimensions (contrast length != q)
  wrong_cvec <- c(1, -1, 0)  # 3 elements when q=4
  fit <- dkge_fit(dat$betas, dat$designs, K = dat$K, rank = 2)
  expect_error(
    dkge_pipeline(fit = fit, contrasts = wrong_cvec),
    regexp = "length|contrast"
  )
})

test_that("dkge_pipeline handles mismatched subject counts gracefully",
{
  # Create data with subjects having different cluster counts
  set.seed(789)
  P_vec <- c(15, 18, 12, 20, 16)
  q <- 4
  S <- 5
  T <- 20

  betas <- lapply(seq_len(S), function(s) {
    matrix(rnorm(q * P_vec[s]), q, P_vec[s])
  })
  designs <- replicate(S, {
    X <- matrix(rnorm(T * q), T, q)
    colnames(X) <- paste0("eff", seq_len(q))
    qr.Q(qr(X))
  }, simplify = FALSE)
  K <- diag(q)

  centroids <- lapply(seq_len(S), function(s) {
    matrix(rnorm(P_vec[s] * 3), P_vec[s], 3)
  })

  fit <- dkge_fit(betas, designs, K = K, rank = 2)
  cvec <- c(1, -1, 0, 0)

  # With transport to medoid, should handle varying dimensions
  transport_cfg <- list(
    centroids = centroids,
    medoid = 1L,
    betas = betas
  )

  res <- suppressWarnings(dkge_pipeline(
    fit = fit,
    contrasts = cvec,
    transport = transport_cfg
  ))

  expect_s3_class(res$contrasts, "dkge_contrasts")
  expect_false(is.null(res$transport))
  # All transported values should be aligned to medoid dimension
  expect_equal(ncol(res$transport[[1]]$subj_values), P_vec[1])
})

# ---- Task 2: Pipeline outputs for downstream consumption ----

test_that("contrasts from pipeline are valid for as.data.frame()", {
  dat <- make_integration_data(S = 5, q = 4, P = 15)
  cvec <- c(1, -1, 0, 0)

  res <- dkge_pipeline(
    betas = dat$betas,
    designs = dat$designs,
    kernel = dat$K,
    contrasts = cvec,
    rank = 2
  )

  # Extract contrasts and convert to data frame
  df <- as.data.frame(res$contrasts)

  expect_s3_class(df, "data.frame")
  expect_true(all(c("contrast", "subject", "component", "value", "method") %in% names(df)))
  expect_gt(nrow(df), 0)

  # Values should be numeric and finite
  expect_true(all(is.finite(df$value)))
})

test_that("inference results from pipeline are valid", {
  dat <- make_integration_data(S = 6, q = 5, P = 20)
  centroids <- make_integration_centroids(S = 6, P = 20)
  cvec <- c(1, -1, 0, 0, 0)

  fit <- dkge_fit(dat$betas, dat$designs, K = dat$K, rank = 3)

  res <- suppressWarnings(dkge_pipeline(
    fit = fit,
    contrasts = cvec,
    transport = list(centroids = centroids, medoid = 1L),
    inference = list(B = 100)
  ))

  # Verify inference has expected structure
  expect_false(is.null(res$inference))
  infer <- res$inference[[1]]

  # p-values in [0, 1]
  expect_true(all(infer$p >= 0 & infer$p <= 1))

  # statistics are finite
  expect_true(all(is.finite(infer$stat)))

  # Result has expected list elements
  expect_true(all(c("stat", "p") %in% names(infer)))
})

test_that("dkge_pipeline with classification produces valid results", {
  dat <- make_integration_data(S = 8, q = 4, P = 15)
  cvec <- c(1, -1, 0, 0)

  fit <- dkge_fit(dat$betas, dat$designs, K = dat$K, rank = 2)

  # Create a simple classification target (2-class: rows are classes, cols are effects)
  target_matrix <- matrix(c(1, 0, 0, 0,
                            0, 1, 0, 0), nrow = 2, byrow = TRUE)
  rownames(target_matrix) <- c("classA", "classB")

  res <- dkge_pipeline(
    fit = fit,
    contrasts = cvec,
    classification = target_matrix
  )

  # Classification result should be populated
  expect_false(is.null(res$classification))
  expect_s3_class(res$classification, "dkge_classification")

  # Should have results
  expect_false(is.null(res$classification$results))
})

test_that("dkge_pipeline with multiple contrasts produces correct structure", {
  dat <- make_integration_data(S = 6, q = 5, P = 20)

  fit <- dkge_fit(dat$betas, dat$designs, K = dat$K, rank = 3)

  # Multiple contrasts as named list
  contrasts_list <- list(
    main1 = c(1, -1, 0, 0, 0),
    main2 = c(0, 0, 1, -1, 0),
    interaction = c(1, -1, -1, 1, 0)
  )

  res <- dkge_pipeline(
    fit = fit,
    contrasts = contrasts_list
  )

  # Should have 3 contrasts
  expect_equal(length(res$contrasts$contrasts), 3)
  expect_equal(names(res$contrasts$contrasts), c("main1", "main2", "interaction"))

  # Each contrast should have values for all subjects
  for (nm in names(res$contrasts$values)) {
    expect_equal(length(res$contrasts$values[[nm]]), 6)
  }
})

test_that("dkge_pipeline diagnostics contain expected metadata", {
  dat <- make_integration_data(S = 5, q = 4, P = 15)
  cvec <- c(1, -1, 0, 0)

  res <- dkge_pipeline(
    betas = dat$betas,
    designs = dat$designs,
    kernel = dat$K,
    contrasts = cvec,
    rank = 2
  )

  diag <- res$diagnostics

  # Diagnostics should include rank and variance explained info
  expect_equal(diag$rank, 2)
  expect_true(!is.null(diag$n_subjects))
  expect_equal(diag$n_subjects, 5)
})
