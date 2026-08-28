# test-pipeline.R
# Basic smoke test for dkge_pipeline

library(testthat)

make_pipeline_inputs <- function(S = 6, q = 3, P = 5, T = 12) {
  set.seed(41)
  betas <- replicate(S, matrix(rnorm(q * P), q, P), simplify = FALSE)
  designs <- replicate(S, {
    X <- matrix(rnorm(T * q), T, q)
    qr.Q(qr(X))
  }, simplify = FALSE)
  list(betas = betas, designs = designs, K = diag(q))
}

test_that("dkge_pipeline returns fit, diagnostics, and contrasts", {
  dat <- make_pipeline_inputs()
  fit <- dkge_fit(dat$betas, dat$designs, K = dat$K, rank = 2)
  cvec <- c(1, -1, 0)

  res <- dkge_pipeline(fit = fit, contrasts = cvec)

  expect_true(all(c("fit", "diagnostics", "contrasts") %in% names(res)))
  expect_s3_class(res$contrasts, "dkge_contrasts")
  expect_equal(res$diagnostics$rank, fit$rank)
})

test_that("dkge_pipeline builds mapper transport but refuses its static inference", {
  data <- create_mismatched_data()
  fit <- dkge_fit(data$betas, data$designs, K = data$K, rank = 2)

  transport_cfg <- list(
    centroids = data$centroids,
    medoid = 1L,
    mapper = dkge_mapper_spec("ridge", lambda = 1e-2),
    betas = data$betas
  )

  expect_error(
    suppressWarnings(dkge_pipeline(
      fit = fit,
      contrasts = c(1, -1, 0),
      transport = transport_cfg,
      inference = list(B = 100)
    )),
    "cannot compose.*legacy descriptive transport",
    class = "dkge_alignment_ineligible_error"
  )

  res <- suppressWarnings(dkge_pipeline(
    fit = fit,
    contrasts = c(1, -1, 0),
    transport = transport_cfg,
    inference = NULL
  ))

  expect_false(is.null(res$transport))
  expect_equal(ncol(res$transport[[1]]$subj_values), nrow(data$centroids[[1]]))
  expect_null(res$inference)
  expect_s3_class(attr(res$transport, "aligned_maps"), "dkge_aligned_maps")
})

test_that("dkge_pipeline accepts service objects", {
  dat <- make_pipeline_inputs(S = 5, q = 3, P = 4, T = 10)
  fit <- dkge_fit(dat$betas, dat$designs, K = dat$K, rank = 2)
  cvec <- c(1, -1, 0)

  transport_spec <- dkge_transport_spec(centroids = replicate(5, matrix(rnorm(12), 4, 3), simplify = FALSE),
                                        medoid = 1L)
  inference_spec <- dkge_inference_spec(
    B = 500,
    tail = "two.sided",
    allow_approximate_alignment = TRUE
  )

  expect_error(
    dkge_pipeline(
      fit = fit,
      contrasts = cvec,
      transport = dkge_transport_service(transport_spec),
      inference = dkge_inference_service(inference_spec)
    ),
    class = "dkge_alignment_ineligible_error"
  )
  transported <- dkge_pipeline(
    fit = fit,
    contrasts = cvec,
    transport = dkge_transport_service(transport_spec),
    inference = NULL
  )
  inferred <- dkge_pipeline(
    fit = fit,
    contrasts = cvec,
    transport = NULL,
    inference = dkge_inference_service(inference_spec)
  )

  expect_s3_class(transported$contrasts, "dkge_contrasts")
  expect_false(is.null(transported$transport))
  expect_null(transported$inference)
  expect_s3_class(inferred$inference, "dkge_inference")
  expect_equal(length(inferred$inference$p_adjusted),
               length(inferred$contrasts$values))
  expect_identical(inferred$inference$metadata$estimator$status,
                   "approximate")
  expect_true(inferred$inference$metadata$estimator$approximate_override)
})

test_that("pipeline inference fails closed and uses one joint max-T family", {
  dat <- make_pipeline_inputs(S = 6, q = 3, P = 5, T = 12)
  fit <- dkge_fit(dat$betas, dat$designs, K = dat$K, rank = 2)
  contrasts <- list(first = c(1, -1, 0), second = c(0, 1, -1))

  expect_error(
    dkge_pipeline(
      fit = fit,
      contrasts = contrasts,
      inference = dkge_inference_spec(B = 100)
    ),
    class = "dkge_alignment_ineligible_error"
  )

  set.seed(260827L)
  inferred <- dkge_pipeline(
    fit = fit,
    contrasts = contrasts,
    inference = dkge_inference_spec(
      B = 100,
      allow_approximate_alignment = TRUE
    )
  )
  expect_s3_class(inferred$inference, "dkge_inference")
  expect_identical(
    inferred$inference$metadata$family_scope,
    "all_contrasts_and_support_locations"
  )
  expect_true(inferred$inference$metadata$shared_subject_signs)

  set.seed(260827L)
  oracle <- dkge:::.infer_signflip(
    inferred$contrasts, 100L, "maxT"
  )
  expect_equal(inferred$inference$p_adjusted, oracle$p_adjusted, tolerance = 0)
  expect_identical(
    inferred$inference$metadata$sign_flips,
    oracle$metadata$sign_flips
  )
})
