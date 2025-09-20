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

test_that("dkge_pipeline uses mapper transport for mismatched clusters", {
  data <- create_mismatched_data()
  fit <- dkge_fit(data$betas, data$designs, K = data$K, rank = 2)

  transport_cfg <- list(
    centroids = data$centroids,
    medoid = 1L,
    mapper = dkge_mapper_spec("ridge", lambda = 1e-2),
    betas = data$betas
  )

  res <- dkge_pipeline(fit = fit,
                       contrasts = c(1, -1, 0),
                       transport = transport_cfg,
                       inference = list(B = 100))

  expect_false(is.null(res$transport))
  expect_equal(ncol(res$transport[[1]]$subj_values), nrow(data$centroids[[1]]))
  expect_false(is.null(res$inference))
  expect_equal(length(res$inference[[1]]$stat), ncol(res$transport[[1]]$subj_values))
})
