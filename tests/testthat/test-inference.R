# test-inference.R
# Simple diagnostics for dkge_infer helpers

library(testthat)
library(dkge)

make_inference_fixture <- function(S = 5, q = 3, P = 4, T = 60, seed = 5151) {
  set.seed(seed)
  effects <- paste0("eff", seq_len(q))
  betas <- replicate(S, {
    mat <- matrix(rnorm(q * P), q, P)
    rownames(mat) <- effects
    mat
  }, simplify = FALSE)
  designs <- replicate(S, {
    X <- matrix(rnorm(T * q), T, q)
    X <- qr.Q(qr(X))
    colnames(X) <- effects
    X
  }, simplify = FALSE)
  fit <- dkge_fit(dkge_data(betas, designs = designs), K = diag(q), rank = 2)
  list(fit = fit, betas = betas, effects = effects)
}

test_that("dkge_infer returns expected structure", {
  fixture <- make_inference_fixture()
  res <- dkge_infer(fixture$fit, c(1, -1, 0))

  expect_s3_class(res, "dkge_inference")
  expect_equal(res$method, "loso")
  expect_equal(res$inference, "signflip")
  expect_equal(res$correction, "maxT")
  expect_equal(length(res$statistics), 1)
  expect_equal(length(res$p_values), 1)
  expect_false(anyNA(res$p_values[[1]]))
  expect_length(res$significant[[1]], ncol(fixture$betas[[1]]))
})

test_that("dkge_infer errors when cluster counts differ without transport", {
  data <- create_mismatched_data()
  fit <- dkge_fit(data$betas, data$designs, K = data$K, rank = 2)

  expect_error(
    dkge_infer(fit, c(1, -1, 0)),
    "Subject cluster counts differ",
    fixed = FALSE
  )
})

test_that("dkge_infer applies mapper-based transport", {
  data <- create_mismatched_data()
  fit <- dkge_fit(data$betas, data$designs, K = data$K, rank = 2)

  transport_cfg <- list(
    centroids = data$centroids,
    medoid = 1L,
    mapper = dkge_mapper_spec("ridge", lambda = 1e-2),
    betas = data$betas
  )

  res <- dkge_infer(fit, c(1, -1, 0), transport = transport_cfg)

  expect_true(is.list(res$transport))
  expect_equal(ncol(res$transport[[1]]$subj_values), nrow(data$centroids[[1]]))
  expect_false(anyNA(res$p_values[[1]]))
})
