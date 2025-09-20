# test-inference-transport.R
# Ensure inference stack respects transport metadata

library(testthat)

set.seed(2024)

test_that("dkge_infer errors without transport when cluster sizes differ", {
  data <- create_mismatched_data()
  fit <- dkge_fit(data$betas, data$designs, K = data$K, rank = 2)
  expect_error(
    dkge_infer(fit, c(1, -1, 0), method = "loso", n_perm = 100),
    "Subject cluster counts differ",
    fixed = FALSE
  )
})

test_that("dkge_infer applies mapper-based transports", {
  data <- create_mismatched_data()
  fit <- dkge_fit(data$betas, data$designs, K = data$K, rank = 2)

  transport_cfg <- list(
    centroids = data$centroids,
    medoid = 1L,
    mapper = dkge_mapper_spec("ridge", lambda = 1e-2),
    betas = data$betas
  )

  res <- dkge_infer(fit, c(1, -1, 0), method = "loso",
                    n_perm = 100,
                    transport = transport_cfg)

  expect_s3_class(res, "dkge_inference")
  expect_true(is.list(res$transport))
  expect_equal(ncol(res$transport[[1]]$subj_values), nrow(data$centroids[[1]]))
})
