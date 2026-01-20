# test-component-write.R
# Coverage for dkge_write_component_stats export helper

library(testthat)

set.seed(777)

test_that("dkge_write_component_stats writes tidy CSV output", {
  fx <- make_small_fit(S = 3, q = 3, P = 4, T = 16, rank = 2, seed = 29)
  fit <- fx$fit
  betas <- fx$betas
  centroids <- replicate(length(betas), matrix(runif(ncol(betas[[1]]) * 3), ncol = 3), simplify = FALSE)
  fit$centroids <- centroids

  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp), add = TRUE)

  res <- dkge_write_component_stats(fit,
                                    file = tmp,
                                    mapper = dkge_mapper_spec("ridge", lambda = 1e-3),
                                    centroids = centroids,
                                    components = 1,
                                    inference = list(type = "parametric"))

  expect_true(file.exists(tmp))
  expect_gt(file.size(tmp), 0)
  expect_true(is.list(res))
  expect_true("summary" %in% names(res))

  csv_data <- utils::read.csv(tmp, stringsAsFactors = FALSE)
  expect_true(all(c("component", "cluster", "stat", "p", "p_adj", "significant") %in% names(csv_data)))
  expect_equal(nrow(csv_data), nrow(res$summary))
})
