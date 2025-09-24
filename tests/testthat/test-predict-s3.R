# test-predict-s3.R

library(testthat)

test_that("predict.dkge matches dkge_predict", {
  fixture <- make_small_fit()
  new_betas <- replicate(2, matrix(rnorm(fixture$q * 3), fixture$q, 3), simplify = FALSE)
  contrasts <- diag(fixture$q)

  direct <- dkge_predict(fixture$fit, B_list = new_betas, contrasts = contrasts)
  via_method <- predict(fixture$fit,
                        newdata = list(betas = new_betas, contrasts = contrasts))

  expect_equal(via_method$values, direct$values)
  expect_equal(via_method$A_list, direct$A_list)

  via_args <- predict(fixture$fit, B_list = new_betas, contrasts = contrasts,
                      return_loadings = FALSE)
  expect_false("A_list" %in% names(via_args))

  expect_error(predict(fixture$fit, newdata = list(betas = new_betas)),
               "Provide `contrasts`", fixed = TRUE)
})

test_that("predict.dkge_model dispatches to dkge_predict", {
  fixture <- make_small_fit()
  model <- dkge_freeze(fixture$fit)
  new_betas <- fixture$betas
  contrasts <- diag(fixture$q)[, 1:2, drop = FALSE]

  expected <- dkge_predict(model, B_list = new_betas, contrasts = contrasts)
  via_method <- predict(model, newdata = list(betas = new_betas, contrasts = contrasts))

  expect_equal(via_method$values, expected$values)
  expect_equal(via_method$A_list, expected$A_list)
})
