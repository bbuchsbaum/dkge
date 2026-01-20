# test-predict-stream.R
# Coverage for streaming prediction helper

library(testthat)

set.seed(2024)

test_that("dkge_predict_stream maps loader outputs to contrast values", {
  fx <- make_small_fit(S = 2, q = 3, P = 4, T = 15, rank = 2, seed = 17)
  fit <- fx$fit
  betas <- fx$betas

  loader <- list(
    n = function() length(betas),
    B = function(s) betas[[s]]
  )

  q <- nrow(fit$U)
  contrasts <- diag(q)[, 1:2, drop = FALSE]
  colnames(contrasts) <- c("c1", "c2")

  res <- dkge_predict_stream(fit, loader, contrasts)

  expect_length(res$values, loader$n())
  expect_equal(names(res$values), paste0("subj", seq_len(loader$n())))
  expect_length(res$A_list, loader$n())

  first_values <- res$values[[1]]
  expect_true(is.matrix(first_values))
  expect_equal(ncol(first_values), ncol(contrasts))
  expect_equal(colnames(first_values), colnames(contrasts))
  expect_equal(nrow(first_values), ncol(betas[[1]]))

  first_A <- res$A_list[[1]]
  expect_true(is.matrix(first_A))
  expect_equal(dim(first_A), c(ncol(betas[[1]]), ncol(fit$U)))

  all_names_ok <- vapply(res$values, function(mat) identical(colnames(mat), colnames(contrasts)), logical(1))
  expect_true(all(all_names_ok))
})
