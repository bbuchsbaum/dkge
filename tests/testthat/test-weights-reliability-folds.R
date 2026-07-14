# test-weights-reliability-folds.R
# Regression: reliability weights must survive fold/LOSO subsetting.

library(testthat)

test_that("reliability weights do not crash under fold subsetting", {
  fit <- toy_fold_fit(nsub = 4L, Q = 6L, V = 8L)
  run2 <- toy_betas(nsub = 4L, Q = 6L, V = 8L, seed = 999)
  fit$weight_spec <- dkge_weights(adapt = "reliability", B_list2 = run2)

  train_ids <- c(1L, 2L, 3L)
  # Previously errored: B_list2 (length 4) vs B_train (length 3) length mismatch.
  expect_no_error(dkge:::.dkge_fold_weight_context(fit, train_ids))
  ctx <- dkge:::.dkge_fold_weight_context(fit, train_ids)

  w <- ctx$weights$total
  expect_length(w, 8L)
  expect_true(all(is.finite(w)))
  expect_true(all(w >= 0))

  # The fold context must use the subset second run, matching a direct call on
  # the training subjects only.
  expected <- dkge:::.dkge_adapt_weights(
    fit$Btil[train_ids],
    adapt = "reliability",
    B_list2 = run2[train_ids],
    winsor = fit$weight_spec$shrink$winsor
  )
  expect_equal(unname(ctx$weights$adapt), unname(expected), tolerance = 1e-8)
})
