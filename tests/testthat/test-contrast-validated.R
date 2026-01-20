library(testthat)
library(dkge)

make_validated_fixture <- function(S = 3, q = 3, P = 4, T = 24, seed = 55) {
  set.seed(seed)
  effects <- paste0("eff", seq_len(q))
  betas <- replicate(S, {
    mat <- matrix(rnorm(q * P), q, P)
    rownames(mat) <- effects
    mat
  }, simplify = FALSE)
  designs <- replicate(S, {
    X <- matrix(rnorm(T * q), T, q)
    qr.Q(qr(X))
  }, simplify = FALSE)
  for (i in seq_along(designs)) colnames(designs[[i]]) <- effects
  list(betas = betas, designs = designs)
}

test_that("dkge_contrast_validated returns dual paths and coverage", {
  fx <- make_validated_fixture()
  data_bundle <- dkge_data(fx$betas, designs = fx$designs)
  kernel <- diag(data_bundle$q)
  fit <- dkge_fit(data_bundle, K = kernel, rank = 2, keep_X = TRUE)
  folds <- dkge_define_folds(fit, type = "subject", k = 2)
  contrast <- rep(0, data_bundle$q)
  contrast[1] <- 1
  res <- dkge_contrast_validated(fit,
                                 contrasts = list(main = contrast),
                                 folds = folds,
                                 ridge = 0,
                                 verbose = FALSE)
  expect_s3_class(res, "dkge_contrast_validated")
  expect_true(all(c("observed", "completed", "summary", "provenance") %in% names(res)))
  expect_equal(nrow(res$summary), 1)
  expect_true(all(c("estimate_observed", "estimate_completed", "sensitivity") %in% names(res$summary)))
  expect_equal(length(res$provenance$folds), folds$k)
  expect_true(is.matrix(res$provenance$folds[[1]]$pair_counts_observed))
  expect_true(is.matrix(res$provenance$folds[[1]]$pair_counts_completed))
})

test_that("dkge_contrast_validated handles zero-sum subject weights", {
  fx <- make_validated_fixture(S = 4)
  data_bundle <- dkge_data(fx$betas, designs = fx$designs)
  kernel <- diag(data_bundle$q)
  fit <- dkge_fit(data_bundle, K = kernel, rank = 2, keep_X = TRUE)
  fit$weights <- rep(0, length(fit$Btil))
  folds <- dkge_define_folds(fit, type = "subject", k = 2)
  contrast <- rep(0, data_bundle$q)
  contrast[1] <- 1
  res <- dkge_contrast_validated(fit,
                                 contrasts = list(main = contrast),
                                 folds = folds,
                                 ridge = 0,
                                 verbose = FALSE)
  expect_true(is.finite(res$summary$estimate_observed))
  subject_ids <- fit$subject_ids
  obs_scores <- dkge:::.dkge_validated_subject_means(res$observed$values, subject_ids)
  expected <- mean(obs_scores[, "main"], na.rm = TRUE)
  expect_equal(res$summary$estimate_observed, expected)
})
