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

test_that("dkge_contrast_validated reports estimability metadata", {
  fx <- make_validated_fixture(S = 4)
  data_bundle <- dkge_data(fx$betas, designs = fx$designs)
  fit <- dkge_fit(data_bundle, K = diag(data_bundle$q), rank = 2,
                  keep_X = TRUE)
  fit$kernel_info <- list(term_scope = c(main = "mixed"))
  folds <- dkge_define_folds(fit, type = "subject", k = 2)
  contrast <- rep(0, data_bundle$q)
  contrast[1] <- 1

  expect_warning(
    dkge_contrast_validated(fit,
                            contrasts = list(main = contrast),
                            folds = folds,
                            ridge = 0,
                            verbose = FALSE),
    "between/mixed effects"
  )
  res <- suppressWarnings(
    dkge_contrast_validated(fit,
                            contrasts = list(main = contrast),
                            folds = folds,
                            ridge = 0,
                            verbose = FALSE)
  )

  expect_true(all(c("estimability", "recommended_inference") %in% names(res$summary)))
  expect_equal(res$summary$estimability, "mixed")
  expect_match(res$summary$recommended_inference, "subject-label permutation")
})

test_that("duplicate contrast names yield one summary row each", {
  fx <- make_validated_fixture(S = 4)
  data_bundle <- dkge_data(fx$betas, designs = fx$designs)
  fit <- dkge_fit(data_bundle, K = diag(data_bundle$q), rank = 2, keep_X = TRUE)
  folds <- dkge_define_folds(fit, type = "subject", k = 2)

  c1 <- c(1, 0, 0)
  c2 <- c(0, 0, 1)
  res <- dkge_contrast_validated(fit,
                                 contrasts = list(dup = c1, dup = c2),
                                 folds = folds,
                                 ridge = 0,
                                 verbose = FALSE)

  # merge() on a duplicated `contrast` key would return 4 rows here.
  expect_equal(nrow(res$summary), 2L)
  expect_equal(res$summary$contrast, c("dup", "dup"))
  # The two rows must describe the two different contrasts, not the first twice.
  expect_false(isTRUE(all.equal(res$summary$estimate_observed[1],
                                res$summary$estimate_observed[2])))

  obs <- dkge:::.dkge_validated_subject_means(res$observed$values, fit$subject_ids)
  expect_equal(ncol(obs), 2L)
  expect_false(isTRUE(all.equal(unname(obs[, 1]), unname(obs[, 2]))))
})

test_that("unnamed contrast collections get default names end to end", {
  fx <- make_validated_fixture(S = 4)
  data_bundle <- dkge_data(fx$betas, designs = fx$designs)
  fit <- dkge_fit(data_bundle, K = diag(data_bundle$q), rank = 2, keep_X = TRUE)
  folds <- dkge_define_folds(fit, type = "subject", k = 2)

  c1 <- c(1, 0, 0)
  c2 <- c(0, 0, 1)
  # An unnamed list used to produce a NULL summary (and then an
  # "argument is of length zero" failure when it was column-bound).
  res <- dkge_contrast_validated(fit,
                                 contrasts = list(c1, c2),
                                 folds = folds,
                                 ridge = 0,
                                 verbose = FALSE)
  expect_equal(nrow(res$summary), 2L)
  expect_equal(res$summary$contrast, c("contrast1", "contrast2"))
  expect_true(all(is.finite(res$summary$estimate_observed)))

  # Naming the contrasts must give the same numbers under different labels.
  named <- dkge_contrast_validated(fit,
                                   contrasts = list(a = c1, b = c2),
                                   folds = folds,
                                   ridge = 0,
                                   verbose = FALSE)
  expect_equal(named$summary$contrast, c("a", "b"))
  expect_equal(named$summary$estimate_observed, res$summary$estimate_observed)

  # Partially named lists fill only the blanks.
  partial <- dkge:::.normalize_contrasts(stats::setNames(list(c1, c2), c("", "b")),
                                         fit)
  expect_equal(names(partial), c("contrast1", "b"))

  # `dkge_contrast()` agrees on the same default labels.
  loso <- dkge_contrast(fit, contrasts = list(c1, c2), method = "loso")
  expect_equal(names(loso$values), c("contrast1", "contrast2"))
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
