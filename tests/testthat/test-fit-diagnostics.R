library(testthat)

make_diagnostic_fit <- function() {
  effects <- c("e1", "e2", "e3")
  make_subject <- function(id, observed) {
    beta <- matrix(seq_len(12) + id, 3, 4,
                   dimnames = list(effects, paste0("v", 1:4)))
    counts <- stats::setNames(c(5, 7, 9), effects)
    counts[-observed] <- 0
    design <- diag(3)
    colnames(design) <- effects
    suppressWarnings(dkge_subject(
      beta, design, id = paste0("s", id),
      observed_rows = observed, effect_n = counts
    ))
  }
  kernel <- diag(3)
  dimnames(kernel) <- list(effects, effects)
  dkge_fit(
    dkge_data(list(make_subject(1, c(1, 2)), make_subject(2, c(2, 3)))),
    K = kernel, rank = 1,
    w_method = "none", effect_scaling = "none",
    effect_weights = dkge_effect_weights("count"),
    missingness = "rescale"
  )
}

test_that("fit summaries expose coverage, ESS, spectrum, and representation", {
  fit <- make_diagnostic_fit()
  diagnosed <- summary(fit)

  expect_s3_class(diagnosed, "summary.dkge")
  expect_equal(diagnosed$n_subjects, 2L)
  expect_equal(diagnosed$n_effects, 3L)
  expect_equal(unname(diagnosed$effect_coverage), c(0.5, 0.5, 1))
  expect_equal(diagnosed$zero_coverage_pairs, 2L)
  expect_true(all(is.finite(diagnosed$pair_ess)))
  expect_equal(diagnosed$effect_weighting, "count")
  expect_equal(diagnosed$representation, "qspace_moment")
  expect_equal(diagnosed$moment_estimator, "observed_second_moment")
  expect_equal(diagnosed$error_models$estimator[["not_recorded"]], 2L)

  output <- capture.output(print(fit))
  output <- paste(output, collapse = "\n")
  expect_match(output, "effect coverage")
  expect_match(output, "pair ESS")
  expect_match(output, "negative spectral mass")
  expect_match(output, "representation: qspace_moment")
  expect_match(output, "trial dependence")
  returned <- NULL
  capture.output(returned <- print(fit))
  expect_identical(returned, fit)
})

test_that("analytic estimator names distinguish IID and covariance-aware fits", {
  set.seed(19501)
  cell <- rep(1:2, each = 8)
  X <- model.matrix(~ 0 + factor(cell, levels = 1:2))
  colnames(X) <- c("e1", "e2")
  truth <- matrix(rnorm(2 * 8), 2, 8)
  make_y <- function() X %*% truth + matrix(rnorm(16 * 8, sd = 0.2), 16, 8)

  iid <- lapply(1:2, function(s) dkge_trial_subject(make_y(), X, id = paste0("i", s)))
  covariance <- diag(16)
  covariance[row(covariance) != col(covariance)] <-
    0.2^abs(row(covariance)[row(covariance) != col(covariance)] -
               col(covariance)[row(covariance) != col(covariance)])
  correlated <- lapply(1:2, function(s) dkge_trial_subject(
    make_y(), X, id = paste0("c", s), trial_covariance = covariance
  ))

  iid_fit <- dkge_fit(dkge_data(iid), K = diag(2), rank = 1,
                      w_method = "none", debias = "analytic")
  covariance_fit <- dkge_fit(dkge_data(correlated), K = diag(2), rank = 1,
                             w_method = "none", debias = "analytic")
  expect_equal(iid_fit$moment_estimator, "analytic_iid_error_correction")
  expect_equal(covariance_fit$moment_estimator,
               "analytic_covariance_aware_correction")
})

test_that("split estimator names expose independence provenance", {
  set.seed(19503)
  cell <- rep(1:2, each = 8)
  X <- model.matrix(~ 0 + factor(cell, levels = 1:2))
  colnames(X) <- c("e1", "e2")
  truth <- matrix(rnorm(2 * 12), 2, 12)
  make_y <- function() X %*% truth + matrix(rnorm(16 * 12, sd = 0.1), 16, 12)
  exploratory <- lapply(1:2, function(s) suppressWarnings(
    dkge_trial_subject(make_y(), X, id = paste0("e", s), split = "within_cell")
  ))
  independent <- lapply(1:2, function(s) dkge_trial_subject(
    make_y(), X, id = paste0("i", s), split = "within_cell",
    split_independent = TRUE
  ))

  exploratory_fit <- suppressWarnings(dkge_fit(
    dkge_data(exploratory), K = diag(2), rank = 1,
    w_method = "none", debias = "split_half"
  ))
  independent_fit <- suppressWarnings(dkge_fit(
    dkge_data(independent), K = diag(2), rank = 1,
    w_method = "none", debias = "split_half"
  ))
  expect_equal(exploratory_fit$moment_estimator,
               "exploratory_split_cross_moment")
  expect_equal(independent_fit$moment_estimator,
               "independent_split_cross_moment")
})
