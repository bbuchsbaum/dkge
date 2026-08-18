library(testthat)

make_resampling_moment_fixture <- function(heldout_shift = 0) {
  effects <- c("e1", "e2", "e3")
  subjects <- lapply(seq_len(5), function(s) {
    B <- matrix(c(
      0.4 * s, 1 + 0.1 * s,
      1.2 - 0.08 * s, -0.2 * s,
      -0.5 + 0.12 * s, 0.7 + 0.05 * s
    ), 3, 2, byrow = TRUE,
    dimnames = list(effects, c("v1", "v2")))
    X <- matrix(c(
      1, 0, 0,
      0, 1, 0,
      0, 0, 1,
      1, 1, 0,
      1, 0, 1,
      0, 1, 1,
      1 + 0.1 * s, 1, 1
    ), 7, 3, byrow = TRUE, dimnames = list(NULL, effects))
    if (s == 5L && heldout_shift != 0) {
      B <- B + heldout_shift * matrix(c(1, -2, 3, 1, -1, 2), 3, 2)
      X <- heldout_shift * X + diag(3)[rep(1:3, length.out = 7), ]
      colnames(X) <- effects
    }
    counts <- stats::setNames(c(4 + s, 15 - s, 3 + 2 * s), effects)
    Lambda <- diag(1 / counts, 3)
    dimnames(Lambda) <- list(effects, effects)
    observed <- if (s == 2L) c(1L, 2L) else seq_along(effects)
    suppressWarnings(dkge_subject(
      B, X, id = paste0("s", s),
      omega = matrix(c(1.2, 0.15, 0.15, 0.9), 2, 2),
      observed_rows = observed,
      effect_n = counts,
      effect_noise_cov = Lambda,
      residual_variance = c(0.12 + 0.01 * s, 0.18 + 0.005 * s),
      residual_df = 8
    ))
  })
  K <- matrix(c(
    1.4, 0.25, 0.1,
    0.25, 1.2, 0.2,
    0.1, 0.2, 1.1
  ), 3, 3, dimnames = list(effects, effects))
  list(subjects = subjects, K = K)
}

fit_resampling_moment_fixture <- function(fixture) {
  dkge_fit(
    dkge_data(fixture$subjects), K = fixture$K, rank = 2,
    w_method = "none", effect_scaling = "pooled_design",
    effect_weights = dkge_effect_weights(
      "random_effects", within = "count", max_ratio = 6
    ),
    debias = "analytic", missingness = "shrink",
    miss_args = list(gamma = 0.5)
  )
}

test_that("fold contexts equal complete direct training-data refits", {
  fixture <- make_resampling_moment_fixture()
  fit <- fit_resampling_moment_fixture(fixture)
  train <- c(1L, 2L, 4L)
  context <- dkge:::.dkge_fold_weight_context(fit, train)

  direct <- fit_resampling_moment_fixture(list(
    subjects = fixture$subjects[train], K = fixture$K
  ))
  expect_equal(context$R, direct$R, tolerance = 1e-12,
               ignore_attr = TRUE,
               info = "the fold ruler must be estimated from training subjects")
  expect_equal(context$Chat, direct$Chat, tolerance = 1e-11,
               ignore_attr = TRUE,
               info = "fold Chat must equal a complete raw-data refit")
  expect_equal(context$pair_counts, direct$pair_counts,
               tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(context$pair_ess, direct$pair_ess,
               tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(context$effect_precision, direct$effect_precision,
               tolerance = 1e-12, ignore_attr = TRUE,
               info = "random-effects precision must be re-estimated in-fold")
  expect_equal(
    context$effect_precision_diagnostics$tau2,
    direct$effect_precision_diagnostics$tau2,
    tolerance = 1e-12
  )
})

test_that("held-out subjects cannot influence any training estimator quantity", {
  clean_fixture <- make_resampling_moment_fixture()
  shifted_fixture <- make_resampling_moment_fixture(heldout_shift = 200)
  clean_fit <- fit_resampling_moment_fixture(clean_fixture)
  shifted_fit <- fit_resampling_moment_fixture(shifted_fixture)
  train <- 1:4

  clean <- dkge:::.dkge_fold_weight_context(clean_fit, train)
  shifted <- dkge:::.dkge_fold_weight_context(shifted_fit, train)
  expect_equal(shifted$R, clean$R, tolerance = 1e-12,
               ignore_attr = TRUE)
  expect_equal(shifted$Chat, clean$Chat, tolerance = 1e-11,
               ignore_attr = TRUE)
  expect_equal(shifted$effect_precision, clean$effect_precision,
               tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(shifted$pair_ess, clean$pair_ess, tolerance = 1e-12,
               ignore_attr = TRUE)
  expect_gt(norm(shifted_fit$R - clean_fit$R, "F"), 1)

  leaked_clean <- clean_fit
  leaked_shifted <- shifted_fit
  leaked_clean$subjects <- NULL
  leaked_shifted$subjects <- NULL
  leaked_clean$debias <- leaked_shifted$debias <- "none"
  leaked_clean$effect_weight_spec <-
    leaked_shifted$effect_weight_spec <- dkge_effect_weights("none")
  leaked_clean$missingness <- leaked_shifted$missingness <- "none"
  leaked_clean_context <- dkge:::.dkge_fold_weight_context(
    leaked_clean, train
  )
  leaked_shifted_context <- dkge:::.dkge_fold_weight_context(
    leaked_shifted, train
  )
  expect_gt(
    norm(leaked_shifted_context$Chat - leaked_clean_context$Chat, "F"),
    1e-4
  )
})

test_that("fold bases and contrast metadata expose the refitted estimator", {
  fixture <- make_resampling_moment_fixture()
  fit <- fit_resampling_moment_fixture(fixture)
  assignments <- list(5L, 1:4)
  folds <- dkge:::.dkge_build_fold_bases(
    fit, assignments, align = FALSE, loader_scope = "heldout"
  )
  direct <- dkge:::.dkge_refit_training_subjects(fit, 1:4)
  expect_equal(
    tcrossprod(folds$folds[[1]]$basis),
    tcrossprod(direct$U),
    tolerance = 1e-9,
    ignore_attr = TRUE,
    info = "fold eigenspace must come from the direct training refit"
  )
  expect_equal(folds$folds[[1]]$R, direct$R, tolerance = 1e-12,
               ignore_attr = TRUE)
  expect_equal(folds$folds[[1]]$estimator$effect_weights, "random_effects")
  expect_equal(folds$folds[[1]]$estimator$debias, "analytic")
  expect_equal(folds$folds[[1]]$pair_ess, direct$pair_ess,
               tolerance = 1e-12, ignore_attr = TRUE)

  contrast <- dkge_contrast(
    fit, c(e1 = 1, e2 = -1, e3 = 0), method = "kfold",
    folds = structure(list(
      type = "custom", k = 2L, assignments = assignments,
      metadata = list(), align = FALSE
    ), class = "dkge_folds"), align = FALSE
  )
  expect_length(contrast$metadata$estimators, 2L)
  expect_length(contrast$metadata$pair_counts, 2L)
  expect_length(contrast$metadata$pair_ess, 2L)
})

test_that("Poisson multiplicities reproduce literal subject duplication", {
  fixture <- make_resampling_moment_fixture()
  fit <- fit_resampling_moment_fixture(fixture)
  multiplicity <- c(2, 0, 1, 2, 0)
  repooled <- dkge:::.dkge_repool_fit(fit, sample_weights = multiplicity)
  expanded <- rep(seq_along(multiplicity), multiplicity)
  literal <- dkge:::.dkge_refit_training_subjects(fit, expanded)

  expect_equal(repooled$refit, "literal_subject_multiset")
  expect_equal(repooled$R, literal$R, tolerance = 1e-12,
               ignore_attr = TRUE)
  expect_equal(repooled$Chat, literal$Chat, tolerance = 1e-11,
               ignore_attr = TRUE)
  expect_equal(repooled$pair_counts, literal$pair_counts,
               tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(repooled$effect_precision, literal$effect_precision,
               tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(repooled$subject_weights, literal$weights,
               tolerance = 1e-12, ignore_attr = TRUE)
  expect_null(repooled$fit)
})

test_that("run-split moments and reliability remain fold-local", {
  set.seed(19011)
  rows <- expand.grid(repeat_id = 1:2, effect = 1:2, run = 1:4)
  X <- model.matrix(~ 0 + factor(rows$effect, levels = 1:2))
  colnames(X) <- c("e1", "e2")
  run <- paste0("run", rows$run)
  subjects <- lapply(1:4, function(s) {
    truth <- rbind(
      c(0.5 * s, -0.2, 0.1, 0.3, -0.4, 0.6),
      c(-0.1, 0.3 * s, 0.5, -0.2, 0.7, -0.3)
    )
    Y <- X %*% truth +
      matrix(rnorm(nrow(X) * 6, sd = 0.05), nrow(X), 6)
    dkge_trial_subject(
      Y, X, id = paste0("split", s), split = "run",
      run_labels = run, effect_precision = "split_half"
    )
  })
  fit <- dkge_fit(
    dkge_data(subjects), K = diag(2), rank = 1, w_method = "none",
    effect_weights = dkge_effect_weights("precision"),
    debias = "split_half"
  )
  context <- dkge:::.dkge_fold_weight_context(fit, 1:3)
  direct <- dkge:::.dkge_refit_training_subjects(fit, 1:3)
  expect_equal(context$Chat, direct$Chat, tolerance = 1e-11,
               ignore_attr = TRUE)
  expect_equal(context$effect_precision, direct$effect_precision,
               tolerance = 1e-12, ignore_attr = TRUE)
  expect_true(all(vapply(direct$subjects, function(subject) {
    isTRUE(subject$split_provenance$run_disjoint)
  }, logical(1))))
})

test_that("analytic shortcuts fall back for covariance-aware moments", {
  fixture <- make_resampling_moment_fixture()
  fit <- fit_resampling_moment_fixture(fixture)
  result <- dkge_analytic_loso(fit, 5, c(1, -1, 0))
  expect_equal(result$method, "fallback")
  expect_equal(result$diagnostic$reason, "pair_normalized_pooling")

  precision_fit <- dkge_fit(
    dkge_data(fixture$subjects), K = fixture$K, rank = 2,
    w_method = "none", debias = "analytic"
  )
  covariance_result <- dkge_analytic_loso(
    precision_fit, 5, c(1, -1, 0)
  )
  expect_equal(covariance_result$method, "fallback")
  expect_equal(covariance_result$diagnostic$reason,
               "covariance_aware_moment")
})

test_that("supported sign-flip p-values and projected intervals calibrate under a null", {
  set.seed(19017)
  n_rep <- 40L
  p_value <- numeric(n_rep)
  covered <- logical(n_rep)
  for (i in seq_len(n_rep)) {
    sample <- rnorm(20)
    p_value[[i]] <- dkge_signflip_maxT(
      matrix(sample, ncol = 1), B = 199
    )$p[[1]]
    interval <- dkge_bootstrap_projected(
      lapply(sample, function(x) x), B = 299, seed = 20000 + i,
      return_samples = FALSE
    )$medoid$ci[1, ]
    covered[[i]] <- interval[[1]] <= 0 && interval[[2]] >= 0
  }
  expect_lte(mean(p_value <= 0.05), 0.15)
  expect_gt(stats::median(p_value), 0.2)
  expect_gte(mean(covered), 0.8)
})
