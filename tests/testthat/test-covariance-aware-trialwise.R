library(testthat)

ar1_covariance <- function(n, rho) {
  outer(seq_len(n), seq_len(n), function(i, j) rho^abs(i - j))
}

corrected_mean_diagonal <- function(subject) {
  corrected <- tcrossprod(subject$beta) -
    subject$noise_trace * subject$effect_noise_cov
  mean(diag(corrected))
}

test_that("GLS trial sufficient statistics match a dense independent oracle", {
  set.seed(8101)
  n_trials <- 12L
  P <- 7L
  X <- cbind(intercept = 1, trend = seq(-1, 1, length.out = n_trials))
  C <- ar1_covariance(n_trials, 0.45)
  W <- solve(C)
  Y <- matrix(rnorm(n_trials * P), n_trials, P)
  colnames(Y) <- paste0("v", seq_len(P))

  info <- crossprod(X, W %*% X)
  score <- crossprod(X, W %*% Y)
  beta <- solve(info, score)
  residual <- Y - X %*% beta
  residual_df <- n_trials - ncol(X)
  sigma2 <- colSums(residual * (W %*% residual)) / residual_df

  by_covariance <- dkge_trial_subject(
    Y, X, id = "cov", trial_covariance = C
  )
  by_precision <- dkge_trial_subject(
    Y, X, id = "prec", trial_precision = W
  )
  by_operator <- dkge_trial_subject(
    Y, X, id = "op", trial_precision = function(z) W %*% z
  )

  for (subject in list(by_covariance, by_precision, by_operator)) {
    expect_equal(subject$beta, beta, tolerance = 1e-11,
                 ignore_attr = TRUE)
    expect_equal(subject$effect_information, info, tolerance = 1e-11,
                 ignore_attr = TRUE)
    expect_equal(subject$effect_score, score, tolerance = 1e-11,
                 ignore_attr = TRUE)
    expect_equal(subject$effect_noise_cov, solve(info), tolerance = 1e-11,
                 ignore_attr = TRUE)
    expect_equal(unname(subject$residual_variance), unname(sigma2),
                 tolerance = 1e-11)
    expect_equal(subject$residual_df, residual_df)
    expect_equal(subject$noise_trace, sum(sigma2), tolerance = 1e-11)
    expect_equal(
      solve(subject$effect_information, subject$effect_score),
      subject$beta,
      tolerance = 1e-11,
      ignore_attr = TRUE
    )
    expect_null(subject$y)
    top_level_matrices <- Filter(is.matrix, unclass(subject))
    expect_false(any(vapply(top_level_matrices, function(x) {
      identical(dim(x), c(P, P))
    }, logical(1))))
  }
  expect_equal(by_covariance$error_model$trial_dependence, "covariance")
  expect_equal(by_precision$error_model$trial_dependence, "precision")
  expect_equal(by_operator$error_model$trial_dependence, "precision_operator")
  expect_true(all(vapply(
    list(by_covariance, by_precision, by_operator),
    function(x) identical(x$error_model$estimator, "gls"), logical(1)
  )))

  fit <- dkge_fit(
    dkge_data(list(by_covariance, by_precision)),
    K = diag(ncol(X)), rank = 1, w_method = "none"
  )
  expected_pool <- by_covariance$effect_information +
    by_precision$effect_information + 1e-10 * diag(ncol(X))
  expect_equal(fit$R, chol(expected_pool), tolerance = 1e-11,
               ignore_attr = TRUE)
  expect_equal(fit$ruler_source, "effect_information")
})

test_that("direct effect covariance and uncertainty are retained and audited", {
  set.seed(8107)
  X <- cbind(a = 1, b = rep(c(-1, 1), each = 4))
  Y <- matrix(rnorm(8 * 3), 8, 3)
  direct_cov <- matrix(c(0.4, 0.08, 0.08, 0.3), 2, 2,
                       dimnames = list(colnames(X), colnames(X)))
  supplied_variance <- c(0.5, 0.75, 1.25)
  subjects <- lapply(seq_len(2), function(s) {
    dkge_trial_subject(
      Y + 0.1 * s, X, id = paste0("d", s),
      effect_noise_cov = direct_cov,
      residual_variance = supplied_variance
    )
  })

  expect_equal(subjects[[1]]$effect_noise_cov, direct_cov,
               tolerance = 1e-12)
  expect_equal(unname(subjects[[1]]$residual_variance), supplied_variance,
               tolerance = 1e-12)
  expect_equal(subjects[[1]]$noise_trace, sum(supplied_variance),
               tolerance = 1e-12)
  expect_equal(subjects[[1]]$error_model$effect_covariance_source, "supplied")
  expect_equal(subjects[[1]]$error_model$residual_scale_source, "supplied")

  fit <- dkge_fit(
    dkge_data(subjects), K = diag(2), rank = 1,
    w_method = "none", effect_scaling = "none", debias = "analytic"
  )
  expected <- Reduce(`+`, lapply(subjects, function(subject) {
    tcrossprod(subject$beta) - subject$noise_trace * direct_cov
  }))
  expect_equal(fit$effect_moment, expected, tolerance = 1e-11,
               ignore_attr = TRUE)
  expect_equal(fit$error_models[["d1"]]$effect_covariance_source, "supplied")
  expect_equal(fit$provenance$error_models[["d2"]]$estimator, "ols")
})

test_that("IID analytic correction is unbiased in a null Monte Carlo", {
  set.seed(8111)
  n_trials <- 20L
  P <- 2L
  X <- model.matrix(~ 0 + factor(rep(1:2, each = n_trials / 2)))
  colnames(X) <- c("e1", "e2")
  corrected <- replicate(400, {
    Y <- matrix(rnorm(n_trials * P), n_trials, P)
    corrected_mean_diagonal(dkge_trial_subject(Y, X))
  })

  expect_lt(abs(mean(corrected)), 0.06)
})

test_that("AR1 GLS removes the audit bias left by IID correction", {
  set.seed(8119)
  n_trials <- 20L
  P <- 2L
  rho <- 0.6
  X <- model.matrix(~ 0 + factor(rep(1:2, each = n_trials / 2)))
  colnames(X) <- c("e1", "e2")
  C <- ar1_covariance(n_trials, rho)
  L <- t(chol(C))

  corrected <- replicate(500, {
    Y <- L %*% matrix(rnorm(n_trials * P), n_trials, P)
    iid <- dkge_trial_subject(Y, X)
    gls <- dkge_trial_subject(Y, X, trial_covariance = C)
    c(iid = corrected_mean_diagonal(iid),
      gls = corrected_mean_diagonal(gls))
  })
  average <- rowMeans(corrected)

  expect_equal(unname(average[["iid"]]), 0.495, tolerance = 0.12)
  expect_gt(average[["iid"]], 0.35)
  expect_lt(abs(average[["gls"]]), 0.06)
  expect_lt(abs(average[["gls"]]), abs(average[["iid"]]) / 5)
})

test_that("covariance-aware trial inputs fail closed on invalid uncertainty", {
  Y <- matrix(rnorm(18), 6, 3)
  X <- cbind(a = 1, b = rep(c(-1, 1), each = 3))
  I6 <- diag(6)

  expect_error(
    dkge_trial_subject(
      Y, X, trial_covariance = I6, trial_precision = I6
    ),
    "at most one"
  )
  expect_error(
    dkge_trial_subject(Y, X, trial_covariance = diag(5)),
    "T x T"
  )
  asymmetric <- I6
  asymmetric[1, 2] <- 0.5
  expect_error(
    dkge_trial_subject(Y, X, trial_covariance = asymmetric),
    "symmetric"
  )
  singular <- I6
  singular[6, 6] <- 0
  expect_error(
    dkge_trial_subject(Y, X, trial_precision = singular),
    "positive definite"
  )
  expect_error(
    dkge_trial_subject(Y, X, effect_noise_cov = diag(3)),
    "q x q"
  )
  expect_error(
    dkge_trial_subject(Y, X, residual_variance = c(1, 2)),
    "one entry per design effect"
  )
  expect_error(
    dkge_trial_subject(
      Y, X, trial_precision = function(z) z[-1, , drop = FALSE]
    ),
    "T x T"
  )
  expect_warning(
    exploratory <- dkge_trial_subject(
      Y, X, trial_covariance = I6, split = "alternate"
    ),
    "exploratory"
  )
  expect_false(exploratory$split_provenance$independent)
  expect_error(
    dkge_trial_subject(
      Y, X, trial_covariance = ar1_covariance(6, 0.4),
      split = "alternate", split_independent = TRUE
    ),
    "non-zero covariance across halves"
  )

  saturated_X <- diag(3)
  colnames(saturated_X) <- paste0("e", 1:3)
  expect_error(
    dkge_trial_subject(
      matrix(rnorm(9), 3, 3), saturated_X,
      trial_covariance = diag(3)
    ),
    "positive residual degrees of freedom"
  )
  saturated_subjects <- lapply(seq_len(2), function(s) {
    dkge_trial_subject(
      matrix(rnorm(9), 3, 3), saturated_X, id = paste0("sat", s),
      effect_noise_cov = diag(3), residual_variance = rep(1, 3)
    )
  })
  expect_error(
    dkge_fit(
      dkge_data(saturated_subjects), K = diag(3), rank = 1,
      effect_scaling = "none", debias = "analytic"
    ),
    "positive residual degrees of freedom"
  )
})

test_that("local full-rank GLS effects union-align across empty global cells", {
  set.seed(8123)
  make_local <- function(effects, id) {
    X <- model.matrix(~ 0 + factor(rep(seq_along(effects), each = 4)))
    colnames(X) <- effects
    C <- ar1_covariance(nrow(X), 0.3)
    Y <- t(chol(C)) %*% matrix(rnorm(nrow(X) * 3), nrow(X), 3)
    dkge_trial_subject(Y, X, id = id, trial_covariance = C)
  }
  subjects <- list(
    make_local(c("e1", "e2"), "u1"),
    make_local(c("e2", "e3"), "u2")
  )
  data <- suppressWarnings(dkge_data(subjects))

  expect_equal(data$effects, c("e1", "e2", "e3"))
  expect_equal(data$subjects[[1]]$observed_rows, c(1L, 2L))
  expect_equal(data$subjects[[2]]$observed_rows, c(2L, 3L))
  expect_true(all(data$subjects[[1]]$effect_information[3, ] == 0))
  expect_true(all(data$subjects[[1]]$effect_information[, 3] == 0))
  expect_true(all(data$subjects[[2]]$effect_score[1, ] == 0))

  fit <- dkge_fit(
    data, K = diag(3), rank = 1, w_method = "none",
    effect_scaling = "none", debias = "analytic"
  )
  expect_s3_class(fit, "dkge_qspace")
  expect_true(all(vapply(
    fit$provenance$error_models,
    function(model) identical(model$estimator, "gls"), logical(1)
  )))
  expect_true(all(is.finite(fit$Chat)))
})

test_that("stored unweighted trace is recomputed under spatial weights", {
  set.seed(8129)
  X <- cbind(a = 1, b = rep(c(-1, 1), each = 4))
  Y <- matrix(rnorm(8 * 3), 8, 3)
  omega <- c(0.5, 2, 4)
  subject <- dkge_trial_subject(Y, X, omega = omega)

  expect_equal(subject$noise_trace, sum(subject$residual_variance),
               tolerance = 1e-12)
  expect_equal(
    dkge:::.dkge_noise_trace(subject, Omega = omega),
    sum(omega * subject$residual_variance),
    tolerance = 1e-12
  )
})
