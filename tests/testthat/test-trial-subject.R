library(testthat)

test_that("trialwise constructor matches QR sufficient statistics", {
  set.seed(73)
  cell <- factor(rep(c("a", "b", "c"), times = c(2, 5, 3)),
                 levels = c("a", "b", "c"))
  X <- model.matrix(~ 0 + cell)
  colnames(X) <- levels(cell)
  truth <- matrix(rnorm(3 * 4), 3, 4)
  Y <- X %*% truth + matrix(rnorm(10 * 4, sd = 0.3), 10, 4)
  colnames(Y) <- paste0("v", 1:4)

  subject <- dkge_trial_subject(Y, X, id = "trial-s1")
  expected_beta <- qr.coef(qr(X), Y)
  expected_resid <- Y - X %*% expected_beta

  expect_s3_class(subject, "dkge_subject")
  expect_equal(subject$beta, expected_beta, tolerance = 1e-12)
  expect_equal(unname(subject$effect_n), c(2, 5, 3))
  expect_equal(subject$effect_noise_cov, solve(crossprod(X)), tolerance = 1e-12)
  expect_equal(unname(subject$residual_variance),
               unname(colSums(expected_resid^2) / 7), tolerance = 1e-12)
  expect_equal(subject$residual_df, 7)
})

test_that("analytic debiasing matches the finite-trial noise oracle", {
  set.seed(88)
  cells <- factor(rep(c("e1", "e2"), each = 5), levels = c("e1", "e2"))
  X <- model.matrix(~ 0 + cells)
  colnames(X) <- levels(cells)
  truth1 <- matrix(c(2, 0.5, 1, -0.5), 2, 2)
  truth2 <- matrix(c(1, -1, 0.5, 1.5), 2, 2)
  Y1 <- X %*% truth1 + matrix(rnorm(20, sd = 0.4), 10, 2)
  Y2 <- X %*% truth2 + matrix(rnorm(20, sd = 0.6), 10, 2)
  subjects <- list(
    dkge_trial_subject(Y1, X, id = "s1"),
    dkge_trial_subject(Y2, X, id = "s2")
  )
  data <- dkge_data(subjects)
  fit <- dkge_fit(data, K = diag(2), rank = 1, w_method = "none",
                  effect_scaling = "none", debias = "analytic")

  expected <- Reduce(`+`, lapply(subjects, function(s) {
    tcrossprod(s$beta) - sum(s$residual_variance) * s$effect_noise_cov
  }))
  expect_equal(fit$effect_moment, expected, tolerance = 1e-11)
  expect_equal(unname(fit$Chat), unname(expected), tolerance = 1e-11)
  expect_equal(fit$debias, "analytic")

  # Analytic debiasing is the only path that stores noise moments, and they must
  # actually carry mass -- an all-zero check here would be vacuous otherwise.
  expect_length(fit$noise_moments, length(subjects))
  expect_true(all(vapply(fit$noise_moments,
                         function(M) any(abs(M) > 0), logical(1))))
  for (i in seq_along(subjects)) {
    s <- subjects[[i]]
    expect_equal(unname(fit$noise_moments[[i]]),
                 unname(sum(s$residual_variance) * s$effect_noise_cov),
                 tolerance = 1e-11)
    expect_equal(unname(fit$effect_moments_raw[[i]]),
                 unname(tcrossprod(s$beta)), tolerance = 1e-11)
    expect_equal(unname(fit$effect_moments[[i]]),
                 unname(fit$effect_moments_raw[[i]] - fit$noise_moments[[i]]),
                 tolerance = 1e-11)
  }

  # A noiseless subject contributes an exactly zero noise moment.
  zero_noise <- dkge:::.dkge_effect_noise_moment(
    list(id = "quiet", effect_noise_cov = subjects[[1]]$effect_noise_cov,
         noise_trace = 0, noise_trace_scope = "weighted"),
    obs_mask = c(TRUE, TRUE)
  )
  expect_true(all(zero_noise == 0))
})

test_that("analytic noise trace respects spatial and voxel diagonals", {
  effects <- c("e1", "e2")
  Lambda <- matrix(c(0.5, 0.1, 0.1, 0.25), 2, 2,
                   dimnames = list(effects, effects))
  subject <- list(
    id = "weighted",
    effect_noise_cov = Lambda,
    residual_variance = c(2, 3),
    noise_trace = NULL
  )
  Omega <- matrix(c(4, 0.7, 0.7, 5), 2, 2)
  expected_trace <- 0.5 * 4 * 2 + 2 * 5 * 3
  actual <- dkge:::.dkge_effect_noise_moment(
    subject, Omega = Omega, voxel_weights = c(0.5, 2),
    obs_mask = c(TRUE, TRUE)
  )
  expect_equal(actual, expected_trace * Lambda, tolerance = 1e-12)
})

test_that("split-half debiasing matches the symmetrized cross-moment oracle", {
  set.seed(91)
  cells <- factor(rep(c("e1", "e2"), each = 6), levels = c("e1", "e2"))
  X <- model.matrix(~ 0 + cells)
  colnames(X) <- levels(cells)
  truth <- matrix(c(1, -0.5, 0.25, 2), 2, 2)
  Y1 <- X %*% truth + matrix(rnorm(24, sd = 0.5), 12, 2)
  Y2 <- X %*% (truth + 0.2) + matrix(rnorm(24, sd = 0.5), 12, 2)
  subjects <- list(
    dkge_trial_subject(Y1, X, id = "s1", split = "within_cell"),
    dkge_trial_subject(Y2, X, id = "s2", split = "within_cell")
  )
  fit <- dkge_fit(dkge_data(subjects), K = diag(2), rank = 1,
                  w_method = "none", effect_scaling = "none",
                  debias = "split_half")

  expected <- Reduce(`+`, lapply(subjects, function(s) {
    cross <- s$split_betas[[1]] %*% t(s$split_betas[[2]])
    (cross + t(cross)) / 2
  }))
  expect_equal(unname(fit$effect_moment), unname(expected), tolerance = 1e-12)
  # `noise_moments` is populated by analytic debiasing only; the old
  # `all(vapply(fit$noise_moments, ...))` was vacuously TRUE over NULL.
  expect_null(fit$noise_moments)
  # There is no separate pre-debias moment on this path: the stored moment is
  # already the split cross-product, so `effect_moments_raw` aliases it.
  expect_identical(fit$effect_moments_raw, fit$effect_moments)
})

test_that("split-half debiasing requires stored full-rank halves", {
  X <- model.matrix(~ 0 + factor(rep(1:2, each = 3)))
  colnames(X) <- c("e1", "e2")
  Y <- matrix(rnorm(12), 6, 2)
  s1 <- dkge_trial_subject(Y, X, id = "s1")
  s2 <- dkge_trial_subject(Y, X, id = "s2")
  expect_error(
    dkge_fit(dkge_data(list(s1, s2)), K = diag(2), rank = 1,
             effect_scaling = "none", debias = "split_half"),
    "requires two stored split beta matrices"
  )

  X_sparse <- model.matrix(~ 0 + factor(c(1, 1, 2)))
  colnames(X_sparse) <- c("e1", "e2")
  expect_error(
    dkge_trial_subject(matrix(rnorm(6), 3, 2), X_sparse,
                       split = "within_cell"),
    "does not leave both half-designs full rank"
  )
})

test_that("analytic debiasing handles an indefinite corrected moment", {
  effects <- c("e1", "e2")
  X <- diag(2)
  colnames(X) <- effects
  B <- matrix(c(3, 0, 0, 0), 2, 2,
              dimnames = list(effects, c("v1", "v2")))
  make_subject <- function(id) {
    Lambda <- diag(c(0.1, 1), 2)
    dimnames(Lambda) <- list(effects, effects)
    suppressWarnings(dkge_subject(
      B, X, id = id,
      effect_noise_cov = Lambda,
      residual_variance = c(v1 = 1, v2 = 1), residual_df = 1
    ))
  }
  fit <- dkge_fit(dkge_data(list(make_subject("s1"), make_subject("s2"))),
                  K = diag(2), rank = 1, w_method = "none",
                  effect_scaling = "none", debias = "analytic")

  expect_equal(fit$rank, 1)
  expect_lt(min(fit$eig_values_full), 0)
  expect_equal(fit$moment_diagnostics$transformed$negative_count, 1)
  expect_gt(fit$moment_diagnostics$transformed$negative_mass, 0)
})

test_that("trialwise validation catches rank and residual-df edge cases", {
  Y <- matrix(rnorm(12), 4, 3)
  X_bad <- cbind(a = 1, b = 1)
  X_bad <- X_bad[rep(1, 4), , drop = FALSE]
  expect_error(dkge_trial_subject(Y, X_bad), "full column rank")

  X_sat <- diag(2)
  colnames(X_sat) <- c("e1", "e2")
  s1 <- dkge_trial_subject(matrix(rnorm(4), 2, 2), X_sat, id = "s1")
  s2 <- dkge_trial_subject(matrix(rnorm(4), 2, 2), X_sat, id = "s2")
  expect_error(
    dkge_fit(dkge_data(list(s1, s2)), K = diag(2), rank = 1,
             effect_scaling = "none", debias = "analytic"),
    "positive residual degrees of freedom"
  )
})

test_that("the 3 by 5 by 4 trialwise design stays in 60-dimensional effect space", {
  set.seed(105)
  grid <- expand.grid(A = 1:3, B = 1:5, response = 1:4)
  counts <- sample(2:5, nrow(grid), replace = TRUE)
  cell_id <- rep(seq_len(nrow(grid)), counts)
  X <- model.matrix(~ 0 + factor(cell_id, levels = seq_len(nrow(grid))))
  colnames(X) <- paste(grid$A, grid$B, grid$response, sep = ":")
  Y <- matrix(rnorm(nrow(X) * 12), nrow(X), 12)
  subject <- suppressWarnings(dkge_trial_subject(Y, X, id = "s60"))

  expect_equal(dim(subject$beta), c(60, 12))
  expect_equal(dim(subject$effect_noise_cov), c(60, 60))
  expect_equal(length(subject$residual_variance), 12)
  expect_equal(unname(subject$effect_n), counts)
  expect_false(any(vapply(subject, function(x) {
    is.matrix(x) && identical(dim(x), c(12L, 12L))
  }, logical(1))))
})
