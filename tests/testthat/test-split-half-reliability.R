library(testthat)

make_run_split_fixture <- function(seed = 8301) {
  set.seed(seed)
  rows <- expand.grid(
    repeat_id = seq_len(2), effect = seq_len(2), run = seq_len(4)
  )
  effects <- c("e1", "e2")
  X <- model.matrix(~ 0 + factor(rows$effect, levels = seq_len(2)))
  colnames(X) <- effects
  Y <- matrix(rnorm(nrow(X) * 5), nrow(X), 5)
  list(X = X, Y = Y, run = paste0("run", rows$run), effects = effects)
}

test_that("run and explicit splits record auditable independence provenance", {
  fixture <- make_run_split_fixture()
  subject <- dkge_trial_subject(
    fixture$Y, fixture$X, id = "runs",
    split = "run", run_labels = fixture$run
  )

  first <- subject$split_provenance$half_indices$first
  second <- subject$split_provenance$half_indices$second
  expect_true(subject$split_provenance$independent)
  expect_true(subject$split_provenance$run_disjoint)
  expect_equal(subject$split_provenance$independence_basis, "run_disjoint")
  expect_length(intersect(unique(fixture$run[first]),
                          unique(fixture$run[second])), 0)
  expect_equal(
    subject$split_betas[[1]],
    qr.coef(qr(fixture$X[first, , drop = FALSE]),
            fixture$Y[first, , drop = FALSE]),
    tolerance = 1e-12,
    ignore_attr = TRUE
  )
  expect_equal(
    subject$split_betas[[2]],
    qr.coef(qr(fixture$X[second, , drop = FALSE]),
            fixture$Y[second, , drop = FALSE]),
    tolerance = 1e-12,
    ignore_attr = TRUE
  )
  expect_equal(unname(subject$split_provenance$half_ranks), c(2L, 2L))
  expect_true(all(subject$split_provenance$effect_counts > 0))

  explicit_labels <- ifelse(fixture$run %in% c("run1", "run2"), "A", "B")
  explicit <- dkge_trial_subject(
    fixture$Y, fixture$X, split_labels = explicit_labels,
    run_labels = fixture$run
  )
  expect_equal(explicit$split_provenance$mode, "explicit")
  expect_true(explicit$split_provenance$independent)
  expect_equal(explicit$split_provenance$independence_basis, "run_disjoint")
})

test_that("same-source automatic splits are exploratory unless justified", {
  fixture <- make_run_split_fixture()
  expect_warning(
    exploratory <- dkge_trial_subject(
      fixture$Y, fixture$X, split = "within_cell"
    ),
    "exploratory"
  )
  expect_false(exploratory$split_provenance$independent)
  expect_equal(
    exploratory$split_provenance$independence_basis,
    "exploratory_same_source"
  )

  expect_silent(
    declared <- dkge_trial_subject(
      fixture$Y, fixture$X, split = "within_cell",
      split_independent = TRUE
    )
  )
  expect_true(declared$split_provenance$independent)
  expect_equal(declared$split_provenance$independence_basis,
               "user_declared")

  again <- suppressWarnings(dkge_trial_subject(
    fixture$Y, fixture$X, split = "within_cell"
  ))
  expect_equal(
    exploratory$split_provenance$split_labels,
    again$split_provenance$split_labels
  )
  expect_equal(exploratory$split_betas, again$split_betas,
               tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(unname(exploratory$split_provenance$half_ranks), c(2L, 2L))
  expect_true(all(exploratory$split_provenance$effect_counts > 0))
})

test_that("independent run halves recover the signal cross-moment", {
  set.seed(8311)
  q <- 2L
  P <- 4L
  per_cell_run <- 4L
  cell <- rep(rep(seq_len(q), each = per_cell_run), times = 2)
  run <- rep(c("r1", "r2"), each = q * per_cell_run)
  X <- model.matrix(~ 0 + factor(cell, levels = seq_len(q)))
  colnames(X) <- c("e1", "e2")
  truth <- matrix(c(
    1.0, -0.5, 0.25, 1.5,
    -0.25, 0.75, 1.25, -1.0
  ), q, P, byrow = TRUE)
  target <- tcrossprod(truth)

  average <- Reduce(`+`, lapply(seq_len(300), function(iteration) {
    Y <- X %*% truth + matrix(rnorm(nrow(X) * P, sd = 0.8), nrow(X), P)
    subject <- dkge_trial_subject(
      Y, X, split = "run", run_labels = run
    )
    cross <- subject$split_betas[[1]] %*% t(subject$split_betas[[2]])
    (cross + t(cross)) / 2
  })) / 300

  expect_lt(norm(average - target, "F") / norm(target, "F"), 0.08)
})

test_that("correlated same-run alternation retains the expected cross bias", {
  set.seed(8317)
  n_trials <- 12L
  P <- 2L
  X <- matrix(1, n_trials, 1, dimnames = list(NULL, "mean"))
  covariance <- outer(
    seq_len(n_trials), seq_len(n_trials),
    function(i, j) 0.8^abs(i - j)
  )
  L <- t(chol(covariance))

  biased <- replicate(350, {
    Y <- L %*% matrix(rnorm(n_trials * P), n_trials, P)
    subject <- suppressWarnings(dkge_trial_subject(
      Y, X, split = "within_cell"
    ))
    as.numeric(subject$split_betas[[1]] %*% t(subject$split_betas[[2]]))
  })
  independent <- replicate(350, {
    Y <- matrix(rnorm(n_trials * P), n_trials, P)
    subject <- suppressWarnings(dkge_trial_subject(
      Y, X, split = "within_cell"
    ))
    as.numeric(subject$split_betas[[1]] %*% t(subject$split_betas[[2]]))
  })

  expect_gt(mean(biased), 0.2)
  expect_lt(abs(mean(independent)), 0.08)
  expect_gt(mean(biased), abs(mean(independent)) + 0.15)
})

test_that("split reliability exports bounded finite low-count precision", {
  effects <- c("stable", "constant", "sparse")
  pattern <- c(-2, -1, 0, 1, 2, 3)
  spatial_truth <- rbind(
    stable = pattern,
    constant = rep(2, length(pattern)),
    sparse = c(3, 1, -1, 2, -2, 0)
  )
  per_run_effect <- c(rep(1L, 4), rep(2L, 4), 3L)
  cell <- rep(per_run_effect, times = 2)
  run <- rep(c("r1", "r2"), each = length(per_run_effect))
  X <- model.matrix(~ 0 + factor(cell, levels = seq_along(effects)))
  colnames(X) <- effects
  Y <- X %*% spatial_truth

  subject <- dkge_trial_subject(
    Y, X, split = "run", run_labels = run,
    effect_precision = "split_half"
  )
  precision <- subject$effect_precision
  diagnostics <- subject$split_reliability

  expect_true(all(is.finite(precision)))
  expect_true(all(precision >= 0 & precision <= 1))
  expect_equal(unname(precision[["constant"]]), 0)
  expect_gt(precision[["stable"]], precision[["sparse"]])
  expect_equal(unname(diagnostics$min_half_count), c(4, 4, 1))
  expect_true(diagnostics$independent)

  no_count_shrink <- dkge_split_effect_precision(subject, count_prior = 0)
  expect_equal(unname(no_count_shrink[c("stable", "sparse")]), c(1, 1),
               tolerance = 1e-12)
  expect_equal(unname(no_count_shrink[["constant"]]), 0)

  too_few_features <- subject
  too_few_features$split_betas <- lapply(
    too_few_features$split_betas,
    function(B) B[, 1:2, drop = FALSE]
  )
  low_feature_precision <- dkge_split_effect_precision(too_few_features)
  expect_true(all(is.finite(low_feature_precision)))
  expect_true(all(low_feature_precision == 0))

  subject2 <- dkge_trial_subject(
    Y + 0.01 * X %*% spatial_truth, X,
    split = "run", run_labels = run,
    effect_precision = "split_half"
  )
  fit <- dkge_fit(
    dkge_data(list(subject, subject2)), K = diag(3), rank = 1,
    w_method = "none", effect_scaling = "none",
    effect_weights = dkge_effect_weights("precision"),
    debias = "split_half"
  )
  expect_s3_class(fit, "dkge_qspace")
  expect_true(all(is.finite(fit$Chat)))
})

test_that("split label validation fails closed", {
  fixture <- make_run_split_fixture()
  expect_error(
    dkge_trial_subject(
      fixture$Y, fixture$X, split = "explicit",
      split_labels = rep("one", nrow(fixture$X))
    ),
    "exactly two"
  )
  expect_error(
    dkge_trial_subject(
      fixture$Y, fixture$X, split = "run",
      run_labels = rep("one", nrow(fixture$X))
    ),
    "at least two runs"
  )
  expect_error(
    dkge_trial_subject(
      fixture$Y, fixture$X, split = "run",
      run_labels = fixture$run[-1]
    ),
    "one entry per trial"
  )
  expect_error(
    dkge_trial_subject(
      fixture$Y, fixture$X, effect_precision = "split_half"
    ),
    "requires a split-half"
  )
})
