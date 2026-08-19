library(testthat)

make_reliability_subject <- function(beta, id, counts = NULL,
                                     precision = NULL,
                                     observed_rows = NULL,
                                     residual_variance = NULL,
                                     effect_noise_cov = NULL) {
  effects <- rownames(beta)
  X <- diag(length(effects))
  colnames(X) <- effects
  suppressWarnings(dkge_subject(
    beta, X, id = id, observed_rows = observed_rows,
    effect_n = counts, effect_precision = precision,
    effect_noise_cov = effect_noise_cov,
    residual_variance = residual_variance,
    residual_df = if (is.null(residual_variance)) NULL else 10
  ))
}

test_that("random effects cap a high-count heterogeneous subject", {
  effects <- c("e1", "e2")
  counts <- c(10000, rep(10, 5))
  means <- c(6, -0.5, 0, 0.5, 1, -1)
  subjects <- lapply(seq_along(counts), function(s) {
    B <- rbind(
      e1 = means[[s]] + c(-0.2, 0, 0.1, 0.3),
      e2 = c(-1, 0, 1, 2) + 0.05 * s
    )
    make_reliability_subject(
      B, paste0("s", s),
      counts = stats::setNames(rep(counts[[s]], 2), effects)
    )
  })
  data <- dkge_data(subjects)
  fit_with <- function(spec) {
    dkge_fit(
      data, K = diag(2), rank = 1, w_method = "none",
      effect_scaling = "none", effect_weights = spec
    )
  }
  count_fit <- fit_with(dkge_effect_weights("count"))
  random_fit <- fit_with(dkge_effect_weights(
    "random_effects", within = "count", max_ratio = 10
  ))
  diagnostics <- random_fit$effect_precision_diagnostics

  count_precision <- vapply(count_fit$effect_precision, `[[`, numeric(1), 1)
  random_precision <- vapply(random_fit$effect_precision, `[[`, numeric(1), 1)
  expect_s3_class(random_fit, "dkge_qspace")
  expect_gt(max(count_precision) / sum(count_precision), 0.99)
  expect_lt(max(random_precision) / sum(random_precision), 0.4)
  expect_true(any(diagnostics$capped[, "e2"]))
  expect_gt(diagnostics$tau2[["e1"]], 0)
  expect_equal(diagnostics$within_source, "count")
  expect_equal(diagnostics$tau_method, "DL")
  expect_lt(diagnostics$max_share[["e2"]],
            diagnostics$raw_max_share[["e2"]])
  expect_true(all(is.finite(random_fit$pair_ess)))
})

test_that("noise-based random effects downweight heteroscedastic subjects", {
  effects <- c("e1", "e2")
  noise_scale <- c(0.2, 1, 5)
  subjects <- lapply(seq_along(noise_scale), function(s) {
    B <- rbind(
      e1 = c(-1, 0, 1, 2),
      e2 = c(2, 1, 0, -1)
    ) + 0.01 * s
    Lambda <- diag(c(0.25, 0.5), 2)
    dimnames(Lambda) <- list(effects, effects)
    make_reliability_subject(
      B, paste0("h", s),
      residual_variance = rep(noise_scale[[s]], ncol(B)),
      effect_noise_cov = Lambda
    )
  })
  fit <- dkge_fit(
    dkge_data(subjects), K = diag(2), rank = 1, w_method = "none",
    effect_scaling = "none",
    effect_weights = dkge_effect_weights(
      "random_effects", within = "noise", max_ratio = Inf
    )
  )
  precision <- vapply(fit$effect_precision, `[[`, numeric(1), 1)

  expect_gt(precision[[1]], precision[[2]])
  expect_gt(precision[[2]], precision[[3]])
  expect_equal(fit$effect_precision_diagnostics$within_source, "noise")
  expect_true(all(is.finite(
    fit$effect_precision_diagnostics$within_variance
  )))
})

test_that("random effects handle locally sparse and disjoint effect blocks", {
  make_local <- function(effects, counts, id, offset) {
    B <- matrix(
      seq_len(length(effects) * 3) + offset,
      length(effects), 3,
      dimnames = list(effects, NULL)
    )
    make_reliability_subject(B, id, counts = counts)
  }
  subjects <- list(
    make_local(c("e1", "e2"), c(e1 = 1, e2 = 20), "a1", 0),
    make_local(c("e1", "e2"), c(e1 = 30, e2 = 20), "a2", 1),
    make_local(c("e2", "e3"), c(e2 = 20, e3 = 1), "b1", 2),
    make_local(c("e2", "e3"), c(e2 = 20, e3 = 30), "b2", 3)
  )
  data <- suppressWarnings(dkge_data(subjects))
  fit <- dkge_fit(
    data, K = diag(3), rank = 1, w_method = "none",
    effect_scaling = "none",
    effect_weights = dkge_effect_weights(
      "random_effects", within = "count", max_ratio = 5
    )
  )

  expect_true(all(is.finite(fit$Chat)))
  expect_equal(fit$pair_counts["e1", "e3"], 0)
  expect_equal(fit$pair_weight["e1", "e3"], 0)
  expect_equal(fit$pair_ess["e1", "e3"], 0)
  expect_equal(fit$pair_counts["e1", "e1"], 2)
  expect_equal(fit$pair_counts["e3", "e3"], 2)
  expect_true(all(fit$effect_precision_diagnostics$max_share <= 1))
  expect_equal(
    unname(fit$effect_precision_diagnostics$subjects_per_effect),
    c(2L, 4L, 2L)
  )
})

test_that("random-effects geometry beats count ownership in simulation", {
  set.seed(8423)
  common <- rbind(
    e1 = c(-1, -0.25, 0.5, 1.25),
    e2 = c(1, 0.25, -0.5, -1.25)
  )
  target <- 8 * tcrossprod(common)
  errors <- replicate(50, {
    subjects <- lapply(seq_len(8), function(s) {
      B <- common + matrix(rnorm(8, sd = 0.15), 2, 4)
      n <- 20
      if (s == 1L) {
        B[1, ] <- B[1, ] + 4
        n <- 2000
      }
      make_reliability_subject(
        B, paste0("m", s),
        counts = c(e1 = n, e2 = n)
      )
    })
    data <- dkge_data(subjects)
    count_fit <- dkge_fit(
      data, K = diag(2), rank = 1, w_method = "none",
      effect_scaling = "none",
      effect_weights = dkge_effect_weights("count")
    )
    random_fit <- dkge_fit(
      data, K = diag(2), rank = 1, w_method = "none",
      effect_scaling = "none",
      effect_weights = dkge_effect_weights(
        "random_effects", within = "count", max_ratio = 10
      )
    )
    c(
      count = norm(count_fit$Chat - target, "F"),
      random_effects = norm(random_fit$Chat - target, "F")
    )
  })

  expect_lt(mean(errors["random_effects", ]),
            mean(errors["count", ]) / 2)
  expect_gt(mean(errors["count", ] > errors["random_effects", ]), 0.9)
})

test_that("random-effects validation and automatic source fail closed", {
  expect_error(
    dkge_effect_weights("random_effects", max_ratio = 0.5),
    "at least 1"
  )
  expect_error(
    dkge_effect_weights("random_effects", tau_method = "REML"),
    "must be 'DL'"
  )
  B <- matrix(rnorm(8), 2, 4,
              dimnames = list(c("e1", "e2"), NULL))
  data <- dkge_data(list(
    make_reliability_subject(B, "x1"),
    make_reliability_subject(B, "x2")
  ))
  expect_error(
    dkge_fit(
      data, K = diag(2), rank = 1, effect_scaling = "none",
      effect_weights = dkge_effect_weights("random_effects")
    ),
    "requires complete"
  )
})

test_that("random-effects masks keep the main length check", {
  effects <- c("e1", "e2")
  B <- matrix(rnorm(8), 2, 4, dimnames = list(effects, NULL))
  data <- dkge_data(list(
    make_reliability_subject(B, "s1", counts = c(e1 = 4, e2 = 4)),
    make_reliability_subject(B, "s2", counts = c(e1 = 4, e2 = 4))
  ))
  expect_error(
    dkge:::.dkge_resolve_effect_precision(
      data,
      dkge_effect_weights("random_effects", within = "count"),
      obs_masks = list(TRUE, c(TRUE, TRUE))
    ),
    "one entry per effect"
  )
})
