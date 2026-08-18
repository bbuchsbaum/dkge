library(testthat)

make_effect_weight_subject <- function(beta, n = NULL, precision = NULL, id) {
  q <- nrow(beta)
  effects <- rownames(beta) %||% paste0("e", seq_len(q))
  rownames(beta) <- effects
  X <- diag(q)
  colnames(X) <- effects
  suppressWarnings(dkge_subject(
    beta, X, id = id, effect_n = n, effect_precision = precision
  ))
}

test_that("count weighting matches a hand-computed pairwise oracle", {
  B1 <- matrix(c(1, 2), 2, 1, dimnames = list(c("e1", "e2"), "v"))
  B2 <- matrix(c(3, 4), 2, 1, dimnames = list(c("e1", "e2"), "v"))
  p1 <- c(e1 = 1, e2 = 4)
  p2 <- c(e1 = 9, e2 = 1)
  data <- dkge_data(list(
    make_effect_weight_subject(B1, n = p1, id = "s1"),
    make_effect_weight_subject(B2, n = p2, id = "s2")
  ))

  fit <- dkge_fit(
    data, K = diag(2), rank = 1, w_method = "none",
    effect_scaling = "none", effect_weights = dkge_effect_weights("count")
  )

  h1 <- tcrossprod(sqrt(p1))
  h2 <- tcrossprod(sqrt(p2))
  expected <- 2 * (h1 * tcrossprod(B1) + h2 * tcrossprod(B2)) / (h1 + h2)
  expected_ess <- (h1 + h2)^2 / (h1^2 + h2^2)

  expect_equal(unname(fit$effect_moment), unname(expected), tolerance = 1e-12)
  expect_equal(unname(fit$Chat), unname(expected), tolerance = 1e-12)
  expect_equal(unname(fit$pair_ess), unname(expected_ess), tolerance = 1e-12)
})

test_that("effect precision has the expected scaling and zero-weight invariants", {
  B1 <- matrix(c(1, 2, 3, 4), 2, 2,
               dimnames = list(c("e1", "e2"), c("v1", "v2")))
  B2 <- matrix(c(4, 1, 2, 3), 2, 2,
               dimnames = dimnames(B1))
  p1 <- c(e1 = 1, e2 = 3)
  p2 <- c(e1 = 2, e2 = 5)

  fit_for <- function(scale = 1) {
    data <- dkge_data(list(
      make_effect_weight_subject(B1, precision = scale * p1, id = "s1"),
      make_effect_weight_subject(B2, precision = scale * p2, id = "s2")
    ))
    dkge_fit(data, K = diag(2), rank = 1, w_method = "none",
             effect_scaling = "none",
             effect_weights = dkge_effect_weights("precision"))
  }
  expect_equal(fit_for(1)$Chat, fit_for(100)$Chat, tolerance = 1e-12)

  data_zero <- dkge_data(list(
    make_effect_weight_subject(B1, precision = c(e1 = 1, e2 = 0), id = "s1"),
    make_effect_weight_subject(B2, precision = c(e1 = 1, e2 = 1), id = "s2")
  ))
  fit_zero <- dkge_fit(data_zero, K = diag(2), rank = 1,
                       w_method = "none", effect_scaling = "none",
                       effect_weights = dkge_effect_weights("precision"))
  expect_equal(fit_zero$effect_moment["e2", "e2"],
               2 * tcrossprod(B2)[2, 2], tolerance = 1e-12)
  expect_equal(fit_zero$pair_ess["e2", "e2"], 1, tolerance = 1e-12)
})

test_that("unit effect precision reproduces the legacy full-coverage sum", {
  set.seed(41)
  effects <- paste0("e", 1:3)
  subjects <- lapply(1:3, function(s) {
    B <- matrix(rnorm(12), 3, 4, dimnames = list(effects, NULL))
    make_effect_weight_subject(B, n = stats::setNames(rep(1, 3), effects),
                               id = paste0("s", s))
  })
  data <- dkge_data(subjects)
  fit_none <- dkge_fit(data, K = diag(3), rank = 2, w_method = "none",
                       effect_scaling = "none")
  fit_count <- dkge_fit(data, K = diag(3), rank = 2, w_method = "none",
                        effect_scaling = "none",
                        effect_weights = dkge_effect_weights("count"))
  expect_equal(fit_count$Chat, fit_none$Chat, tolerance = 1e-12)
})

test_that("effective pair count controls weighted masking", {
  B <- matrix(c(2, 1, 1, 2), 2, 2,
              dimnames = list(c("e1", "e2"), NULL))
  data <- dkge_data(list(
    make_effect_weight_subject(B, precision = c(e1 = 100, e2 = 1), id = "s1"),
    make_effect_weight_subject(B, precision = c(e1 = 1, e2 = 100), id = "s2")
  ))
  fit <- dkge_fit(data, K = diag(2), rank = 1, w_method = "none",
                  effect_scaling = "none",
                  effect_weights = dkge_effect_weights("precision"),
                  missingness = "mask", miss_args = list(min_pairs = 1.5))
  expect_lt(fit$pair_ess["e1", "e1"], 1.1)
  expect_lt(fit$pair_ess["e2", "e2"], 1.1)
  expect_equal(unname(diag(fit$effect_moment)), c(0, 0))
  expect_gt(fit$pair_ess["e1", "e2"], 1.9)
})

test_that("count weighting requires count sufficient statistics", {
  B <- matrix(rnorm(6), 2, 3, dimnames = list(c("e1", "e2"), NULL))
  X <- diag(2)
  colnames(X) <- rownames(B)
  data <- suppressWarnings(dkge_data(list(
    dkge_subject(B, X, id = "s1"),
    dkge_subject(B, X, id = "s2")
  )))
  expect_error(
    dkge_fit(data, K = diag(2), effect_scaling = "none",
             effect_weights = dkge_effect_weights("count")),
    "has no `effect_n`"
  )
})

test_that("fold pooling inherits effect reliability and missingness policy", {
  effects <- c("e1", "e2")
  subjects <- lapply(1:3, function(s) {
    B <- matrix(c(s, s + 1, 2 * s, 1), 2, 2,
                dimnames = list(effects, NULL))
    make_effect_weight_subject(
      B,
      precision = stats::setNames(c(s, 4 - s), effects),
      id = paste0("s", s)
    )
  })
  fit <- dkge_fit(
    dkge_data(subjects), K = diag(2), rank = 1, w_method = "none",
    effect_scaling = "none",
    effect_weights = dkge_effect_weights("precision"),
    missingness = "shrink", miss_args = list(gamma = 0.5)
  )

  ctx <- dkge:::.dkge_fold_weight_context(fit, train_ids = c(1, 3))
  oracle <- dkge:::.dkge_repool_fit(fit, indices = c(1, 3))
  expect_equal(ctx$Chat, oracle$Chat, tolerance = 1e-12)
  expect_equal(ctx$pair_ess, oracle$pair_ess, tolerance = 1e-12)
  expect_equal(ctx$missingness, "shrink")

  analytic <- dkge_analytic_loso(fit, 2, c(1, -1))
  expect_equal(analytic$method, "fallback")
  expect_equal(analytic$diagnostic$reason, "pair_normalized_pooling")
})
