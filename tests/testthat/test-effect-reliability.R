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

test_that("missingness='rescale' is not a no-op under effect weighting", {
  set.seed(9)
  effects <- paste0("e", 1:3)
  subjects <- lapply(1:3, function(s) {
    B <- matrix(rnorm(9), 3, 3, dimnames = list(effects, NULL))
    make_effect_weight_subject(
      B, precision = stats::setNames(c(s, 4 - s, 2), effects),
      id = paste0("s", s)
    )
  })
  data <- dkge_data(subjects)
  args <- list(data, K = diag(3), rank = 2, w_method = "none",
               effect_scaling = "none",
               effect_weights = dkge_effect_weights("precision"))

  fit_none <- do.call(dkge_fit, c(args, list(missingness = "none")))
  fit_resc <- do.call(dkge_fit, c(args, list(missingness = "rescale")))

  expect_false(isTRUE(all.equal(fit_none$Chat, fit_resc$Chat)))

  # rescale must undo exactly the cohort-mass factor that scaled the pooled
  # per-pair mean up, i.e. land on the precision-weighted per-pair mean
  # sum_s h_s M_s / sum_s h_s. Dividing by `pair_ess` instead (the old code)
  # normalised twice and inflated every entry by cohort_mass / pair_count.
  h <- lapply(data$effect_precision, function(p) tcrossprod(sqrt(p)))
  M <- lapply(data$betas, tcrossprod)
  per_pair <- Reduce(`+`, Map(`*`, h, M)) / Reduce(`+`, h)
  expect_equal(unname(fit_resc$effect_moment), unname(per_pair),
               tolerance = 1e-12)

  cohort_mass <- sum(fit_none$weights)
  expect_equal(unname(fit_resc$effect_moment),
               unname(fit_none$effect_moment / cohort_mass),
               tolerance = 1e-12)
})

test_that("count-weighted rescale equals the unweighted per-pair mean under unequal coverage", {
  effects <- c("e1", "e2", "e3")
  design <- diag(3)
  colnames(design) <- effects
  betas <- list(
    matrix(c(1, 2, 3, 2, 0, 1), 3, 2, dimnames = list(effects, NULL)),
    matrix(c(0, 1, 4, 3, 2, 2), 3, 2, dimnames = list(effects, NULL)),
    matrix(c(2, 1, 0, 1, 3, 5), 3, 2, dimnames = list(effects, NULL))
  )
  # Unequal coverage: e3 is seen by one subject only, e1 by two.
  rows <- list(c(1L, 2L), c(2L, 3L), c(1L, 2L))
  make_one <- function(i) {
    B <- betas[[i]]
    B[-rows[[i]], ] <- 0
    suppressWarnings(dkge_subject(
      B, design, id = paste0("s", i), observed_rows = rows[[i]],
      effect_n = stats::setNames(rep(1, 3), effects)
    ))
  }
  data <- dkge_data(lapply(1:3, make_one))
  args <- list(data, K = diag(3), rank = 2, w_method = "none",
               effect_scaling = "none", missingness = "rescale")

  fit_plain <- do.call(dkge_fit, args)
  fit_count <- do.call(dkge_fit, c(args,
                                   list(effect_weights = dkge_effect_weights("count"))))

  # Hand-computed per-pair mean: sum over subjects observing the pair, divided
  # by how many observed it.
  masks <- lapply(data$provenance$obs_mask, as.numeric)
  structural <- lapply(masks, tcrossprod)
  num <- Reduce(`+`, Map(`*`, structural, lapply(data$betas, tcrossprod)))
  den <- Reduce(`+`, structural)
  oracle <- matrix(0, 3, 3)
  oracle[den > 0] <- num[den > 0] / den[den > 0]

  expect_equal(unname(den), unname(fit_plain$pair_counts), tolerance = 1e-12)
  expect_equal(unname(fit_plain$effect_moment), unname(oracle), tolerance = 1e-12)
  expect_equal(unname(fit_count$effect_moment), unname(oracle), tolerance = 1e-12)
  expect_equal(fit_count$Chat, fit_plain$Chat, tolerance = 1e-12)
  # The e1/e3 pair is observed by nobody, so it stays exactly zero.
  expect_equal(unname(oracle[1, 3]), 0)
})

test_that("rescale/shrink divide by subject-weight mass under mfa_sigma1", {
  effects <- c("e1", "e2", "e3")
  design <- diag(3)
  colnames(design) <- effects
  set.seed(17)
  rows <- list(1:3, 1:2, 2:3)
  scale <- c(1, 1, 8)
  subjects <- lapply(seq_along(rows), function(i) {
    B <- matrix(rnorm(12), 3, 4, dimnames = list(effects, NULL)) * scale[[i]]
    B[-rows[[i]], ] <- 0
    suppressWarnings(dkge_subject(B, design, id = paste0("s", i),
                                  observed_rows = rows[[i]]))
  })
  data <- dkge_data(subjects)
  K <- diag(3)
  dimnames(K) <- list(effects, effects)

  fit_resc <- dkge_fit(data, K = K, rank = 2, w_method = "mfa_sigma1",
                       effect_scaling = "none", missingness = "rescale")
  w <- fit_resc$weights
  expect_false(isTRUE(all.equal(w, rep(mean(w), 3))))
  masks <- lapply(data$provenance$obs_mask, as.numeric)
  structural <- lapply(masks, tcrossprod)
  moments <- lapply(data$betas, tcrossprod)
  num <- Reduce(`+`, Map(function(ws, M) ws * M, as.list(w), moments))
  den <- Reduce(`+`, Map(function(ws, S) ws * S, as.list(w), structural))
  oracle <- matrix(0, 3, 3)
  oracle[den > 0] <- num[den > 0] / den[den > 0]
  oracle <- (oracle + t(oracle)) / 2
  expect_equal(unname(fit_resc$effect_moment), unname(oracle),
               tolerance = 1e-12)

  fit_shrink <- dkge_fit(data, K = K, rank = 2, w_method = "mfa_sigma1",
                         effect_scaling = "none", missingness = "shrink",
                         miss_args = list(gamma = 1))
  pc <- Reduce(`+`, structural)
  rel <- pc / max(pc)
  shrink_oracle <- rel * oracle + (1 - rel) * diag(diag(oracle), 3, 3)
  shrink_oracle <- (shrink_oracle + t(shrink_oracle)) / 2
  expect_equal(unname(fit_shrink$effect_moment), unname(shrink_oracle),
               tolerance = 1e-12)
})

test_that("unweighted and precision-weighted rescale agree under unit precision", {
  set.seed(11)
  effects <- paste0("e", 1:3)
  subjects <- lapply(1:3, function(s) {
    B <- matrix(rnorm(12), 3, 4, dimnames = list(effects, NULL))
    make_effect_weight_subject(B, n = stats::setNames(rep(1, 3), effects),
                               id = paste0("s", s))
  })
  data <- dkge_data(subjects)
  for (miss in c("none", "rescale", "mask", "shrink")) {
    fit_plain <- dkge_fit(data, K = diag(3), rank = 2, w_method = "none",
                          effect_scaling = "none", missingness = miss,
                          miss_args = list(min_pairs = 2, gamma = 1))
    fit_count <- dkge_fit(data, K = diag(3), rank = 2, w_method = "none",
                          effect_scaling = "none",
                          effect_weights = dkge_effect_weights("count"),
                          missingness = miss,
                          miss_args = list(min_pairs = 2, gamma = 1))
    expect_equal(fit_count$Chat, fit_plain$Chat, tolerance = 1e-12, info = miss)
  }
})

test_that("unit-precision count weights do not inflate sparsely observed pairs", {
  fx_plain <- make_partial_fit(seed = 1L, w_method = "none",
                               effect_scaling = "none")
  fx_count <- make_partial_fit(seed = 1L, w_method = "none",
                               effect_scaling = "none",
                               effect_weights = dkge_effect_weights("count"))
  plain <- fx_plain$fit$effect_moment
  count <- fx_count$fit$effect_moment
  expect_equal(unname(count), unname(plain), tolerance = 1e-12)

  # The e4 row/col is observed by 1 of 6 subjects. The old cohort-mass
  # restoration inflated those entries by S/k = 6.
  e4 <- 4L
  expect_equal(unname(count[e4, ] / plain[e4, ]), rep(1, 4), tolerance = 1e-12)
})

test_that("unit-precision count weights agree on a partial-coverage miss loop", {
  fx <- make_partial_fit(seed = 2L, w_method = "none", effect_scaling = "none")
  data <- fx$data
  K <- diag(length(fx$effects))
  dimnames(K) <- list(fx$effects, fx$effects)
  for (miss in c("none", "rescale", "mask", "shrink")) {
    fit_plain <- dkge_fit(data, K = K, rank = 2, w_method = "none",
                          effect_scaling = "none", missingness = miss,
                          miss_args = list(min_pairs = 2, gamma = 1))
    fit_count <- dkge_fit(data, K = K, rank = 2, w_method = "none",
                          effect_scaling = "none",
                          effect_weights = dkge_effect_weights("count"),
                          missingness = miss,
                          miss_args = list(min_pairs = 2, gamma = 1))
    expect_equal(unname(fit_count$effect_moment), unname(fit_plain$effect_moment),
                 tolerance = 1e-12, info = miss)
  }
})

test_that("effect weighting combines with partial effect coverage", {
  effects <- c("e1", "e2", "e3")
  design <- diag(3)
  colnames(design) <- effects
  make_one <- function(id, rows, prec) {
    B <- matrix(c(1, 2, 3, 4, 5, 6), 3, 2, dimnames = list(effects, NULL))
    B[-rows, ] <- 0
    suppressWarnings(dkge_subject(
      B, design, id = id, observed_rows = rows,
      effect_precision = stats::setNames(prec, effects)
    ))
  }
  data <- dkge_data(list(
    make_one("s1", c(1L, 2L), c(2, 3, 0)),
    make_one("s2", c(2L, 3L), c(0, 1, 4)),
    make_one("s3", 1:3, c(1, 1, 1))
  ))
  fit <- dkge_fit(data, K = diag(3), rank = 2, w_method = "none",
                  effect_scaling = "none",
                  effect_weights = dkge_effect_weights("precision"))

  # Hand-computed pooled entry for the pair (e1, e3), observed only by s3.
  # Scale is observed subject mass, not the full cohort (S/k inflation).
  h <- lapply(seq_len(3), function(s) {
    p <- data$effect_precision[[s]]
    p[!data$provenance$obs_mask[[s]]] <- 0
    tcrossprod(sqrt(p))
  })
  M <- lapply(data$betas, tcrossprod)
  num <- Reduce(`+`, Map(`*`, h, M))
  den <- Reduce(`+`, h)
  obs_mass <- fit$weights[[3]]
  expect_equal(unname(fit$effect_moment[1, 3]),
               obs_mass * num[1, 3] / den[1, 3], tolerance = 1e-12)

  # Pairs the precision zeroes out carry no mass.
  expect_equal(unname(fit$pair_ess[1, 3]), 1, tolerance = 1e-12)
  expect_gt(fit$pair_ess[2, 2], 2)
})

test_that("an all-unobserved mask yields a zero effect moment, not the unmasked one", {
  B <- matrix(c(1, 2, 3, 4, 5, 6), 3, 2,
              dimnames = list(c("e1", "e2", "e3"), NULL))
  zero <- dkge:::.dkge_effect_moment(B, obs_mask = rep(FALSE, 3))
  expect_equal(zero, matrix(0, 3, 3,
                            dimnames = list(rownames(B), rownames(B))))
  expect_false(isTRUE(all.equal(zero, dkge:::.dkge_effect_moment(B))))

  Lambda <- diag(3)
  dimnames(Lambda) <- list(rownames(B), rownames(B))
  subject <- list(id = "s", effect_noise_cov = Lambda, noise_trace = 5)
  expect_equal(
    dkge:::.dkge_effect_noise_moment(subject, obs_mask = rep(FALSE, 3)),
    matrix(0, 3, 3, dimnames = dimnames(Lambda))
  )

  splits <- list(B, B)
  expect_equal(
    dkge:::.dkge_split_effect_moment(splits, obs_mask = rep(FALSE, 3)),
    matrix(0, 3, 3, dimnames = list(rownames(B), rownames(B)))
  )
})
