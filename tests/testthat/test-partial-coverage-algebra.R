library(testthat)

make_partial_coverage_pair <- function(omega = NULL, reorder_local = FALSE) {
  effects <- paste0("e", 1:4)
  B1 <- matrix(c(
    1, 2, 3,
    2, 4, 1,
    999, -777, 555,
    3, 1, 5
  ), 4, 3, byrow = TRUE, dimnames = list(effects, paste0("v", 1:3)))
  X1 <- matrix(c(
    1, 0, 91, 0,
    0, 1, 92, 0,
    1, 1, 93, 0,
    2, 1, 94, 1,
    1, 2, 95, 2,
    3, 1, 96, 1
  ), 6, 4, byrow = TRUE, dimnames = list(NULL, effects))

  local_effects <- c("e2", "e3", "e4")
  B2 <- matrix(c(
    4, 2, 1,
    1, 3, 2,
    5, 2, 4
  ), 3, 3, byrow = TRUE,
  dimnames = list(local_effects, paste0("v", 1:3)))
  X2 <- matrix(c(
    1, 0, 0,
    0, 1, 0,
    1, 1, 0,
    2, 1, 1,
    1, 2, 1,
    1, 1, 2
  ), 6, 3, byrow = TRUE, dimnames = list(NULL, local_effects))
  if (reorder_local) {
    local_order <- c("e4", "e2", "e3")
    B2 <- B2[local_order, , drop = FALSE]
    X2 <- X2[, local_order, drop = FALSE]
  }

  s1 <- suppressWarnings(dkge_subject(
    B1, X1, id = "s1", omega = omega, observed_rows = c(1, 2, 4)
  ))
  s2_local <- suppressWarnings(dkge_subject(B2, X2, id = "s2", omega = omega))
  heterogeneous <- dkge_data(list(s1, s2_local))

  B2_global <- matrix(8e5, 4, 3,
                      dimnames = list(effects, paste0("v", 1:3)))
  B2_global[match(rownames(B2), effects), ] <- B2
  X2_global <- matrix(-6e5, 6, 4, dimnames = list(NULL, effects))
  X2_global[, match(colnames(X2), effects)] <- X2
  s2_global <- suppressWarnings(dkge_subject(
    B2_global, X2_global, id = "s2", omega = omega,
    observed_rows = match(rownames(B2), effects)
  ))
  identical_labels <- dkge_data(list(s1, s2_global))

  list(heterogeneous = heterogeneous, identical = identical_labels)
}

partial_kernel <- function(effects) {
  K <- matrix(c(
    1.4, 0.3, 0.2, 0.1,
    0.3, 1.5, 0.4, 0.2,
    0.2, 0.4, 1.3, 0.35,
    0.1, 0.2, 0.35, 1.2
  ), 4, 4, byrow = TRUE)
  dimnames(K) <- list(effects, effects)
  K
}

fit_partial_pair <- function(pair, K) {
  spec <- dkge_weights(prior = c(0.7, 1.1, 1.6), adapt = "none")
  lapply(pair, function(data) {
    dkge_fit(
      data, K = K, rank = 2, w_method = "energy", w_tau = 0,
      weights = spec, missingness = "none", keep_X = TRUE
    )
  })
}

test_that("identical-label masks and heterogeneous local effects are equivalent", {
  omega_cases <- list(
    diagonal = c(0.5, 1.2, 1.8),
    full = matrix(c(1.2, 0.2, 0.1,
                    0.2, 1.5, 0.25,
                    0.1, 0.25, 1.1), 3, 3, byrow = TRUE)
  )

  for (omega in omega_cases) {
    pair <- make_partial_coverage_pair(omega)
    K <- partial_kernel(pair$heterogeneous$effects)
    fits <- fit_partial_pair(pair, K)
    a <- fits$heterogeneous
    b <- fits$identical

    expect_equal(a$R, b$R, tolerance = 1e-12, ignore_attr = TRUE)
    expect_equal(a$weights, b$weights, tolerance = 1e-12)
    expect_equal(a$Chat, b$Chat, tolerance = 1e-11, ignore_attr = TRUE)
    expect_equal(a$contribs, b$contribs, tolerance = 1e-11,
                 ignore_attr = TRUE)
    expect_equal(a$X_concat, b$X_concat, tolerance = 1e-11,
                 ignore_attr = TRUE)

    ctx_a <- dkge:::.dkge_fold_weight_context(a, train_ids = 1:2)
    ctx_b <- dkge:::.dkge_fold_weight_context(b, train_ids = 1:2)
    expect_equal(ctx_a$subject_weights, ctx_b$subject_weights,
                 tolerance = 1e-12)
    expect_equal(ctx_a$Chat, ctx_b$Chat, tolerance = 1e-11,
                 ignore_attr = TRUE)
  }
})

test_that("unobserved garbage cannot affect ruler, weights, moment, or folds", {
  pair <- make_partial_coverage_pair(
    matrix(c(1.2, 0.2, 0.1,
             0.2, 1.5, 0.25,
             0.1, 0.25, 1.1), 3, 3, byrow = TRUE)
  )
  clean <- pair$identical
  dirty <- clean
  dirty$betas[[1]][3, ] <- c(1e12, -2e12, 3e12)
  dirty$designs[[1]][, 3] <- seq(1e9, 6e9, length.out = 6)
  dirty$subjects[[1]]$beta <- dirty$betas[[1]]
  dirty$subjects[[1]]$design <- dirty$designs[[1]]
  dirty$betas[[2]][1, ] <- c(-4e11, 5e11, -6e11)
  dirty$designs[[2]][, 1] <- seq(-6e8, -1e8, length.out = 6)
  dirty$subjects[[2]]$beta <- dirty$betas[[2]]
  dirty$subjects[[2]]$design <- dirty$designs[[2]]

  K <- partial_kernel(clean$effects)
  clean_fit <- fit_partial_pair(list(clean = clean), K)$clean
  dirty_fit <- fit_partial_pair(list(dirty = dirty), K)$dirty

  expect_equal(dirty_fit$R, clean_fit$R, tolerance = 1e-12,
               ignore_attr = TRUE)
  expect_equal(dirty_fit$weights, clean_fit$weights, tolerance = 1e-12)
  expect_equal(dirty_fit$Chat, clean_fit$Chat, tolerance = 1e-11,
               ignore_attr = TRUE)
  expect_equal(dirty_fit$U %*% t(dirty_fit$U),
               clean_fit$U %*% t(clean_fit$U), tolerance = 1e-9,
               ignore_attr = TRUE)
  expect_equal(dirty_fit$Btil, clean_fit$Btil, tolerance = 1e-12,
               ignore_attr = TRUE)

  dirty_ctx <- dkge:::.dkge_fold_weight_context(dirty_fit, train_ids = 1:2)
  clean_ctx <- dkge:::.dkge_fold_weight_context(clean_fit, train_ids = 1:2)
  expect_equal(dirty_ctx$Chat, clean_ctx$Chat, tolerance = 1e-11,
               ignore_attr = TRUE)

  assignments <- list(1L, 2L)
  dirty_folds <- dkge:::.dkge_build_fold_bases(
    dirty_fit, assignments, align = FALSE, loader_scope = "all"
  )
  clean_folds <- dkge:::.dkge_build_fold_bases(
    clean_fit, assignments, align = FALSE, loader_scope = "all"
  )
  for (f in seq_along(assignments)) {
    expect_equal(dirty_folds$folds[[f]]$basis %*%
                   t(dirty_folds$folds[[f]]$basis),
                 clean_folds$folds[[f]]$basis %*%
                   t(clean_folds$folds[[f]]$basis),
                 tolerance = 1e-9, ignore_attr = TRUE)
    expect_equal(dirty_folds$folds[[f]]$loaders,
                 clean_folds$folds[[f]]$loaders,
                 tolerance = 1e-10, ignore_attr = TRUE)
  }
})

test_that("disjoint missing blocks follow a raw-effect pair-count oracle", {
  make_local <- function(effects, offset, id) {
    B <- matrix(offset + seq_len(6), 2, 3,
                dimnames = list(effects, paste0("v", 1:3)))
    X <- matrix(c(1, 0, 0, 1, 1, 1, 2, 1), 4, 2, byrow = TRUE,
                dimnames = list(NULL, effects))
    suppressWarnings(dkge_subject(B, X, id = id))
  }
  data <- dkge_data(list(
    make_local(c("e1", "e2"), 0, "s1"),
    make_local(c("e3", "e4"), 10, "s2")
  ))
  K <- partial_kernel(data$effects)
  fit <- dkge_fit(data, K = K, rank = 2, w_method = "none",
                  missingness = "rescale")

  masks <- dkge:::.dkge_obs_masks_from_provenance(
    fit$provenance, fit$subject_ids, 4
  )
  moments <- lapply(seq_along(fit$Braw), function(s) {
    B <- fit$Braw[[s]]
    B[!masks[[s]], ] <- 0
    tcrossprod(B)
  })
  numerator <- Reduce(`+`, moments)
  pair_counts <- Reduce(`+`, lapply(masks, function(mask) {
    tcrossprod(as.numeric(mask))
  }))
  pooled <- numerator / pmax(pair_counts, 1)
  pooled[pair_counts == 0] <- 0
  oracle <- fit$Khalf %*% t(fit$R) %*% pooled %*% fit$R %*% fit$Khalf
  oracle <- (oracle + t(oracle)) / 2

  expect_equal(fit$effect_moment, pooled, tolerance = 1e-12,
               ignore_attr = TRUE)
  expect_equal(fit$Chat, oracle, tolerance = 1e-11, ignore_attr = TRUE)
  expect_equal(pair_counts[1:2, 3:4], matrix(0, 2, 2))
})

test_that("fold rebuilds recompute subject energy on the same masked blocks", {
  Omega <- matrix(c(1.2, 0.2, 0.1,
                    0.2, 1.5, 0.25,
                    0.1, 0.25, 1.1), 3, 3, byrow = TRUE)
  pair <- make_partial_coverage_pair(Omega)
  subjects <- pair$identical$subjects
  s3 <- subjects[[1]]
  s3$id <- "s3"
  s3$beta <- 1.7 * s3$beta
  data <- dkge_data(c(subjects, list(s3)))
  K <- partial_kernel(data$effects)
  spec <- dkge_weights(prior = c(0.7, 1.1, 1.6), adapt = "none")
  fit <- dkge_fit(data, K = K, rank = 2, w_method = "energy", w_tau = 0,
                  weights = spec, missingness = "none")

  train <- c(1, 3)
  ctx <- dkge:::.dkge_fold_weight_context(fit, train_ids = train)
  masks <- dkge:::.dkge_obs_masks_from_provenance(
    fit$provenance, fit$subject_ids, nrow(K)
  )
  eo <- eigen(Omega, symmetric = TRUE)
  Ohalf <- eo$vectors %*% diag(sqrt(eo$values), 3) %*% t(eo$vectors)
  voxel_weight <- ctx$weights$total
  blocks <- lapply(train, function(s) {
    B <- fit$Braw[[s]]
    B[!masks[[s]], ] <- 0
    B <- sweep(B, 2, sqrt(voxel_weight), "*") %*% Ohalf
    fit$Khalf %*% t(ctx$R) %*% B
  })
  raw_weights <- vapply(blocks, function(X) 1 / (sum(X * X) + 1e-12),
                          numeric(1))
  expected_weights <- raw_weights / mean(raw_weights)
  expected_chat <- Reduce(`+`, Map(function(X, w) w * tcrossprod(X),
                                   blocks, expected_weights))

  expect_equal(ctx$subject_weights, expected_weights, tolerance = 1e-11)
  expect_equal(ctx$Chat, expected_chat, tolerance = 1e-10,
               ignore_attr = TRUE)
})

test_that("full-coverage Chat matches an independent dense factor oracle", {
  set.seed(7301)
  q <- 4
  P <- 3
  effects <- paste0("e", 1:q)
  betas <- lapply(1:3, function(s) {
    matrix(rnorm(q * P), q, P,
           dimnames = list(effects, paste0("v", 1:P)))
  })
  designs <- lapply(1:3, function(s) {
    X <- matrix(rnorm(8 * q), 8, q)
    colnames(X) <- effects
    X
  })
  Omega <- list(
    c(0.7, 1.2, 1.5),
    matrix(c(1.3, 0.2, 0.1,
             0.2, 1.1, 0.15,
             0.1, 0.15, 1.4), 3, 3, byrow = TRUE),
    NULL
  )
  A <- matrix(rnorm(q * q), q, q)
  K <- crossprod(A) + 0.4 * diag(q)
  dimnames(K) <- list(effects, effects)
  fit <- suppressWarnings(dkge_fit(
    betas, designs, K = K, Omega_list = Omega, rank = 3,
    w_method = "none", keep_X = TRUE
  ))

  G <- Reduce(`+`, lapply(designs, crossprod))
  diag(G) <- diag(G) + 1e-10
  R <- chol(G)
  ek <- eigen((K + t(K)) / 2, symmetric = TRUE)
  Khalf <- ek$vectors %*% diag(sqrt(ek$values), q) %*% t(ek$vectors)
  right_root <- function(O) {
    if (is.null(O)) return(diag(P))
    if (is.vector(O)) return(diag(sqrt(O), P))
    eo <- eigen((O + t(O)) / 2, symmetric = TRUE)
    eo$vectors %*% diag(sqrt(pmax(eo$values, 0)), P) %*% t(eo$vectors)
  }
  Xstar <- do.call(cbind, lapply(seq_along(betas), function(s) {
    Khalf %*% t(R) %*% betas[[s]] %*% right_root(Omega[[s]])
  }))
  rel_error <- norm(fit$Chat - tcrossprod(Xstar), "F") /
    max(1, norm(fit$Chat, "F"))

  expect_lt(rel_error, 1e-12)
  expect_equal(fit$X_concat, Xstar, tolerance = 1e-11,
               ignore_attr = TRUE)
})

test_that("subject and local-effect ordering do not change the pooled result", {
  standard <- make_partial_coverage_pair(c(0.8, 1.1, 1.4), reorder_local = FALSE)
  reordered <- make_partial_coverage_pair(c(0.8, 1.1, 1.4), reorder_local = TRUE)
  K <- partial_kernel(standard$heterogeneous$effects)
  fit_standard <- fit_partial_pair(list(x = standard$heterogeneous), K)$x
  fit_local_reordered <- fit_partial_pair(list(x = reordered$heterogeneous), K)$x

  global_reversed <- standard$identical
  global_reversed$subjects <- rev(global_reversed$subjects)
  global_reversed$betas <- rev(global_reversed$betas)
  global_reversed$designs <- rev(global_reversed$designs)
  global_reversed$omega <- rev(global_reversed$omega)
  global_reversed$subject_ids <- rev(global_reversed$subject_ids)
  global_reversed$effect_n <- rev(global_reversed$effect_n)
  global_reversed$effect_precision <- rev(global_reversed$effect_precision)
  global_reversed$effect_noise_cov <- rev(global_reversed$effect_noise_cov)
  global_reversed$residual_variance <- rev(global_reversed$residual_variance)
  global_reversed$residual_df <- rev(global_reversed$residual_df)
  global_reversed$noise_trace <- rev(global_reversed$noise_trace)
  global_reversed$split_betas <- rev(global_reversed$split_betas)
  global_reversed$provenance$obs_mask <- rev(global_reversed$provenance$obs_mask)
  names(global_reversed$provenance$obs_mask) <- global_reversed$subject_ids
  fit_reversed <- fit_partial_pair(list(x = global_reversed), K)$x

  expect_equal(fit_local_reordered$Chat, fit_standard$Chat,
               tolerance = 1e-11, ignore_attr = TRUE)
  expect_equal(fit_reversed$Chat, fit_standard$Chat,
               tolerance = 1e-11, ignore_attr = TRUE)
  expect_equal(sort(fit_reversed$weights), sort(fit_standard$weights),
               tolerance = 1e-12)
})
