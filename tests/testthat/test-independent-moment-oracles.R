library(testthat)

make_combined_oracle_fixture <- function() {
  global_effects <- c("e1", "e2", "e3")
  local_effects <- list(global_effects, c("e1", "e2"), c("e2", "e3"))
  betas <- list(
    matrix(c(
      2.0, 1.0, -0.5,
      0.5, 1.5, 2.0,
      -1.0, 0.25, 1.25
    ), 3, 3, byrow = TRUE),
    matrix(c(
      1.5, -0.5, 0.75,
      2.0, 1.0, -1.0
    ), 2, 3, byrow = TRUE),
    matrix(c(
      0.75, 1.25, -0.25,
      1.0, -1.5, 2.5
    ), 2, 3, byrow = TRUE)
  )
  designs <- list(
    matrix(c(
      1, 0, 1,
      1, 1, 0,
      0, 1, 1,
      2, 1, 0,
      0, 1, 2
    ), 5, 3, byrow = TRUE),
    matrix(c(
      1, 0,
      1, 1,
      0, 1,
      2, 1,
      1, 2
    ), 5, 2, byrow = TRUE),
    matrix(c(
      1, 0,
      1, 1,
      0, 1,
      2, 1,
      1, 3
    ), 5, 2, byrow = TRUE)
  )
  omega <- list(
    matrix(c(1.4, 0.2, 0.1,
             0.2, 1.1, 0.15,
             0.1, 0.15, 0.9), 3, 3),
    matrix(c(1.2, 0.1, 0.05,
             0.1, 0.8, 0.12,
             0.05, 0.12, 1.3), 3, 3),
    matrix(c(0.9, 0.08, 0.03,
             0.08, 1.5, 0.2,
             0.03, 0.2, 1.1), 3, 3)
  )
  precision <- list(
    c(e1 = 1, e2 = 4, e3 = 2),
    c(e1 = 9, e2 = 1),
    c(e2 = 3, e3 = 7)
  )
  residual_variance <- list(
    c(0.12, 0.18, 0.15),
    c(0.20, 0.16, 0.14),
    c(0.10, 0.22, 0.17)
  )
  effect_covariance <- list(
    matrix(c(0.08, 0.01, 0.005,
             0.01, 0.06, 0.008,
             0.005, 0.008, 0.07), 3, 3),
    matrix(c(0.09, 0.012,
             0.012, 0.05), 2, 2),
    matrix(c(0.07, 0.01,
             0.01, 0.08), 2, 2)
  )
  subjects <- lapply(seq_len(3), function(s) {
    rownames(betas[[s]]) <- local_effects[[s]]
    colnames(betas[[s]]) <- paste0("v", seq_len(3))
    colnames(designs[[s]]) <- local_effects[[s]]
    dimnames(effect_covariance[[s]]) <- list(local_effects[[s]],
                                              local_effects[[s]])
    suppressWarnings(dkge_subject(
      betas[[s]], designs[[s]], id = paste0("s", s),
      omega = omega[[s]], effect_precision = precision[[s]],
      effect_noise_cov = effect_covariance[[s]],
      residual_variance = residual_variance[[s]], residual_df = 5
    ))
  })
  K <- matrix(c(
    1.5, 0.3, 0.15,
    0.3, 1.2, 0.25,
    0.15, 0.25, 1.1
  ), 3, 3, dimnames = list(global_effects, global_effects))
  list(
    subjects = subjects,
    effects = global_effects,
    K = K,
    local_effects = local_effects
  )
}

fit_combined_oracle_fixture <- function(fixture) {
  dkge_fit(
    dkge_data(fixture$subjects), K = fixture$K, rank = 2,
    w_method = "none", effect_scaling = "pooled_design",
    effect_weights = dkge_effect_weights("precision"),
    debias = "analytic"
  )
}

dense_combined_oracle <- function(fixture) {
  effects <- fixture$effects
  q <- length(effects)
  S <- length(fixture$subjects)
  embed_rows <- function(x, local, columns = ncol(x)) {
    out <- matrix(0, q, columns,
                  dimnames = list(effects, colnames(x)))
    out[match(local, effects), ] <- x
    out
  }
  embed_square <- function(x, local) {
    out <- matrix(0, q, q, dimnames = list(effects, effects))
    idx <- match(local, effects)
    out[idx, idx] <- x
    out
  }

  G <- matrix(0, q, q)
  moments <- vector("list", S)
  noise <- vector("list", S)
  precision <- vector("list", S)
  masks <- vector("list", S)
  B_global <- vector("list", S)
  for (s in seq_len(S)) {
    subject <- fixture$subjects[[s]]
    local <- subject$effects
    X_global <- matrix(0, nrow(subject$design), q)
    X_global[, match(local, effects)] <- subject$design
    G <- G + crossprod(X_global)

    B <- embed_rows(subject$beta, local)
    Lambda <- embed_square(subject$effect_noise_cov, local)
    trace <- sum(diag(subject$omega) * subject$residual_variance)
    raw <- B %*% subject$omega %*% t(B)
    noise[[s]] <- trace * Lambda
    moments[[s]] <- (raw - noise[[s]] + t(raw - noise[[s]])) / 2
    p <- setNames(numeric(q), effects)
    p[local] <- subject$effect_precision
    precision[[s]] <- p
    masks[[s]] <- p > 0
    B_global[[s]] <- B
  }
  diag(G) <- diag(G) + 1e-10
  R <- chol(G)

  numerator <- matrix(0, q, q)
  pair_weight <- matrix(0, q, q)
  pair_weight_sq <- matrix(0, q, q)
  pair_counts <- matrix(0, q, q)
  for (s in seq_len(S)) {
    structural <- tcrossprod(as.numeric(masks[[s]]))
    h <- tcrossprod(sqrt(precision[[s]])) * structural
    numerator <- numerator + h * moments[[s]]
    pair_weight <- pair_weight + h
    pair_weight_sq <- pair_weight_sq + h^2
    pair_counts <- pair_counts + structural
  }
  pooled <- matrix(0, q, q)
  valid <- pair_weight > 0
  pooled[valid] <- S * numerator[valid] / pair_weight[valid]
  pair_ess <- matrix(0, q, q)
  pair_ess[valid] <- pair_weight[valid]^2 / pair_weight_sq[valid]

  kernel_eigen <- eigen(fixture$K, symmetric = TRUE)
  Khalf <- kernel_eigen$vectors %*%
    diag(sqrt(kernel_eigen$values), q) %*%
    t(kernel_eigen$vectors)
  Kihalf <- kernel_eigen$vectors %*%
    diag(1 / sqrt(kernel_eigen$values), q) %*%
    t(kernel_eigen$vectors)
  Chat <- Khalf %*% t(R) %*% pooled %*% R %*% Khalf
  Chat <- (Chat + t(Chat)) / 2
  eig <- eigen(Chat, symmetric = TRUE)
  U <- Kihalf %*% eig$vectors[, 1:2, drop = FALSE]

  list(
    G = G, R = R, Khalf = Khalf, Kihalf = Kihalf,
    B = B_global, noise = noise, moments = moments,
    pooled = pooled, pair_counts = pair_counts,
    pair_weight = pair_weight, pair_ess = pair_ess,
    Chat = Chat, eig = eig, U = U
  )
}

keigen_sqrt <- function(x) {
  eig <- eigen((x + t(x)) / 2, symmetric = TRUE)
  eig$vectors %*% diag(sqrt(pmax(eig$values, 0)), nrow(x)) %*%
    t(eig$vectors)
}

test_that("combined moment obeys the full dense q-space contract", {
  fixture <- make_combined_oracle_fixture()
  fit <- fit_combined_oracle_fixture(fixture)
  oracle <- dense_combined_oracle(fixture)

  expect_equal(fit$R, oracle$R, tolerance = 1e-11,
               ignore_attr = TRUE,
               info = "pooled ruler must equal chol(sum X_s'X_s + jitter I)")
  expect_equal(fit$noise_moments, oracle$noise, tolerance = 1e-11,
               ignore_attr = TRUE,
               info = "analytic subtraction must use tr(Omega Sigma) Lambda")
  expect_equal(fit$effect_moments, oracle$moments, tolerance = 1e-11,
               ignore_attr = TRUE,
               info = "raw full-Omega moments must be debiased before pooling")
  expect_equal(fit$effect_moment, oracle$pooled, tolerance = 1e-11,
               ignore_attr = TRUE,
               info = "unequal effect precision requires pair-normalized pooling")
  expect_equal(fit$pair_counts, oracle$pair_counts, tolerance = 1e-12,
               ignore_attr = TRUE,
               info = "pair counts must reflect effective observed precision")
  expect_equal(fit$pair_weight, oracle$pair_weight, tolerance = 1e-12,
               ignore_attr = TRUE,
               info = "pair weights must use geometric-mean precision")
  expect_equal(fit$pair_ess, oracle$pair_ess, tolerance = 1e-12,
               ignore_attr = TRUE,
               info = "pair ESS must be sum(w)^2/sum(w^2)")
  expect_equal(fit$Chat, oracle$Chat, tolerance = 1e-10,
               ignore_attr = TRUE,
               info = "kernel/ruler transform must occur after raw-effect pooling")
  expect_equal(fit$eig_values_full, oracle$eig$values, tolerance = 1e-10,
               ignore_attr = TRUE,
               info = "stored spectrum must be the dense Chat eigenspectrum")
  expect_equal(
    tcrossprod(fit$Khalf %*% fit$U),
    tcrossprod(oracle$Khalf %*% oracle$U),
    tolerance = 1e-9,
    ignore_attr = TRUE,
    info = "K-whitened leading eigenspaces must agree"
  )

  manual_loadings <- lapply(oracle$B, function(B) {
    t(B) %*% oracle$R %*% fixture$K %*% oracle$U
  })
  observed_loadings <- dkge_predict_loadings(fit, fit$Braw)
  expect_equal(observed_loadings, manual_loadings, tolerance = 1e-9,
               ignore_attr = TRUE,
               info = "subject projection must be B_raw' R K U")
})

test_that("combined estimator is invariant to subject and effect order", {
  fixture <- make_combined_oracle_fixture()
  fit <- fit_combined_oracle_fixture(fixture)

  subject_permutation <- c(3, 1, 2)
  reordered <- fixture
  reordered$subjects <- fixture$subjects[subject_permutation]
  fit_subject <- fit_combined_oracle_fixture(reordered)
  expect_equal(fit_subject$Chat, fit$Chat, tolerance = 1e-10,
               ignore_attr = TRUE,
               info = "subject order cannot change the pooled coordinate")
  expect_equal(fit_subject$pair_ess, fit$pair_ess, tolerance = 1e-12,
               ignore_attr = TRUE)

  effect_permutation <- c(3, 1, 2)
  permuted_effects <- fixture$effects[effect_permutation]
  permuted <- fixture
  permuted$effects <- permuted_effects
  permuted$K <- fixture$K[permuted_effects, permuted_effects]
  permuted$subjects <- lapply(fixture$subjects, function(subject) {
    order_local <- permuted_effects[permuted_effects %in% subject$effects]
    dkge_subject(
      subject$beta[order_local, , drop = FALSE],
      subject$design[, order_local, drop = FALSE],
      id = subject$id, omega = subject$omega,
      effect_precision = subject$effect_precision[order_local],
      effect_noise_cov = subject$effect_noise_cov[order_local, order_local,
                                                  drop = FALSE],
      residual_variance = subject$residual_variance,
      residual_df = subject$residual_df
    )
  })
  fit_effect <- suppressWarnings(fit_combined_oracle_fixture(permuted))
  expect_equal(
    fit_effect$Chat,
    fit$Chat,
    tolerance = 1e-10,
    ignore_attr = TRUE,
    info = "local effect-row permutations cannot change a labelled fit"
  )
})

test_that("precision scale, duplication, and feature units obey metamorphic laws", {
  fixture <- make_combined_oracle_fixture()
  fit <- fit_combined_oracle_fixture(fixture)

  precision_scaled <- fixture
  precision_scaled$subjects <- lapply(fixture$subjects, function(subject) {
    subject$effect_precision <- 37 * subject$effect_precision
    subject
  })
  fit_precision <- fit_combined_oracle_fixture(precision_scaled)
  expect_equal(fit_precision$Chat, fit$Chat, tolerance = 1e-10,
               ignore_attr = TRUE,
               info = "common precision scaling must cancel in pair normalization")
  expect_equal(fit_precision$pair_ess, fit$pair_ess, tolerance = 1e-12,
               ignore_attr = TRUE)

  duplicated <- fixture
  duplicated$subjects <- c(fixture$subjects, lapply(fixture$subjects, function(subject) {
    subject$id <- paste0(subject$id, "_copy")
    subject
  }))
  fit_duplicated <- fit_combined_oracle_fixture(duplicated)
  expect_equal(fit_duplicated$effect_moment, 2 * fit$effect_moment,
               tolerance = 1e-10, ignore_attr = TRUE,
               info = "literal cohort duplication doubles the pooled raw moment")
  expect_equal(fit_duplicated$Chat, 4 * fit$Chat, tolerance = 1e-9,
               ignore_attr = TRUE,
               info = "cohort duplication doubles both raw moment and pooled ruler mass")
  expect_equal(fit_duplicated$pair_ess, 2 * fit$pair_ess,
               tolerance = 1e-12, ignore_attr = TRUE)

  unit_changed <- fixture
  unit_changed$subjects <- lapply(seq_along(fixture$subjects), function(s) {
    subject <- fixture$subjects[[s]]
    scale <- c(0.5, 2, 3 + s)
    D <- diag(scale)
    Dinv <- diag(1 / scale)
    subject$beta <- subject$beta %*% D
    subject$omega <- Dinv %*% subject$omega %*% Dinv
    subject$residual_variance <- subject$residual_variance * scale^2
    subject
  })
  fit_units <- fit_combined_oracle_fixture(unit_changed)
  expect_equal(fit_units$Chat, fit$Chat, tolerance = 1e-9,
               ignore_attr = TRUE,
               info = "congruent feature-unit changes must preserve signal and noise moments")
  expect_equal(fit_units$effect_moment, fit$effect_moment,
               tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("full coverage remains identical to the pre-tier dense formula", {
  set.seed(8519)
  S <- 3L
  q <- 3L
  P <- 4L
  effects <- paste0("e", seq_len(q))
  betas <- lapply(seq_len(S), function(s) {
    matrix(rnorm(q * P), q, P, dimnames = list(effects, NULL))
  })
  designs <- lapply(seq_len(S), function(s) {
    X <- matrix(rnorm(7 * q), 7, q)
    colnames(X) <- effects
    X
  })
  omega <- lapply(seq_len(S), function(s) {
    A <- matrix(rnorm(P * P), P, P)
    crossprod(A) + diag(P)
  })
  Kbase <- matrix(rnorm(q * q), q, q)
  K <- crossprod(Kbase) + diag(q)
  fit <- dkge_fit(
    betas, designs, K = K, Omega_list = omega,
    rank = 2, w_method = "none"
  )

  G <- Reduce(`+`, lapply(designs, crossprod))
  diag(G) <- diag(G) + 1e-10
  R <- chol(G)
  ke <- eigen(K, symmetric = TRUE)
  Khalf <- ke$vectors %*% diag(sqrt(ke$values), q) %*% t(ke$vectors)
  dense_blocks <- lapply(seq_len(S), function(s) {
    Khalf %*% t(R) %*% betas[[s]] %*%
      keigen_sqrt(omega[[s]])
  })
  Xstar <- do.call(cbind, dense_blocks)

  expect_equal(fit$Chat, tcrossprod(Xstar), tolerance = 1e-9,
               ignore_attr = TRUE,
               info = "ordinary full coverage must retain the legacy block formula")
})
