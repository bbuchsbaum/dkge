library(testthat)

line_coords <- function(P, labels = NULL) {
  x <- cbind(x = seq_len(P) - 1, y = 0, z = 0)
  if (!is.null(labels)) rownames(x) <- labels
  x
}

line_spatial <- function(P, lambda, labels = NULL) {
  dkge_spatial_regularizer(
    coords = line_coords(P, labels),
    lambda = lambda,
    dthresh = 1.01,
    nnk = 3,
    weight_mode = "binary",
    normalized = FALSE,
    handle_isolates = "keep_zero"
  )
}

test_that("coordinate construction dogfoods the adjoin Laplacian", {
  labels <- paste0("v", seq_len(8))
  coords <- line_coords(8, labels)
  spatial <- line_spatial(8, lambda = 0.5, labels = labels)
  direct <- adjoin::spatial_laplacian(
    coords,
    dthresh = 1.01,
    nnk = 3,
    weight_mode = "binary",
    sigma = 1.01 / 2,
    normalized = FALSE,
    stochastic = FALSE,
    handle_isolates = "keep_zero"
  )

  expect_s3_class(spatial, "dkge_spatial_regularizer")
  expect_equal(as.matrix(spatial$laplacians[[1]]), as.matrix(direct),
               tolerance = 0, ignore_attr = TRUE)
  expect_identical(spatial$construction$function_name,
                   "adjoin::spatial_laplacian")
  expect_identical(spatial$construction$adjoin_version,
                   as.character(packageVersion("adjoin")))
  expect_equal(as.numeric(spatial$laplacians[[1]] %*% rep(1, 8)),
               rep(0, 8), tolerance = 1e-14)
  expect_true(inherits(spatial$laplacians[[1]], "sparseMatrix"))
})

test_that("spatial constructor and domain contracts fail closed", {
  coords <- line_coords(5)
  L <- adjoin::spatial_laplacian(
    coords, dthresh = 1.01, nnk = 3, weight_mode = "binary",
    normalized = FALSE, handle_isolates = "keep_zero"
  )

  expect_error(dkge_spatial_regularizer(lambda = 1),
               class = "dkge_spatial_spec_error")
  expect_error(
    dkge_spatial_regularizer(coords = coords, laplacian = L, lambda = 1),
    class = "dkge_spatial_spec_error"
  )
  expect_error(dkge_spatial_regularizer(coords, lambda = -1),
               class = "dkge_spatial_spec_error")
  expect_error(
    dkge_spatial_regularizer(coords, lambda = 1, normalized = TRUE),
    class = "dkge_spatial_spec_error"
  )

  asymmetric <- as.matrix(L)
  asymmetric[1, 2] <- asymmetric[1, 2] / 2
  expect_error(
    dkge_spatial_regularizer(laplacian = asymmetric, lambda = 1),
    class = "dkge_spatial_laplacian_error"
  )

  bad_constant <- as.matrix(L)
  bad_constant[1, 1] <- bad_constant[1, 1] + 0.25
  expect_error(
    dkge_spatial_regularizer(laplacian = bad_constant, lambda = 1),
    class = "dkge_spatial_laplacian_error"
  )

  B <- replicate(3, matrix(rnorm(2 * 5), 2, 5), simplify = FALSE)
  colnames(B[[1]]) <- paste0("wrong", seq_len(5))
  X <- replicate(3, diag(2), simplify = FALSE)
  named_spatial <- line_spatial(5, lambda = 1,
                                labels = paste0("v", seq_len(5)))
  expect_error(
    dkge_fit(B, X, diag(2), rank = 1, w_method = "none",
             effect_scaling = "none", spatial = named_spatial),
    class = "dkge_spatial_domain_error"
  )

  ordered_B <- replicate(3, {
    out <- matrix(rnorm(2 * 5), 2, 5)
    colnames(out) <- rev(paste0("v", seq_len(5)))
    out
  }, simplify = FALSE)
  ordered_fit <- dkge_fit(
    ordered_B, X, diag(2), rank = 1, w_method = "none",
    effect_scaling = "none", spatial = named_spatial
  )
  expect_identical(
    ordered_fit$spatial$operators[[1]]$domain,
    colnames(ordered_B[[1]])
  )
  expect_identical(
    rownames(ordered_fit$spatial$operators[[1]]$L),
    colnames(ordered_B[[1]])
  )

  mixed_order <- ordered_B
  mixed_order[[2]] <- mixed_order[[2]][, rev(seq_len(5)), drop = FALSE]
  expect_error(
    dkge_fit(
      mixed_order, X, diag(2), rank = 1, w_method = "none",
      effect_scaling = "none", spatial = named_spatial
    ),
    class = "dkge_spatial_domain_error"
  )
})

test_that("regularized moment and maps match a dense resolvent oracle", {
  set.seed(7401)
  S <- 3L
  q <- 3L
  P <- 7L
  labels <- paste0("v", seq_len(P))
  B <- lapply(seq_len(S), function(s) {
    out <- matrix(rnorm(q * P), q, P,
                  dimnames = list(paste0("e", seq_len(q)), labels))
    out
  })
  X0 <- diag(q)
  colnames(X0) <- rownames(B[[1]])
  X <- replicate(S, X0, simplify = FALSE)
  K <- matrix(c(1.4, 0.2, 0.1,
                0.2, 1.1, 0.15,
                0.1, 0.15, 0.9), q, q,
              dimnames = list(rownames(B[[1]]), rownames(B[[1]])))
  omega <- lapply(seq_len(S), function(s) seq(0.6, 1.4, length.out = P))
  lambda <- 0.75
  spatial <- line_spatial(P, lambda, labels)

  fit <- dkge_fit(
    B, X, K,
    Omega_list = omega,
    rank = 2,
    w_method = "none",
    effect_scaling = "none",
    keep_X = TRUE,
    spatial = spatial
  )

  L <- as.matrix(spatial$laplacians[[1]])
  H <- solve(diag(P) + lambda * L)
  Bmodel <- lapply(B, function(Bs) Bs %*% t(H))
  moments <- Map(function(Bs, om) {
    W <- sweep(Bs, 2L, sqrt(om), "*")
    tcrossprod(W)
  }, Bmodel, omega)
  expected_chat <- fit$Khalf %*% Reduce(`+`, moments) %*% fit$Khalf
  expected_blocks <- Map(function(Bs, om) {
    fit$Khalf %*% sweep(Bs, 2L, sqrt(om), "*")
  }, Bmodel, omega)
  expected_X <- do.call(cbind, expected_blocks)

  expect_equal(fit$effect_moments, moments, tolerance = 1e-11,
               ignore_attr = TRUE)
  expect_equal(fit$Chat, expected_chat, tolerance = 1e-10,
               ignore_attr = TRUE)
  expect_equal(fit$X_concat, expected_X, tolerance = 1e-10,
               ignore_attr = TRUE)
  expect_equal(fit$Chat, tcrossprod(fit$X_concat), tolerance = 1e-10,
               ignore_attr = TRUE)

  projected <- dkge_project_btil(fit, fit$Btil)
  expected_projected <- lapply(Bmodel, function(Bs) {
    t(Bs) %*% fit$K %*% fit$U
  })
  expect_equal(projected, expected_projected, tolerance = 1e-10,
               ignore_attr = TRUE)
  expect_error(
    dkge_project_btil(
      fit, fit$Btil[[1]][, rev(seq_len(P)), drop = FALSE], subject = 1
    ),
    class = "dkge_spatial_domain_error"
  )

  projected_clusters <- dkge_project_clusters(
    fit, B[[1]], omega_vec = omega[[1]], subject = 1
  )
  expected_clusters <- t(fit$Khalf %*% Bmodel[[1]] %*%
                           diag(sqrt(omega[[1]]), P)) %*%
    sweep(fit$Khalf %*% fit$U, 2, fit$sdev, "/")
  expect_equal(projected_clusters, expected_clusters, tolerance = 1e-10,
               ignore_attr = TRUE)

  op <- fit$spatial$operators[[1]]
  constant_fields <- matrix(c(2, -1, 0.5), q, P)
  expect_equal(
    .dkge_spatial_apply_betas(constant_fields, op), constant_fields,
    tolerance = 1e-12
  )
  expect_null(op$H)
  expect_true(inherits(op$L, "sparseMatrix"))
  expect_true(inherits(op$A, "sparseMatrix"))
  expect_s4_class(op$factor, "CHMfactor")
  expect_identical(fit$provenance$spatial_regularization$method,
                   "laplacian_resolvent")
})

test_that("lambda zero is an exact identity fit", {
  set.seed(7402)
  S <- 4L
  q <- 3L
  P <- 8L
  B <- replicate(S, matrix(rnorm(q * P), q, P), simplify = FALSE)
  X <- replicate(S, diag(q), simplify = FALSE)

  baseline <- dkge_fit(B, X, diag(q), rank = 2, w_method = "none",
                       effect_scaling = "none", keep_X = TRUE)
  identity_fit <- dkge_fit(
    B, X, diag(q), rank = 2, w_method = "none",
    effect_scaling = "none", keep_X = TRUE,
    spatial = line_spatial(P, lambda = 0)
  )

  expect_equal(identity_fit$Chat, baseline$Chat, tolerance = 0)
  expect_equal(identity_fit$U, baseline$U, tolerance = 0)
  expect_equal(identity_fit$X_concat, baseline$X_concat, tolerance = 0)
  expect_false(identity_fit$spatial$active)

  cvec <- c(1, -1, 0)
  c0 <- dkge_contrast(baseline, cvec, method = "loso")
  c1 <- dkge_contrast(identity_fit, cvec, method = "loso")
  expect_equal(c1$values, c0$values, tolerance = 0)
})

test_that("spatial regularization changes the learned solution and reduces roughness", {
  P <- 30L
  high_frequency <- 2 * (-1)^seq_len(P)
  constant <- rep(1, P)
  B0 <- rbind(high_frequency = high_frequency, constant = constant)
  B <- replicate(4, B0, simplify = FALSE)
  X <- replicate(4, diag(2), simplify = FALSE)

  raw <- dkge_fit(B, X, diag(2), rank = 1, w_method = "none",
                  effect_scaling = "none")
  smooth <- dkge_fit(
    B, X, diag(2), rank = 1, w_method = "none",
    effect_scaling = "none",
    spatial = line_spatial(P, lambda = 5)
  )

  expect_gt(abs(raw$U[1, 1]), 0.99)
  expect_gt(abs(smooth$U[2, 1]), 0.99)
  expect_lt(abs(smooth$U[1, 1]), 1e-8)

  raw_map <- dkge_project_btil(raw, raw$Btil[[1]])[, 1]
  smooth_map <- dkge_project_btil(smooth, smooth$Btil[[1]])[, 1]
  expect_lt(sum(diff(smooth_map)^2), sum(diff(raw_map)^2))
  expect_true(all(is.finite(smooth_map)))
})

test_that("LOSO maps use the same fitted spatial operator", {
  set.seed(7403)
  S <- 4L
  q <- 2L
  P <- 9L
  B <- replicate(S, matrix(rnorm(q * P), q, P), simplify = FALSE)
  X <- replicate(S, diag(q), simplify = FALSE)
  spatial <- line_spatial(P, lambda = 0.6)
  fit <- dkge_fit(B, X, diag(q), rank = 1, w_method = "none",
                  effect_scaling = "none", spatial = spatial)

  held <- dkge_loso_contrast(fit, 1, c(1, 0))
  H <- solve(diag(P) + 0.6 * as.matrix(spatial$laplacians[[1]]))
  Bmodel <- B[[1]] %*% t(H)
  alpha <- as.numeric(t(held$basis) %*% c(1, 0))
  expected <- as.numeric(t(Bmodel) %*% held$basis %*% alpha)

  expect_equal(held$v, expected, tolerance = 1e-10)
})

test_that("prediction reuses shared domains and requires new subject-specific domains", {
  set.seed(7404)
  q <- 2L
  X <- list(s1 = diag(q), s2 = diag(q), s3 = diag(q))
  B <- lapply(c(s1 = 6L, s2 = 6L, s3 = 6L), function(P) {
    matrix(rnorm(q * P), q, P)
  })
  shared_fit <- dkge_fit(
    B, X, diag(q), rank = 1, w_method = "none",
    effect_scaling = "none", spatial = line_spatial(6, lambda = 0.4)
  )
  predicted <- dkge_predict_loadings(shared_fit, B)
  expected <- dkge_project_btil(shared_fit, shared_fit$Btil)
  expect_identical(names(expected), shared_fit$subject_ids)
  expect_equal(predicted, expected, tolerance = 1e-10, ignore_attr = TRUE)
  frozen <- dkge_freeze(shared_fit)
  expect_null(frozen$spatial$operators)
  expect_equal(
    dkge_predict_loadings(frozen, B), expected,
    tolerance = 1e-10, ignore_attr = TRUE
  )

  widths <- c(s1 = 5L, s2 = 6L, s3 = 7L)
  B_subject <- lapply(widths, function(P) matrix(rnorm(q * P), q, P))
  subject_spatial <- dkge_spatial_regularizer(
    coords = lapply(rev(widths), line_coords),
    lambda = 0.4,
    dthresh = 1.01,
    nnk = 3,
    weight_mode = "binary"
  )
  subject_data <- dkge_data(B_subject, X, subject_ids = names(B_subject))
  fit_subject <- dkge_fit(
    subject_data, K = diag(q), rank = 1, w_method = "none",
    effect_scaling = "none", spatial = subject_spatial
  )
  expect_error(
    dkge_project_btil(fit_subject, fit_subject$Btil[[1]]),
    class = "dkge_spatial_domain_error"
  )
  expect_error(
    dkge_predict_loadings(fit_subject, B_subject),
    class = "dkge_spatial_prediction_error"
  )
  expect_length(
    dkge_predict_loadings(fit_subject, B_subject,
                          spatial = subject_spatial),
    3L
  )
})

test_that("analytic noise debiasing is rejected but split-half moments are regularized", {
  set.seed(7405)
  q <- 2L
  P <- 6L
  B <- replicate(3, matrix(rnorm(q * P), q, P), simplify = FALSE)
  X <- replicate(3, diag(q), simplify = FALSE)
  spatial <- line_spatial(P, lambda = 0.5)
  expect_error(
    dkge_fit(B, X, diag(q), rank = 1, w_method = "none",
             effect_scaling = "none", debias = "analytic",
             spatial = spatial),
    class = "dkge_spatial_debias_error"
  )

  subjects <- lapply(seq_len(3), function(s) {
    split <- list(B[[s]] + 0.1, B[[s]] - 0.1)
    dkge_subject(B[[s]], X[[s]], id = paste0("s", s), split_betas = split)
  })
  split_fit <- dkge_fit(
    dkge_data(subjects), K = diag(q), rank = 1,
    w_method = "none", effect_scaling = "none",
    debias = "split_half", spatial = spatial
  )
  expect_true(all(is.finite(split_fit$Chat)))
  expect_true(split_fit$spatial$active)
  H <- solve(diag(P) + 0.5 * as.matrix(spatial$laplacians[[1]]))
  split1 <- subjects[[1]]$split_betas[[1]] %*% t(H)
  split2 <- subjects[[1]]$split_betas[[2]] %*% t(H)
  expected <- 0.5 * (split1 %*% t(split2) + split2 %*% t(split1))
  expect_equal(split_fit$effect_moments[[1]], expected, tolerance = 1e-10,
               ignore_attr = TRUE)
})

test_that("kernel-factor inputs cannot masquerade as a spatial domain", {
  spatial <- line_spatial(2, lambda = 0.5)
  expect_error(
    dkge_fit_from_kernels(
      list(s1 = diag(2), s2 = diag(2)),
      effect_ids = c("a", "b"),
      spatial = spatial
    ),
    class = "dkge_spatial_domain_error"
  )
})

test_that("spatial CV scores raw held-out blocks and prefers smoother one-SE fits", {
  set.seed(7406)
  S <- 4L
  P <- 8L
  B <- replicate(S, matrix(rnorm(P), 1, P), simplify = FALSE)
  X <- replicate(S, matrix(1, 1, 1), simplify = FALSE)
  template <- line_spatial(P, lambda = 1)

  expect_warning(
    cv <- dkge_cv_spatial_grid(
      B, X, matrix(1, 1, 1), template,
      lambdas = c(0, 0.25, 1), rank = 1,
      w_method = "none", effect_scaling = "none"
    ),
    class = "dkge_cv_saturation_warning"
  )
  expect_equal(cv$pick, 1)
  expect_equal(cv$best, 0)
  expect_identical(cv$heldout_geometry, "raw_beta_block")
  expect_identical(cv$selection_rule, "largest_lambda_within_one_se")
  expect_identical(cv$fit_settings$effect_scaling, "none")
  expect_equal(cv$spatial$lambda, 1)
  expect_true(all(cv$table$admissible))
})

test_that("spatial CV keeps Omega but omits candidate smoothing from held-out scores", {
  set.seed(7407)
  S <- 4L
  q <- 2L
  P <- 7L
  B <- replicate(S, matrix(rnorm(q * P), q, P), simplify = FALSE)
  X <- replicate(S, diag(q), simplify = FALSE)
  Omega <- lapply(seq_len(S), function(s) seq(0.5, 1.5, length.out = P))
  template <- line_spatial(P, lambda = 1)

  cv <- dkge_cv_spatial_grid(
    B, X, diag(q), template,
    lambdas = c(0, 0.5), rank = 1,
    Omega_list = Omega,
    w_method = "none", effect_scaling = "none"
  )

  row <- which(cv$raw$lambda == 0.5 & cv$raw$subject == 1)[[1]]
  lambda <- cv$raw$lambda[[row]]
  candidate <- .dkge_spatial_with_lambda(template, lambda)
  base <- dkge_fit(
    B, X, diag(q), Omega_list = Omega,
    rank = 1, w_method = "none", effect_scaling = "none",
    spatial = candidate
  )
  ctx <- .dkge_fold_weight_context(base, 2:S)
  fold <- .dkge_cv_fold_basis(ctx$Chat, base, rank = 1)
  heldout_fixed <- sweep(base$Btil[[1]], 2L, sqrt(Omega[[1]]), "*")
  expected <- .dkge_cv_score_fixed(
    heldout_fixed, fold$basis,
    .dkge_cv_validation_geometry(NULL, q)$roots
  )
  heldout_smoothed <- .dkge_apply_fit_spatial(
    base, base$Btil[[1]], subject = 1
  )
  heldout_smoothed <- sweep(
    heldout_smoothed, 2L, sqrt(Omega[[1]]), "*"
  )
  self_graded <- .dkge_cv_score_fixed(
    heldout_smoothed, fold$basis,
    .dkge_cv_validation_geometry(NULL, q)$roots
  )

  expect_equal(cv$raw$score[[row]], expected, tolerance = 1e-12)
  expect_false(isTRUE(all.equal(expected, self_graded, tolerance = 1e-8)))
  expect_identical(cv$heldout_spatial_metric, "Omega_list")
})
