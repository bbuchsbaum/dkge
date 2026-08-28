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

spatial_cycle_laplacian <- function(P) {
  stopifnot(P >= 3L)
  i <- seq_len(P)
  j <- c(i[-1L], 1L)
  adjacency <- Matrix::sparseMatrix(
    i = c(i, j), j = c(j, i), x = 1,
    dims = c(P, P)
  )
  Matrix::Diagonal(x = Matrix::rowSums(adjacency)) - adjacency
}

spatial_path_laplacian <- function(P) {
  stopifnot(P >= 2L)
  i <- seq_len(P - 1L)
  adjacency <- Matrix::sparseMatrix(
    i = c(i, i + 1L), j = c(i + 1L, i), x = 1,
    dims = c(P, P)
  )
  Matrix::Diagonal(x = Matrix::rowSums(adjacency)) - adjacency
}

irregular_3d_coords <- function() {
  rbind(
    a = c(x = 0, y = 0, z = 0),
    b = c(x = 0.15, y = 0.10, z = 0.75),
    c = c(x = 0.75, y = 0.15, z = 0.10),
    d = c(x = 0.65, y = 0.75, z = 0.25),
    e = c(x = 0.10, y = 0.85, z = 0.35),
    f = c(x = 0.75, y = 0.85, z = 0.90),
    far_z = c(x = 0, y = 0, z = 2.20),
    isolate = c(x = 5, y = 5, z = 5)
  )
}

irregular_3d_spatial <- function(lambda, coords = irregular_3d_coords()) {
  dkge_spatial_regularizer(
    coords = coords,
    lambda = lambda,
    dthresh = 1.05,
    nnk = nrow(coords),
    weight_mode = "heat",
    sigma = 0.5,
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

test_that("mutable spatial specifications are revalidated at fit and CV boundaries", {
  set.seed(7415)
  P <- 6L
  q <- 2L
  labels <- paste0("parcel", seq_len(P))
  B <- replicate(3L, {
    value <- matrix(rnorm(q * P), q, P)
    colnames(value) <- labels
    value
  }, simplify = FALSE)
  X <- replicate(3L, diag(q), simplify = FALSE)
  base <- line_spatial(P, lambda = 0.4, labels = labels)

  mutations <- list(
    diagonal_drift = list(
      class = "dkge_spatial_laplacian_error",
      apply = function(x) {
        x$laplacians[[1L]] <-
          x$laplacians[[1L]] + Matrix::Diagonal(P, x = rep(0.1, P))
        x
      }
    ),
    asymmetry = list(
      class = "dkge_spatial_laplacian_error",
      apply = function(x) {
        L <- Matrix::Matrix(as.matrix(x$laplacians[[1L]]), sparse = TRUE)
        L[1L, 2L] <- L[1L, 2L] + 0.1
        x$laplacians[[1L]] <- L
        x
      }
    ),
    domain = list(
      class = "dkge_spatial_domain_error",
      apply = function(x) {
        x$domains[[1L]][[1L]] <- x$domains[[1L]][[2L]]
        x
      }
    ),
    domain_removed = list(
      class = "dkge_spatial_domain_error",
      apply = function(x) {
        x$domains[1L] <- list(NULL)
        x
      }
    ),
    foreign_entry = list(
      class = "dkge_spatial_provenance_error",
      apply = function(x) { x$foreign_entry <- TRUE; x }
    ),
    construction = list(
      class = "dkge_spatial_provenance_error",
      apply = function(x) {
        x$construction$function_name <- "foreign::laplacian"
        x
      }
    )
  )

  for (mutation in mutations) {
    candidate <- mutation$apply(base)
    expect_error(
      dkge_fit(
        B, X, diag(q), rank = 1L, w_method = "none",
        effect_scaling = "none", spatial = candidate
      ),
      class = mutation$class
    )
    expect_error(
      dkge_cv_spatial_grid(
        B, X, diag(q), candidate,
        lambdas = c(0, 0.4), rank = 1L,
        w_method = "none", effect_scaling = "none"
      ),
      class = mutation$class
    )
    expect_error(
      .dkge_spatial_with_lambda(candidate, 0.2),
      class = mutation$class
    )
  }

  invalid_display <- mutations$diagonal_drift$apply(base)
  expect_error(
    print(invalid_display),
    class = "dkge_spatial_laplacian_error"
  )

  fitted <- dkge_fit(
    B, X, diag(q), rank = 1L, w_method = "none",
    effect_scaling = "none", spatial = base
  )
  status_spoof <- fitted
  status_spoof$spatial$status <- "inert"
  for (audit in list(print, dkge_diagnostics)) {
    expect_error(
      audit(status_spoof),
      class = "dkge_spatial_provenance_error"
    )
  }
  diagonal_spoof <- fitted
  diagonal_spoof$spatial$spec$laplacians[[1L]] <-
    diagonal_spoof$spatial$spec$laplacians[[1L]] +
    Matrix::Diagonal(P, x = rep(0.1, P))
  for (audit in list(print, dkge_diagnostics)) {
    expect_error(
      audit(diagonal_spoof),
      class = "dkge_spatial_laplacian_error"
    )
  }
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

test_that("the resolvent has the analytic Laplacian eigenmode response", {
  P <- 12L
  lambda <- 0.7
  L <- spatial_cycle_laplacian(P)
  eig <- eigen(as.matrix(L), symmetric = TRUE)
  modes <- t(eig$vectors)
  spatial <- dkge_spatial_regularizer(laplacian = L, lambda = lambda)
  resolved <- .dkge_resolve_spatial(
    spatial, list(modes), subject_ids = "s1"
  )

  observed <- .dkge_spatial_apply_betas(
    modes, resolved$operators[[1]]
  )
  gains <- 1 / (1 + lambda * eig$values)
  expected <- sweep(modes, 1L, gains, "*")
  observed_gains <- sqrt(rowSums(observed^2) / rowSums(modes^2))

  expect_equal(observed, expected, tolerance = 1e-10)
  expect_equal(observed_gains, gains, tolerance = 1e-10)
  expect_true(all(diff(observed_gains) >= -1e-12))

  constant_mode <- which.min(abs(eig$values))
  expect_equal(
    observed[constant_mode, ], modes[constant_mode, ],
    tolerance = 1e-12
  )
})

test_that("irregular three-dimensional geometry smooths within graph components", {
  coords <- irregular_3d_coords()
  spatial <- irregular_3d_spatial(lambda = 1.5, coords = coords)
  L <- spatial$laplacians[[1]]

  # These checks require all three coordinates: b is a genuine 3-D neighbour
  # of a, while far_z shares a's x/y location but lies beyond the radius.
  expect_lt(L["a", "b"], 0)
  expect_equal(L["a", "far_z"], 0)
  expect_equal(
    Matrix::rowSums(abs(L[c("far_z", "isolate"), , drop = FALSE])),
    c(far_z = 0, isolate = 0)
  )

  raw <- matrix(
    c(2, -1, 3, -2, 1, 0.5, 7, 9),
    nrow = 1L,
    dimnames = list("effect", rownames(coords))
  )
  resolved <- .dkge_resolve_spatial(
    spatial, list(raw), subject_ids = "s1"
  )
  smoothed <- .dkge_spatial_apply_betas(raw, resolved$operators[[1]])
  raw_energy <- as.numeric(raw %*% L %*% t(raw))
  smoothed_energy <- as.numeric(smoothed %*% L %*% t(smoothed))
  connected <- c("a", "b", "c", "d", "e", "f")
  isolates <- c("far_z", "isolate")

  expect_gt(raw_energy, 0)
  expect_lt(smoothed_energy, raw_energy)
  expect_equal(
    sum(smoothed[, connected]), sum(raw[, connected]),
    tolerance = 1e-12
  )
  expect_equal(
    smoothed[, isolates], raw[, isolates],
    tolerance = 1e-12
  )
  expect_identical(spatial$construction$weight_mode, "heat")
})

test_that("fitting and projection are equivariant to spatial permutation", {
  set.seed(7409)
  coords <- irregular_3d_coords()
  labels <- rownames(coords)
  S <- 4L
  q <- 3L
  B <- replicate(S, {
    out <- matrix(rnorm(q * length(labels)), q, length(labels))
    colnames(out) <- labels
    out
  }, simplify = FALSE)
  X <- replicate(S, diag(q), simplify = FALSE)
  K <- matrix(c(1.3, 0.2, 0.1,
                0.2, 1.1, 0.15,
                0.1, 0.15, 0.9), q, q)
  spatial <- irregular_3d_spatial(lambda = 0.8, coords = coords)
  reference <- dkge_fit(
    B, X, K, rank = 2, w_method = "none",
    effect_scaling = "none", spatial = spatial
  )

  permutation <- c(4L, 1L, 8L, 2L, 6L, 3L, 7L, 5L)
  permuted_B <- lapply(B, function(block) {
    block[, permutation, drop = FALSE]
  })
  permuted_spatial <- irregular_3d_spatial(
    lambda = 0.8,
    coords = coords[permutation, , drop = FALSE]
  )
  permuted <- dkge_fit(
    permuted_B, X, K, rank = 2, w_method = "none",
    effect_scaling = "none", spatial = permuted_spatial
  )
  inverse <- order(permutation)

  expect_equal(
    as.matrix(permuted_spatial$laplacians[[1]][inverse, inverse]),
    as.matrix(spatial$laplacians[[1]]),
    tolerance = 1e-14
  )
  expect_equal(permuted$Chat, reference$Chat, tolerance = 1e-11)
  expect_equal(permuted$sdev, reference$sdev, tolerance = 1e-11)
  expect_equal(
    tcrossprod(permuted$U), tcrossprod(reference$U),
    tolerance = 1e-10
  )

  reference_maps <- dkge_project_btil(reference, reference$Btil)
  permuted_maps <- dkge_project_btil(permuted, permuted$Btil)
  first_restored <- permuted_maps[[1]][inverse, , drop = FALSE]
  signs <- sign(colSums(reference_maps[[1]] * first_restored))
  signs[signs == 0] <- 1
  for (s in seq_len(S)) {
    restored <- permuted_maps[[s]][inverse, , drop = FALSE]
    restored <- sweep(restored, 2L, signs, "*")
    expect_equal(restored, reference_maps[[s]], tolerance = 1e-10)
  }
})

test_that("spatial CV recovers a shared smooth truth over subject-specific rough noise", {
  set.seed(7410)
  S <- 6L
  q <- S + 1L
  P <- 24L
  position <- seq_len(P)
  smooth_truth <- sqrt(2) * sin(2 * pi * position / P)
  rough_noise <- (-1)^position
  B <- lapply(seq_len(S), function(s) {
    out <- 1e-6 * matrix(rnorm(q * P), q, P)
    out[1L, ] <- out[1L, ] + smooth_truth
    out[s + 1L, ] <- out[s + 1L, ] + 3 * rough_noise
    out
  })
  X <- replicate(S, diag(q), simplify = FALSE)
  L <- spatial_cycle_laplacian(P)
  spatial <- dkge_spatial_regularizer(laplacian = L, lambda = 1)

  smooth_energy <- as.numeric(t(smooth_truth) %*% L %*% smooth_truth)
  rough_energy <- as.numeric(t(rough_noise) %*% L %*% rough_noise)
  expect_lt(smooth_energy, rough_energy)

  cv <- dkge_cv_spatial_grid(
    B, X, diag(q), spatial,
    lambdas = c(0, 0.25, 1), rank = 1,
    w_method = "none", effect_scaling = "none"
  )
  zero_score <- cv$table$mean[cv$table$lambda == 0]
  regularized_scores <- cv$table$mean[cv$table$lambda > 0]

  expect_gt(min(regularized_scores), zero_score + 0.09)
  expect_equal(cv$pick, 1)

  raw_fit <- dkge_fit(
    B, X, diag(q), rank = 1, w_method = "none",
    effect_scaling = "none"
  )
  selected_fit <- dkge_fit(
    B, X, diag(q), rank = 1, w_method = "none",
    effect_scaling = "none", spatial = cv$spatial
  )
  expect_lt(abs(raw_fit$U[1L, 1L]), 1e-4)
  expect_gt(abs(selected_fit$U[1L, 1L]), 0.9999)
})

test_that("large spatial domains retain sparse linear storage", {
  P <- 10000L
  L <- spatial_path_laplacian(P)
  B <- rbind(
    rough = (-1)^seq_len(P),
    smooth = sin(seq_len(P) / 17)
  )
  spatial <- dkge_spatial_regularizer(laplacian = L, lambda = 0.75)
  resolved <- .dkge_resolve_spatial(
    spatial, list(B), subject_ids = "s1"
  )
  operator <- resolved$operators[[1]]
  smoothed <- .dkge_spatial_apply_betas(B, operator)

  expect_s4_class(operator$L, "sparseMatrix")
  expect_s4_class(operator$A, "sparseMatrix")
  expect_s4_class(operator$factor, "CHMfactor")
  expect_equal(Matrix::nnzero(operator$L), 3L * P - 2L)
  expect_equal(Matrix::nnzero(operator$A), 3L * P - 2L)
  expect_null(operator$H)

  dense_operator_bytes <- 8 * as.double(P)^2
  expect_lt(
    as.numeric(object.size(operator)),
    dense_operator_bytes / 100
  )
  expect_equal(dim(smoothed), dim(B))
  expect_true(all(is.finite(smoothed)))
  expect_lt(sum(diff(smoothed[1L, ])^2), sum(diff(B[1L, ])^2))
})

test_that("an edgeless graph is reported rather than silently inert", {
  # The defaults assume unit-spaced voxel indices. Millimetre coordinates with
  # the default dthresh isolate every unit, giving L = 0 and H = I: the fit is
  # then bit-identical to an unregularized one despite lambda > 0. That must be
  # visible at construction, not discovered later.
  mm <- cbind(x = seq(0, 33, by = 3), y = 0, z = 0)

  expect_warning(
    spatial <- dkge_spatial_regularizer(coords = mm, lambda = 2),
    "no edges",
    class = "dkge_spatial_inert_warning"
  )
  expect_equal(dkge:::.dkge_spatial_edge_count(spatial$laplacians[[1]]), 0)
  expect_identical(spatial$topology$status, "inert")
  expect_true(spatial$topology$requested)
  expect_false(spatial$topology$effective)
  expect_false(spatial$topology$fully_effective)
  expect_output(print(spatial), "edges\\s*: 0")
  expect_output(print(spatial), "status\\s*: inert")
  expect_output(print(spatial), "no effect")

  # A threshold on the right scale for 3 mm spacing connects the chain.
  expect_silent(
    connected <- dkge_spatial_regularizer(coords = mm, lambda = 2,
                                          dthresh = 3.1, nnk = 6)
  )
  expect_equal(dkge:::.dkge_spatial_edge_count(connected$laplacians[[1]]),
               nrow(mm) - 1)
  expect_identical(connected$topology$status, "active")
  expect_true(connected$topology$effective)
  expect_true(connected$topology$fully_effective)

  # lambda = 0 is already an explicit no-op, so it must not warn.
  expect_silent(inactive <- dkge_spatial_regularizer(coords = mm, lambda = 0))
  expect_identical(inactive$topology$status, "inactive")
  expect_false(inactive$topology$requested)
  expect_false(inactive$topology$effective)
})

test_that("edgeless status survives fitting and stale positive retunes re-warn", {
  set.seed(7412)
  mm <- cbind(x = seq(0, 21, by = 3), y = 0, z = 0)
  S <- 3L
  q <- 2L
  P <- nrow(mm)
  B <- replicate(S, matrix(rnorm(q * P), q, P), simplify = FALSE)
  X <- replicate(S, diag(q), simplify = FALSE)
  raw <- dkge_fit(B, X, diag(q), rank = 1, w_method = "none",
                  effect_scaling = "none")

  expect_warning(
    spatial <- dkge_spatial_regularizer(mm, lambda = 2),
    class = "dkge_spatial_inert_warning"
  )
  expect_silent(
    fit <- dkge_fit(B, X, diag(q), rank = 1, w_method = "none",
                    effect_scaling = "none", spatial = spatial)
  )
  diagnostics <- dkge_diagnostics(fit)$spatial
  expect_equal(fit$Chat, raw$Chat, tolerance = 0)
  expect_identical(diagnostics$status, "inert")
  expect_true(diagnostics$requested)
  expect_false(diagnostics$active)
  expect_false(diagnostics$effective)
  expect_false(diagnostics$fully_effective)
  expect_false(any(diagnostics$diagnostics$effective))
  expect_output(print(fit), "Spatial smoothing: inert", fixed = TRUE)

  expect_silent(template <- dkge_spatial_regularizer(mm, lambda = 0))
  retuned <- .dkge_spatial_with_lambda(template, 2)
  expect_identical(retuned$topology$status, "inert")
  expect_warning(
    dkge_fit(B, X, diag(q), rank = 1, w_method = "none",
             effect_scaling = "none", spatial = retuned),
    class = "dkge_spatial_inert_warning"
  )
})

test_that("spatial CV fails closed when positive lambdas cannot change the fit", {
  set.seed(7413)
  mm <- cbind(x = seq(0, 21, by = 3), y = 0, z = 0)
  B <- replicate(3L, matrix(rnorm(2L * nrow(mm)), 2L, nrow(mm)),
                 simplify = FALSE)
  X <- replicate(3L, diag(2L), simplify = FALSE)
  expect_silent(template <- dkge_spatial_regularizer(mm, lambda = 0))

  expect_error(
    dkge_cv_spatial_grid(
      B, X, diag(2L), template,
      lambdas = c(0, 1, 2), rank = 1L,
      w_method = "none", effect_scaling = "none"
    ),
    "cannot distinguish",
    class = "dkge_cv_spatial_inert_error"
  )

  # An all-zero grid is an explicit request for the identity operator, not an
  # attempt to tune an unidentifiable positive penalty.
  expect_no_error(
    cv_zero <- suppressWarnings(dkge_cv_spatial_grid(
      B, X, diag(2L), template,
      lambdas = 0, rank = 1L,
      w_method = "none", effect_scaling = "none"
    ))
  )
  expect_identical(cv_zero$pick, 0)
  expect_identical(cv_zero$spatial$topology$status, "inactive")
})

test_that("subject-specific edgeless graphs are reported as partially active", {
  set.seed(7414)
  P <- 6L
  coords <- list(
    sub01 = line_coords(P),
    sub02 = cbind(x = seq(0, by = 3, length.out = P), y = 0, z = 0)
  )
  expect_warning(
    spatial <- dkge_spatial_regularizer(coords, lambda = 1,
                                        dthresh = 1.42, nnk = 6),
    "sub02",
    class = "dkge_spatial_partial_warning"
  )
  expect_identical(spatial$topology$status, "partial")
  expect_identical(spatial$topology$empty_labels, "sub02")

  B <- lapply(coords, function(x) matrix(rnorm(2L * nrow(x)), 2L, nrow(x)))
  X <- lapply(coords, function(x) diag(2L))
  data <- dkge_data(B, X, subject_ids = names(B))
  expect_silent(
    fit <- dkge_fit(data, K = diag(2L), rank = 1L,
                    w_method = "none", effect_scaling = "none",
                    spatial = spatial)
  )
  diagnostics <- dkge_diagnostics(fit)$spatial
  expect_identical(diagnostics$status, "partial")
  expect_true(diagnostics$active)
  expect_true(diagnostics$effective)
  expect_false(diagnostics$fully_effective)
  expect_identical(diagnostics$diagnostics$effective, c(TRUE, FALSE))

  observed_empty <- dkge:::.dkge_spatial_apply_betas(
    B$sub02, fit$spatial$operators$sub02
  )
  expect_identical(observed_empty, B$sub02)

  warning_classes <- list()
  cv <- withCallingHandlers(
    dkge_cv_spatial_grid(
      B, X, diag(2L), spatial,
      lambdas = c(0, 0.25), rank = 1L,
      w_method = "none", effect_scaling = "none"
    ),
    warning = function(w) {
      warning_classes[[length(warning_classes) + 1L]] <<- class(w)
      invokeRestart("muffleWarning")
    }
  )
  partial_warning_n <- sum(vapply(
    warning_classes,
    function(classes) "dkge_spatial_partial_warning" %in% classes,
    logical(1)
  ))
  expect_equal(partial_warning_n, 1L)
  expect_identical(cv$spatial$topology$status, "partial")
})

test_that("voxel weights smooth inside the resolvent and Omega outside it", {
  # Pins the documented moment algebra, including the asymmetry between the two
  # per-unit weightings:
  #   M = (Btil W^1/2) H Omega H (Btil W^1/2)'
  # Voxel weights enter before the smoother, Omega after it.
  set.seed(7411)
  P <- 8L
  q <- 3L
  B <- matrix(rnorm(q * P), q, P)
  spatial <- line_spatial(P, lambda = 1.5)
  fit <- dkge_fit(
    list(B, matrix(rnorm(q * P), q, P)),
    list(diag(q), diag(q)), diag(q),
    rank = 1, spatial = spatial,
    w_method = "none", effect_scaling = "none"
  )
  op <- fit$spatial$operators[[1]]
  H <- solve(diag(P) + op$lambda * as.matrix(op$L))
  w <- seq(0.3, 1.7, length.out = P)
  omega <- seq(0.5, 1.5, length.out = P)

  observed <- dkge:::.dkge_effect_moment(B, Omega = omega, voxel_weights = w,
                                         spatial = op)
  inner <- (B %*% diag(sqrt(w))) %*% H %*% diag(sqrt(omega))
  expect_equal(unname(observed), unname(inner %*% t(inner)), tolerance = 1e-10)

  # The alternative orderings must not coincide, or the assertion above is
  # vacuous and the documented distinction is untestable.
  swapped <- (B %*% H) %*% diag(sqrt(w * omega))
  expect_false(isTRUE(all.equal(unname(observed), unname(swapped %*% t(swapped)),
                                tolerance = 1e-10)))

  # H really does appear twice: the moment is not Btil H Btil'.
  plain <- dkge:::.dkge_effect_moment(B, spatial = op)
  expect_equal(unname(plain), unname(B %*% H %*% H %*% t(B)), tolerance = 1e-10)
  expect_false(isTRUE(all.equal(unname(plain), unname(B %*% H %*% t(B)),
                                tolerance = 1e-10)))
})
