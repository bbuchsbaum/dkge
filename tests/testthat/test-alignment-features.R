library(testthat)

make_alignment_feature_fixture <- function(seed = 260827L,
                                           S = 6L,
                                           q = 6L,
                                           P = 12L,
                                           kernel_rank = 5L,
                                           fit_rank = 2L) {
  set.seed(seed)
  effects <- paste0("e", seq_len(q))
  parcels <- paste0("p", seq_len(P))
  eigenvalues <- c(seq(kernel_rank, 1, length.out = kernel_rank),
                   rep(0, q - kernel_rank))
  K <- diag(eigenvalues, q)
  dimnames(K) <- list(effects, effects)
  covariance_seed <- matrix(rnorm(q * q), q, q)
  Lambda <- crossprod(covariance_seed) + diag(q)
  dimnames(Lambda) <- list(effects, effects)
  subjects <- lapply(seq_len(S), function(s) {
    B <- matrix(rnorm(q * P), q, P,
                dimnames = list(effects, parcels))
    X <- matrix(rnorm((q + 5L) * q), q + 5L, q,
                dimnames = list(NULL, effects))
    dkge_subject(
      B, X, id = paste0("s", s),
      effect_noise_cov = Lambda,
      residual_variance = stats::setNames(rep(1, P), parcels)
    )
  })
  fit <- suppressWarnings(dkge_fit(
    dkge_data(subjects), K = K, rank = fit_rank,
    w_method = "none", effect_scaling = "pooled_design"
  ))
  contrasts <- cbind(
    first = c(1, -1, rep(0, q - 2L)),
    second = c(0, 0, 1, -1, rep(0, q - 4L))
  )
  contrast_obj <- suppressWarnings(dkge_contrast(
    fit, contrasts, method = "loso", align = FALSE
  ))
  independent_betas <- lapply(fit$Braw, function(B) {
    B + matrix(rnorm(length(B), sd = 0.35), nrow(B), ncol(B),
               dimnames = dimnames(B))
  })
  names(independent_betas) <- fit$subject_ids
  centroids <- lapply(seq_len(S), function(s) {
    cbind(x = seq_len(P), y = s / 10, z = 0)
  })
  list(
    fit = fit,
    contrast = contrast_obj,
    independent_betas = independent_betas,
    centroids = centroids,
    K = K,
    Lambda = Lambda
  )
}

test_that("compact image(K) coordinates preserve full-root feature geometry", {
  set.seed(260827)
  q <- 8L
  P <- 17L
  Q <- qr.Q(qr(matrix(rnorm(q * q), q, q)))
  eigenvalues <- c(5, 2, 0.7, 0.2, rep(0, 4))
  K <- Q %*% diag(eigenvalues) %*% t(Q)
  Btil <- matrix(rnorm(q * P), q, P)
  compact <- dkge:::.dkge_compact_kernel_factor(K, 1e-10)
  roots <- dkge:::.dkge_kernel_geometry(K, tol = 1e-10)
  G <- t(Btil) %*% compact$factor
  Z <- t(Btil) %*% roots$Khalf

  expect_equal(ncol(G), 4L)
  expect_equal(compact$factor %*% t(compact$factor), K,
               tolerance = 1e-10)
  expect_equal(tcrossprod(G), tcrossprod(Z), tolerance = 1e-10)
  expect_equal(
    Z,
    G %*% t(compact$support_vectors),
    tolerance = 1e-10
  )
})

test_that("conditional residualization matches a hand-matrix oracle", {
  set.seed(260828)
  P <- 19L
  k <- 7L
  m <- 2L
  G <- matrix(rnorm(P * k), P, k)
  A <- matrix(rnorm(k * k), k, k)
  Sigma <- crossprod(A) + diag(k)
  gamma <- matrix(rnorm(k * m), k, m)
  actual <- dkge_residualize_alignment_features(G, gamma, Sigma)

  C <- crossprod(gamma, Sigma %*% gamma)
  coefficient <- solve(C, crossprod(gamma, Sigma))
  expected <- G - (G %*% gamma) %*% coefficient
  expected_covariance <- Sigma - Sigma %*% gamma %*% coefficient

  expect_equal(actual$features, expected, tolerance = 1e-11)
  expect_equal(actual$residual_covariance, expected_covariance,
               tolerance = 1e-11)
  expect_equal(actual$values, G %*% gamma, tolerance = 0)
  expect_lt(max(abs(actual$covariance_cross_oracle)), 1e-10)
  expect_match(
    actual$diagnostics$conditional_independence,
    "conditional on a fixed generating direction",
    fixed = TRUE
  )
})

test_that("the separable oracle covers every cross-parcel covariance", {
  set.seed(260829)
  P <- 5L
  k <- 6L
  m <- 2L
  G <- matrix(rnorm(80 * k), 80, k)
  effect_seed <- matrix(rnorm(k * k), k, k)
  Sigma <- crossprod(effect_seed) + diag(k)
  gamma <- matrix(rnorm(k * m), k, m)
  spatial_seed <- matrix(rnorm(P * P), P, P)
  rho <- stats::cov2cor(crossprod(spatial_seed) + diag(P))
  result <- dkge_residualize_alignment_features(G, gamma, Sigma)

  parcel_pair_covariances <- lapply(seq_len(P), function(p) {
    lapply(seq_len(P), function(q) {
      rho[p, q] * result$residual_covariance %*% gamma
    })
  })
  expect_lt(max(abs(unlist(parcel_pair_covariances))), 1e-9)
  expect_equal(result$diagnostics$contrast_span_rank, m)
  expect_equal(result$diagnostics$available_dimension, k - m)
})

test_that("residualization is equivariant to feature scale and permutation", {
  set.seed(260830)
  P <- 21L
  k <- 6L
  G <- matrix(rnorm(P * k), P, k)
  A <- matrix(rnorm(k * k), k, k)
  Sigma <- crossprod(A) + diag(k)
  gamma <- matrix(rnorm(k * 2L), k, 2L)
  base <- dkge_residualize_alignment_features(G, gamma, Sigma)

  permutation <- c(4, 1, 6, 2, 5, 3)
  transform <- diag(c(0.5, 2, 1.5, 0.75, 3, 1))[permutation, , drop = FALSE]
  transformed <- dkge_residualize_alignment_features(
    G %*% transform,
    solve(transform, gamma),
    t(transform) %*% Sigma %*% transform
  )
  expect_equal(transformed$values, base$values, tolerance = 1e-10)
  expect_equal(transformed$features, base$features %*% transform,
               tolerance = 1e-9)

  parcel_permutation <- sample(seq_len(P))
  permuted <- dkge_residualize_alignment_features(
    G[parcel_permutation, , drop = FALSE], gamma, Sigma
  )
  expect_equal(permuted$features,
               base$features[parcel_permutation, , drop = FALSE],
               tolerance = 1e-11)
})

test_that("production residual features retain rank(K) minus family span", {
  fx <- make_alignment_feature_fixture()
  features <- dkge_alignment_features(
    fx$fit, fx$contrast,
    feature_source = "same_data_residualized"
  )
  expect_s3_class(features, "dkge_alignment_features")
  expect_identical(features$feature_source, "same_data_residualized")
  expect_identical(features$eligibility$status, "approximate")
  expect_equal(features$kernel_rank, 5L)
  expect_equal(ncol(fx$fit$U), 2L)
  expect_true(all(vapply(features$diagnostics, function(x) {
    x$available_dimension == 3L &&
      x$residual_covariance_rank == 3L &&
      x$reconstruction_error < 1e-8 &&
      x$covariance_orthogonality_error < 1e-8
  }, logical(1))))
  expect_true(all(vapply(features$diagnostics, function(x) {
    x$full_square_root_equivalence_error < 1e-10
  }, logical(1))))
  L <- features$kernel_factor
  expected_covariance <- t(L) %*% t(fx$fit$R) %*% fx$Lambda %*%
    fx$fit$R %*% L
  expect_equal(features$feature_covariances[[1]], expected_covariance,
               tolerance = 1e-10)
  expect_match(
    features$provenance$inferential_contract,
    "only full null-action re-estimation is exact",
    fixed = TRUE
  )
})

test_that("independent is the default verified feature source", {
  fx <- make_alignment_feature_fixture(seed = 260831)
  features <- dkge_alignment_features(
    fx$fit, fx$contrast,
    independent_betas = fx$independent_betas,
    independent_data_hash = "independent-acquisition-260831"
  )
  expect_identical(features$feature_source, "independent")
  expect_true(features$source_verified)
  expect_true(all(vapply(features$features, ncol, integer(1)) == 5L))
  expect_true(all(vapply(features$diagnostics, function(x) {
    x$available_dimension == 5L &&
      x$full_square_root_equivalence_error < 1e-10
  }, logical(1))))
  expected_model <- t(fx$fit$R) %*% fx$independent_betas[[1]]
  expect_equal(
    features$features[[1]],
    t(expected_model) %*% features$kernel_factor,
    tolerance = 1e-11
  )
  expect_identical(features$eligibility$status, "approximate")

  reversed <- fx$independent_betas[rev(fx$fit$subject_ids)]
  reordered <- dkge_alignment_features(
    fx$fit, fx$contrast,
    independent_betas = reversed,
    independent_data_hash = "independent-acquisition-260831"
  )
  expect_identical(reordered$features, features$features)

  invalid_name_sets <- list(
    partial = replace(fx$fit$subject_ids, 2L, ""),
    duplicate = replace(
      fx$fit$subject_ids, length(fx$fit$subject_ids), fx$fit$subject_ids[[1L]]
    ),
    foreign = replace(
      fx$fit$subject_ids, length(fx$fit$subject_ids), "foreign-subject"
    )
  )
  for (bad_names in invalid_name_sets) {
    bad_betas <- fx$independent_betas
    names(bad_betas) <- bad_names
    expect_error(
      dkge_alignment_features(
        fx$fit, fx$contrast,
        independent_betas = bad_betas,
        independent_data_hash = "invalid-subject-names"
      ),
      "must match fitted subject IDs exactly",
      class = "dkge_alignment_feature_error"
    )
  }

  expect_error(
    dkge_alignment_features(
      fx$fit, fx$contrast,
      independent_betas = fx$fit$Braw,
      independent_data_hash = "false-independence"
    ),
    "identical to the fitted data",
    class = "dkge_alignment_feature_error"
  )
})

test_that("typed over-ranked features flow into transport without loose matrices", {
  fx <- make_alignment_feature_fixture(seed = 260832, P = 8L)
  features <- dkge_alignment_features(
    fx$fit, fx$contrast,
    independent_betas = fx$independent_betas,
    independent_data_hash = "independent-acquisition-260832"
  )
  transported <- dkge_transport_contrasts_to_reference(
    fx$fit, fx$contrast,
    reference_subject = 2L,
    selection_method = "explicit",
    centroids = fx$centroids,
    alignment_features = features,
    mapper = dkge_mapper_spec("ridge", lambda = 1e-3)
  )
  fitted <- attr(transported, "fitted_alignment")
  expect_identical(attr(transported, "alignment_mode"), "independent")
  expect_identical(attr(transported, "alignment_features_hash"),
                   features$structural_hash)
  expect_identical(fitted$feature_source, "independent")
  expect_identical(fitted$eligibility$status, "approximate")

  expect_error(
    suppressWarnings(dkge_transport_contrasts_to_medoid(
      fx$fit, fx$contrast, medoid = 2L,
      centroids = fx$centroids,
      alignment_features = features,
      loadings = features$features,
      mapper = dkge_mapper_spec("ridge", lambda = 1e-3)
    )),
    class = "dkge_alignment_feature_error"
  )
  expect_error(
    suppressWarnings(dkge_transport_contrasts_to_medoid(
      fx$fit, fx$contrast, medoid = 2L,
      centroids = fx$centroids,
      alignment_mode = "independent",
      mapper = dkge_mapper_spec("ridge", lambda = 1e-3)
    )),
    "requires a typed",
    class = "dkge_alignment_feature_error"
  )
})

test_that("feature objects fail closed under mutation and mixed contrast provenance", {
  fx <- make_alignment_feature_fixture(seed = 260833)
  features <- dkge_alignment_features(
    fx$fit, fx$contrast,
    independent_betas = fx$independent_betas,
    independent_data_hash = "independent-acquisition-260833"
  )
  mutated <- features
  mutated$features[[1]][1, 1] <- mutated$features[[1]][1, 1] + 0.1
  expect_error(
    print(mutated),
    "mutated",
    class = "dkge_alignment_feature_error"
  )
  other_contrast <- fx$contrast
  other_contrast$values[[1]][[1]][1] <-
    other_contrast$values[[1]][[1]][1] + 0.1
  expect_error(
    suppressWarnings(dkge_transport_contrasts_to_medoid(
      fx$fit, other_contrast, medoid = 1L,
      centroids = fx$centroids,
      alignment_features = features,
      mapper = dkge_mapper_spec("ridge", lambda = 1e-3)
    )),
    "different contrast result",
    class = "dkge_alignment_feature_error"
  )
})

test_that("feature provenance binds the exact estimator and rejects CPCA or JD", {
  fx <- make_alignment_feature_fixture(seed = 260856)
  features <- dkge_alignment_features(
    fx$fit, fx$contrast,
    independent_betas = fx$independent_betas,
    independent_data_hash = "independent-acquisition-260856"
  )
  base_binding <- dkge:::.dkge_alignment_fit_binding(fx$fit)

  rank_changed <- fx$fit
  rank_changed$rank <- (rank_changed$rank %||% ncol(rank_changed$U)) + 1L
  weight_changed <- fx$fit
  weight_changed$weights[[1L]] <- weight_changed$weights[[1L]] + 0.01
  basis_changed <- fx$fit
  basis_changed$U[1L, 1L] <- basis_changed$U[1L, 1L] + 1e-6
  kernel_root_changed <- fx$fit
  kernel_root_changed$Kihalf[1L, 2L] <-
    kernel_root_changed$Kihalf[1L, 2L] + 0.4
  support_changed <- fx$fit
  support_changed$kernel_support_projector[1L, 1L] <-
    support_changed$kernel_support_projector[1L, 1L] - 0.1
  pool_cache_changed <- fx$fit
  pool_cache_changed$pool_cache$structural[[1L]][1L, 1L] <-
    pool_cache_changed$pool_cache$structural[[1L]][1L, 1L] + 0.25
  provenance_changed <- fx$fit
  provenance_changed$provenance$alignment_adversary <- TRUE
  for (changed in list(
      rank_changed, weight_changed, basis_changed, kernel_root_changed,
      support_changed, pool_cache_changed, provenance_changed
  )) {
    expect_false(identical(
      dkge:::.dkge_alignment_fit_binding(changed), base_binding
    ))
    expect_error(
      dkge:::.dkge_validate_alignment_features(features, fit = changed),
      "mixed provenance",
      class = "dkge_alignment_feature_error"
    )
  }

  jd_changed <- fx$fit
  jd_changed$solver <- "jd"
  cpca_changed <- fx$fit
  cpca_changed$cpca <- list(part = "design")
  for (changed in list(jd_changed, cpca_changed)) {
    expect_false(identical(
      dkge:::.dkge_alignment_fit_binding(changed), base_binding
    ))
    expect_error(
      dkge:::.dkge_validate_alignment_features(features, fit = changed),
      "ordinary pooled eigensolve",
      class = "dkge_crossfit_estimator_error"
    )
    expect_error(
      dkge_alignment_features(
        changed, fx$contrast,
        independent_betas = fx$independent_betas,
        independent_data_hash = "independent-acquisition-260856"
      ),
      "ordinary pooled eigensolve",
      class = "dkge_crossfit_estimator_error"
    )
  }
})

test_that("rank, condition, energy, and covariance gates fail closed", {
  set.seed(260834)
  G <- matrix(rnorm(60), 10, 6)
  Sigma <- diag(6)
  gamma <- matrix(rnorm(12), 6, 2)

  expect_error(
    dkge_residualize_alignment_features(G, cbind(gamma[, 1], gamma[, 1]), Sigma),
    "estimable span",
    class = "dkge_alignment_rank_error"
  )
  expect_error(
    dkge_residualize_alignment_features(
      G[, 1:3], gamma = diag(3)[, 1:2, drop = FALSE], Sigma_G = diag(3)
    ),
    "Residual feature dimension",
    class = "dkge_alignment_rank_error"
  )
  ill_covariance <- diag(c(1, 1e-12, rep(1, 4)))
  expect_error(
    dkge_residualize_alignment_features(
      G, gamma = diag(6)[, 1:2, drop = FALSE], Sigma_G = ill_covariance
    ),
    class = "dkge_alignment_covariance_error"
  )
  collapse_G <- cbind(rnorm(20), matrix(rnorm(20 * 4, sd = 1e-8), 20, 4))
  expect_error(
    dkge_residualize_alignment_features(
      collapse_G, gamma = matrix(c(1, 0, 0, 0, 0), 5, 1),
      Sigma_G = diag(5),
      control = dkge_alignment_feature_control(
        min_feature_rank = 1L, min_residual_rank = 1L,
        min_effective_rank = 1, min_retained_energy = 0.1
      )
    ),
    class = "dkge_alignment_feature_collapse_error"
  )
  row_collapse_G <- rbind(
    cbind(matrix(rnorm(4), 4, 1), matrix(0, 4, 3)),
    cbind(matrix(rnorm(4), 4, 1), matrix(rnorm(12), 4, 3))
  )
  expect_error(
    dkge_residualize_alignment_features(
      row_collapse_G,
      gamma = matrix(c(1, 0, 0, 0), 4, 1),
      Sigma_G = diag(4),
      control = dkge_alignment_feature_control(
        min_feature_rank = 1L, min_residual_rank = 1L,
        min_effective_rank = 1, min_retained_energy = 0,
        max_row_collapse_fraction = 0.25
      )
    ),
    "collapsed",
    class = "dkge_alignment_feature_collapse_error"
  )
  expect_error(
    dkge_residualize_alignment_features(
      matrix(rnorm(120), 30, 4),
      gamma = matrix(c(1, 0, 0, 0), 4, 1),
      Sigma_G = diag(c(1, 1, 1e-8, 1e-8)),
      control = dkge_alignment_feature_control(
        rank_tolerance = 1e-12,
        min_feature_rank = 1L, min_residual_rank = 1L,
        min_effective_rank = 1.5, min_retained_energy = 0
      )
    ),
    "effective rank",
    class = "dkge_alignment_rank_error"
  )
  expect_error(
    dkge_residualize_alignment_features(G, gamma, Sigma + upper.tri(Sigma)),
    class = "dkge_alignment_covariance_error"
  )
})

test_that("missing covariance and unverified preprocessing are hard failures", {
  fx <- make_alignment_feature_fixture(seed = 260835)
  missing_fit <- fx$fit
  missing_fit$subjects[[1]]$effect_noise_cov <- NULL
  expect_error(
    dkge_alignment_features(
      missing_fit, fx$contrast,
      feature_source = "same_data_residualized"
    ),
    "Missing effect covariance",
    class = "dkge_alignment_covariance_error"
  )

  unverified <- fx$contrast
  unverified$metadata$alignment_receipts[[1]]$preprocessing$effect_separability$transform_verified <- FALSE
  expect_error(
    dkge_alignment_features(
      fx$fit, unverified,
      feature_source = "same_data_residualized"
    ),
    "Unverified alignment preprocessing",
    class = "dkge_alignment_preprocessing_error"
  )
})

test_that("kernel-image residualization has a bounded moderate-size runtime", {
  skip_on_cran()
  set.seed(260836)
  P <- 400L
  k <- 80L
  m <- 4L
  G <- matrix(rnorm(P * k), P, k)
  A <- matrix(rnorm(k * k), k, k)
  Sigma <- crossprod(A) + diag(k)
  gamma <- matrix(rnorm(k * m), k, m)
  elapsed <- system.time({
    result <- dkge_residualize_alignment_features(G, gamma, Sigma)
  })[["elapsed"]]
  expect_true(result$diagnostics$gates_passed)
  expect_lt(elapsed, 10)
})
