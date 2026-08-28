library(testthat)

set.seed(100)
S <- 3
q <- 3
P <- 4
T <- 40

betas <- replicate(S, matrix(rnorm(q * P), q, P), simplify = FALSE)
designs <- replicate(S, qr.Q(qr(matrix(rnorm(T * q), T, q))), simplify = FALSE)
centroids <- replicate(S, matrix(runif(P * 3), P, 3), simplify = FALSE)

fit <- dkge(betas, designs = designs, K = diag(q), rank = 2)
fit$centroids <- centroids

transport_loadings <- suppressWarnings(dkge_transport_loadings_to_medoid(fit,
                                                        medoid = 1,
                                                        centroids = centroids,
                                                        mapper = dkge_mapper_spec("sinkhorn", epsilon = 0.05)))
legacy_cache <- transport_loadings$cache

contrasts <- c(1, -1, 0)

contrast_obj <- suppressWarnings(dkge_contrast(
  fit, contrasts, method = "loso", align = FALSE
))
set.seed(101)
independent_betas <- lapply(fit$Braw, function(Bs) {
  matrix(rnorm(length(Bs)), nrow(Bs), ncol(Bs), dimnames = dimnames(Bs))
})
alignment_features <- dkge_alignment_features(
  fit,
  contrast_obj,
  independent_betas = independent_betas,
  independent_data_hash = "bootstrap-independent-fixture"
)
cache <- dkge_prepare_alignment(
  fit,
  alignment_features,
  centroids = centroids,
  mapper = dkge_mapper_spec("ridge", lambda = 1e-3),
  reference_subject = 1L
)
aligned_maps <- dkge:::.dkge_apply_fitted_alignment(
  contrast_obj$values,
  cache,
  contrast_obj = contrast_obj,
  estimand = list(target = "mean aligned cross-fitted contrast"),
  application_context = "bootstrap_test_fixture"
)
subject_maps <- aligned_maps$values[[1L]]
values_medoid <- lapply(seq_len(nrow(subject_maps)), function(i) subject_maps[i, ])

vox_map <- diag(ncol(subject_maps))

test_that("Poisson multiplier replicates condition on a non-empty cohort", {
  calls <- 0L
  zero_then_positive <- function(n, lambda) {
    calls <<- calls + 1L
    if (calls == 1L) rep(0, n) else c(1, rep(0, n - 1L))
  }
  draw <- dkge:::.dkge_bootstrap_multipliers(
    "poisson", 4L, poisson_draw = zero_then_positive
  )
  expect_identical(calls, 2L)
  expect_identical(draw, c(1, 0, 0, 0))
  expect_error(
    dkge:::.dkge_bootstrap_multipliers(
      "poisson", 4L, poisson_draw = function(n, lambda) rep(0, n),
      max_redraws = 2L
    ),
    "non-empty Poisson multiplier cohort",
    class = "dkge_bootstrap_weighting_error"
  )
  expect_error(
    dkge:::.dkge_equal_subject_multiplier_mean(
      matrix(seq_len(12), 4L, 3L), rep(0, 4L)
    ),
    "non-empty subject cohort",
    class = "dkge_bootstrap_weighting_error"
  )
})



test_that("projection bootstrap returns expected shapes", {
  expect_error(
    dkge_bootstrap_projected(values_medoid, B = 20),
    "requires typed",
    class = "dkge_alignment_ineligible_error"
  )
  expect_error(
    dkge_bootstrap_projected(aligned_maps, B = 20),
    class = "dkge_alignment_ineligible_error"
  )
  boot <- dkge_bootstrap_projected(
    aligned_maps,
    B = 20,
    voxel_operator = vox_map,
    allow_approximate_alignment = TRUE
  )
  expect_equal(length(boot$medoid$mean), ncol(subject_maps))
  expect_equal(dim(boot$medoid$boot), c(20, ncol(subject_maps)))
  expect_equal(dim(boot$voxel$boot), c(20, ncol(subject_maps)))
  expect_identical(boot$metadata$alignment_status, "approximate")
  expect_true(boot$metadata$approximate_override)
})


test_that("q-space bootstrap runs with cached transport", {
  expect_error(
    dkge_bootstrap_qspace(
      fit, contrasts = contrasts, B = 2,
      transport_cache = list(
        operators = lapply(fit$Btil, function(Bs) diag(ncol(Bs)))
      ),
      allow_approximate_alignment = TRUE
    ),
    "typed",
    class = "dkge_alignment_ineligible_error"
  )
  expect_error(
    dkge_bootstrap_qspace(
      fit, contrasts = contrasts, B = 2,
      transport_cache = legacy_cache,
      allow_approximate_alignment = TRUE
    ),
    class = "dkge_alignment_ineligible_error"
  )
  mutated_cache <- cache
  mutated_cache$operators[[1L]][1L, 1L] <-
    mutated_cache$operators[[1L]][1L, 1L] + 0.25
  expect_error(
    dkge_bootstrap_qspace(
      fit, contrasts = contrasts, B = 2,
      transport_cache = mutated_cache,
      allow_approximate_alignment = TRUE
    ),
    class = "dkge_alignment_cache_mismatch"
  )
  expect_error(
    dkge_bootstrap_qspace(
      fit, contrasts = contrasts, B = 2, transport_cache = cache
    ),
    class = "dkge_alignment_ineligible_error"
  )
  boot_q <- dkge_bootstrap_qspace(fit,
                                  contrasts = contrasts,
                                  B = 10,
                                  seed = 123,
                                  transport_cache = cache,
                                  medoid = 1,
                                  voxel_operator = vox_map,
                                  scheme = "poisson",
                                  allow_approximate_alignment = TRUE)
  expect_equal(boot_q$B, 10)
  expect_equal(length(boot_q$summary), 1)
  expect_equal(ncol(boot_q$summary[[1]]$boot_medoid), ncol(subject_maps))
  expect_true(boot_q$summary[[1]]$medoid$sd[1] >= 0)
  expect_identical(boot_q$metadata$alignment_status, "approximate")
})

test_that("q-space bootstrap separates moment and equal-subject group weights", {
  expect_gt(diff(range(fit$weights)), 0.5)
  seed <- 260858L
  boot_q <- dkge_bootstrap_qspace(
    fit,
    contrasts = contrasts,
    B = 1L,
    scheme = "exp",
    seed = seed,
    transport_cache = cache,
    allow_approximate_alignment = TRUE
  )

  set.seed(seed)
  xi <- dkge:::.dkge_bootstrap_multipliers("exp", S)
  q_dim <- nrow(fit$U)
  rank <- ncol(fit$U)
  contribution_matrix <- vapply(
    fit$contribs, function(M) as.numeric(M), numeric(q_dim * q_dim)
  )
  Chat_b <- matrix(
    contribution_matrix %*% (as.numeric(fit$weights) * xi),
    q_dim, q_dim
  )
  Chat_b <- (Chat_b + t(Chat_b)) / 2
  eig <- eigen(Chat_b, symmetric = TRUE)
  Ub <- fit$Kihalf %*% eig$vectors[, seq_len(rank), drop = FALSE]
  Ub <- dkge_k_orthonormalize(Ub, fit$K)
  Ub <- dkge_procrustes_K(
    fit$U, Ub, fit$K, allow_reflection = FALSE
  )$U_aligned
  corr_diag <- diag(t(fit$U) %*% fit$K %*% Ub)
  Ub <- sweep(Ub, 2L, ifelse(corr_diag < 0, -1, 1), `*`)

  normalized <- dkge:::.normalize_contrasts(contrasts, fit)[[1L]]
  ctil <- backsolve(fit$R, normalized, transpose = FALSE)
  alpha <- as.numeric(crossprod(Ub, fit$K %*% ctil))
  Bmodel <- lapply(seq_along(fit$Btil), function(s) {
    dkge:::.dkge_apply_fit_spatial(fit, fit$Btil[[s]], subject = s)
  })
  subject_maps_oracle <- do.call(rbind, lapply(seq_len(S), function(s) {
    loading <- t(fit$K %*% Bmodel[[s]]) %*% Ub
    value <- as.numeric(loading %*% alpha)
    as.numeric(t(cache$operators[[s]]) %*% value)
  }))
  equal_subject_oracle <- dkge:::.dkge_equal_subject_multiplier_mean(
    subject_maps_oracle, xi
  )
  leaked_mfa_oracle <- colSums(
    subject_maps_oracle * (xi * as.numeric(fit$weights))
  ) / (sum(xi * as.numeric(fit$weights)) + 1e-12)

  expect_equal(boot_q$summary[[1L]]$boot[1L, ], equal_subject_oracle,
               tolerance = 1e-10)
  expect_false(isTRUE(all.equal(
    boot_q$summary[[1L]]$boot[1L, ], leaked_mfa_oracle,
    tolerance = 1e-8
  )))
  expect_identical(boot_q$metadata$group_weighting, "equal_subject")
  expect_identical(boot_q$metadata$moment_weighting$method, fit$w_method)
  expect_identical(
    boot_q$metadata$moment_weighting$weights,
    unname(fit$weights)
  )

  boot_a <- dkge_bootstrap_analytic(
    fit,
    contrasts = contrasts,
    B = 1L,
    scheme = "exp",
    seed = seed,
    transport_cache = cache,
    allow_approximate_alignment = TRUE
  )
  expect_identical(boot_a$metadata$group_weighting, "equal_subject")
  expect_identical(boot_a$metadata$moment_weighting,
                   boot_q$metadata$moment_weighting)
})

test_that("q-space bootstrap rejects a typed cache bound to another fit", {
  perturbed_betas <- betas
  perturbed_betas[[1L]][1L, 1L] <- perturbed_betas[[1L]][1L, 1L] + 0.5
  other_fit <- dkge(
    perturbed_betas,
    designs = designs,
    K = diag(q),
    rank = 2
  )
  other_fit$centroids <- centroids

  expect_identical(other_fit$subject_ids, fit$subject_ids)
  expect_equal(vapply(other_fit$Btil, ncol, integer(1)),
               vapply(fit$Btil, ncol, integer(1)))
  expect_false(identical(
    dkge:::.dkge_alignment_fit_binding(other_fit),
    dkge:::.dkge_alignment_fit_binding(fit)
  ))
  expect_error(
    dkge_bootstrap_qspace(
      other_fit,
      contrasts = contrasts,
      B = 2,
      transport_cache = cache,
      allow_approximate_alignment = TRUE
    ),
    "does not match the fit bound",
    class = "dkge_alignment_cache_mismatch"
  )

  rank_changed <- dkge(
    betas,
    designs = designs,
    K = diag(q),
    rank = 1
  )
  rank_changed$centroids <- centroids
  expect_false(identical(
    dkge:::.dkge_alignment_fit_binding(rank_changed),
    dkge:::.dkge_alignment_fit_binding(fit)
  ))
  expect_error(
    dkge_bootstrap_qspace(
      rank_changed,
      contrasts = contrasts,
      B = 2,
      transport_cache = cache,
      allow_approximate_alignment = TRUE
    ),
    "does not match the fit bound",
    class = "dkge_alignment_cache_mismatch"
  )

  root_changed <- fit
  root_changed$Kihalf[1L, 2L] <- root_changed$Kihalf[1L, 2L] + 0.4
  pool_cache_changed <- fit
  pool_cache_changed$pool_cache$structural[[1L]][1L, 1L] <-
    pool_cache_changed$pool_cache$structural[[1L]][1L, 1L] + 0.25
  for (changed in list(root_changed, pool_cache_changed)) {
    expect_false(identical(
      dkge:::.dkge_alignment_fit_binding(changed),
      dkge:::.dkge_alignment_fit_binding(fit)
    ))
    expect_error(
      dkge_bootstrap_qspace(
        changed,
        contrasts = contrasts,
        B = 2,
        transport_cache = cache,
        allow_approximate_alignment = TRUE
      ),
      "does not match the fit bound",
      class = "dkge_alignment_cache_mismatch"
    )
  }

  spoofed <- cache
  spoofed$fit_binding <- dkge:::.dkge_alignment_fit_binding(other_fit)
  expect_error(
    dkge_bootstrap_qspace(
      other_fit,
      contrasts = contrasts,
      B = 2,
      transport_cache = spoofed,
      allow_approximate_alignment = TRUE
    ),
    "bindings must live only",
    class = "dkge_alignment_cache_mismatch"
  )
})

test_that("spatial solve factors are validated and bound to bootstrap caches", {
  adjacency <- Matrix::sparseMatrix(
    i = c(1L, 2L, 2L, 3L, 3L, 4L),
    j = c(2L, 1L, 3L, 2L, 4L, 3L),
    x = 1,
    dims = c(P, P)
  )
  laplacian <- Matrix::Diagonal(x = Matrix::rowSums(adjacency)) - adjacency
  spatial_labels <- paste0("parcel", seq_len(P))
  dimnames(laplacian) <- list(spatial_labels, spatial_labels)
  spatial_betas <- lapply(betas, function(Bs) {
    colnames(Bs) <- spatial_labels
    Bs
  })
  spatial_independent_betas <- lapply(independent_betas, function(Bs) {
    colnames(Bs) <- spatial_labels
    Bs
  })
  spatial_fit <- dkge(
    spatial_betas,
    designs = designs,
    K = diag(q),
    rank = 2,
    spatial = dkge_spatial_regularizer(laplacian = laplacian, lambda = 0.4)
  )
  spatial_contrast <- suppressWarnings(dkge_contrast(
    spatial_fit, contrasts, method = "loso", align = FALSE
  ))
  spatial_features <- dkge_alignment_features(
    spatial_fit,
    spatial_contrast,
    independent_betas = spatial_independent_betas,
    independent_data_hash = "bootstrap-spatial-independent-fixture"
  )
  spatial_cache <- dkge_prepare_alignment(
    spatial_fit,
    spatial_features,
    centroids = centroids,
    mapper = dkge_mapper_spec("ridge", lambda = 1e-3),
    reference_subject = 1L
  )

  operator_spoofs <- list(
    effective = function(x) { x$effective <- FALSE; x },
    edge_count = function(x) { x$n_edges <- 0L; x },
    fractional_edge_count = function(x) { x$n_edges <- x$n_edges + 0.9; x },
    domain_permutation = function(x) { x$domain <- rev(x$domain); x },
    status = function(x) { x$status <- "inert"; x }
  )
  for (mutate_operator in operator_spoofs) {
    spoofed_operator <- mutate_operator(spatial_fit$spatial$operators[[1L]])
    expect_error(
      dkge:::.dkge_spatial_apply_betas(
        spatial_fit$Btil[[1L]], spoofed_operator
      ),
      "effectiveness metadata",
      class = "dkge_spatial_provenance_error"
    )
  }

  fit_spoofs <- list(
    active = function(x) { x$spatial$active <- FALSE; x },
    requested = function(x) { x$spatial$requested <- FALSE; x },
    effective = function(x) { x$spatial$effective <- FALSE; x },
    fully_effective = function(x) {
      x$spatial$fully_effective <- FALSE
      x
    },
    partial_status = function(x) { x$spatial$status <- "partial"; x },
    edge_diagnostic = function(x) {
      x$spatial$diagnostics$n_edges[[1L]] <-
        x$spatial$diagnostics$n_edges[[1L]] + 1L
      x
    },
    operator_effective = function(x) {
      x$spatial$operators[[1L]]$effective <- FALSE
      x
    },
    spec_diagonal_drift = function(x) {
      x$spatial$spec$laplacians[[1L]] <-
        x$spatial$spec$laplacians[[1L]] +
        Matrix::Diagonal(P, x = rep(0.1, P))
      x
    },
    spec_asymmetry = function(x) {
      L <- Matrix::Matrix(
        as.matrix(x$spatial$spec$laplacians[[1L]]), sparse = TRUE
      )
      L[1L, 2L] <- L[1L, 2L] + 0.1
      x$spatial$spec$laplacians[[1L]] <- L
      x
    },
    spec_domain = function(x) {
      x$spatial$spec$domains[[1L]][[1L]] <-
        x$spatial$spec$domains[[1L]][[2L]]
      x
    },
    spec_domain_removed = function(x) {
      x$spatial$spec$domains[1L] <- list(NULL)
      x
    },
    spec_foreign_entry = function(x) {
      x$spatial$spec$foreign_entry <- TRUE
      x
    }
  )
  for (mutate_fit in fit_spoofs) {
    spoofed_fit <- mutate_fit(spatial_fit)
    expect_error(
      dkge:::.dkge_alignment_fit_binding(spoofed_fit),
      class = "dkge_error"
    )
    expect_error(
      dkge:::.dkge_validate_alignment_features(
        spatial_features, fit = spoofed_fit,
        contrast_obj = spatial_contrast
      ),
      class = "dkge_error"
    )
    expect_error(
      dkge_bootstrap_qspace(
        spoofed_fit,
        contrasts = contrasts,
        B = 2,
        transport_cache = spatial_cache,
        allow_approximate_alignment = TRUE
      ),
      class = "dkge_error"
    )
  }

  mutated <- spatial_fit
  op <- mutated$spatial$operators[[1L]]
  incompatible_A <- Matrix::forceSymmetric(
    Matrix::Diagonal(P) + 0.9 * op$L,
    uplo = "U"
  )
  mutated$spatial$operators[[1L]]$factor <- Matrix::Cholesky(
    incompatible_A, LDL = FALSE, perm = TRUE
  )

  expect_error(
    dkge:::.dkge_alignment_fit_binding(mutated),
    "inconsistent with `I \\+ lambda \\* L`",
    class = "dkge_spatial_factor_error"
  )
  expect_error(
    dkge:::.dkge_apply_fit_spatial(
      mutated, mutated$Btil[[1L]], subject = 1L
    ),
    class = "dkge_spatial_factor_error"
  )
  expect_error(
    dkge_bootstrap_qspace(
      mutated,
      contrasts = contrasts,
      B = 2,
      transport_cache = spatial_cache,
      allow_approximate_alignment = TRUE
    ),
    class = "dkge_spatial_factor_error"
  )
})

test_that("residualized bootstrap cache is bound to its contrast family", {
  residual_features <- dkge_alignment_features(
    fit,
    contrast_obj,
    feature_source = "same_data_residualized",
    effect_noise_cov = replicate(S, diag(q), simplify = FALSE),
    control = dkge_alignment_feature_control(min_effective_rank = 1.0)
  )
  residual_cache <- dkge_prepare_alignment(
    fit,
    residual_features,
    centroids = centroids,
    mapper = dkge_mapper_spec("ridge", lambda = 1e-3),
    reference_subject = 1L
  )
  expect_error(
    dkge_bootstrap_qspace(
      fit,
      contrasts = c(0, 1, -1),
      B = 2,
      transport_cache = residual_cache,
      allow_approximate_alignment = TRUE
    ),
    "different contrast family",
    class = "dkge_alignment_cache_mismatch"
  )
})

test_that("independent bootstrap cache is bound to its contrast family", {
  other_contrast <- c(0, 1, -1)
  expect_error(
    dkge_bootstrap_qspace(
      fit,
      contrasts = other_contrast,
      B = 2,
      transport_cache = cache,
      allow_approximate_alignment = TRUE
    ),
    "different contrast family",
    class = "dkge_alignment_cache_mismatch"
  )
  expect_error(
    dkge_bootstrap_analytic(
      fit,
      contrasts = other_contrast,
      B = 2,
      transport_cache = cache,
      allow_approximate_alignment = TRUE
    ),
    "different contrast family",
    class = "dkge_alignment_cache_mismatch"
  )
})


test_that("analytic bootstrap falls back gracefully", {
  expect_error(
    dkge_bootstrap_analytic(
      fit,
      contrasts = contrasts,
      B = 2,
      transport_cache = cache
    ),
    class = "dkge_alignment_ineligible_error"
  )
  boot_a <- dkge_bootstrap_analytic(fit,
                                    contrasts = contrasts,
                                    B = 8,
                                    seed = 321,
                                    transport_cache = cache,
                                    medoid = 1,
                                    scheme = "bayes",
                                    voxel_operator = vox_map,
                                    perturb_tol = 0.5,
                                    allow_approximate_alignment = TRUE)
  expect_equal(boot_a$B, 8)
  expect_true(boot_a$fallbacks >= 0)
  expect_equal(ncol(boot_a$summary[[1]]$boot_medoid), ncol(subject_maps))
  expect_identical(boot_a$metadata$estimator_status, "approximate")
})

test_that("q-space bootstrap refuses a typed nonconverged alignment", {
  failed_mapper <- dkge_mapper_spec(
    "sinkhorn", epsilon = 1e-6, lambda_emb = 1, lambda_spa = 1,
    max_iter = 1L, tol = 1e-14, warm_start = FALSE
  )
  expect_error(
    suppressWarnings(dkge_prepare_alignment(
      fit,
      alignment_features,
      centroids = centroids,
      mapper = failed_mapper,
      reference_subject = 1L
    )),
    "numerically invalid|converge",
    class = "dkge_alignment_numerical_error"
  )
})
