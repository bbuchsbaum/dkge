library(testthat)

make_algebra_fixture <- function(S = 3L, q = 3L, P = 4L, seed = 6201) {
  set.seed(seed)
  effects <- paste0("e", seq_len(q))
  betas <- lapply(seq_len(S), function(s) {
    B <- matrix(rnorm(q * P), q, P)
    rownames(B) <- effects
    B
  })
  designs <- lapply(seq_len(S), function(s) {
    X <- matrix(rnorm((q + 4L) * q), q + 4L, q)
    colnames(X) <- effects
    X
  })
  A <- matrix(rnorm(q * q), q, q)
  K <- crossprod(A) + 0.5 * diag(q)
  dimnames(K) <- list(effects, effects)
  list(betas = betas, designs = designs, K = K, effects = effects)
}

test_that("ordinary fits advertise an exact block biprojector", {
  fixture <- make_algebra_fixture()
  omega <- lapply(fixture$betas, function(B) seq(0.5, 1.4, length.out = ncol(B)))
  fit <- dkge_fit(
    fixture$betas, fixture$designs, K = fixture$K, Omega_list = omega,
    rank = 2, w_method = "none", keep_X = TRUE
  )

  expect_equal(fit$representation, "block_biprojector")
  expect_s3_class(fit, "multiblock_biprojector")
  expect_equal(fit$Chat, tcrossprod(fit$X_concat), tolerance = 1e-10,
               ignore_attr = TRUE)
  expect_equal(crossprod(fit$v), diag(fit$rank), tolerance = 1e-9,
               ignore_attr = TRUE)
  expect_equal(crossprod(fit$U, fit$K %*% fit$U), diag(fit$rank),
               tolerance = 1e-9, ignore_attr = TRUE)

  Qr <- fit$eig_vectors_full[, seq_len(fit$rank), drop = FALSE]
  expected_rank_r <- Qr %*% crossprod(Qr, fit$X_concat)
  expect_equal(fit$s %*% t(fit$v), expected_rank_r, tolerance = 1e-9,
               ignore_attr = TRUE)
})

test_that("missing rows are zeroed before non-diagonal R and K mix effects", {
  B1 <- matrix(c(1, 2, 3, 4), 2, 2,
               dimnames = list(c("e1", "e2"), c("v1", "v2")))
  B2 <- matrix(c(5, 1, 2, 6), 2, 2,
               dimnames = list(c("e2", "e3"), c("v1", "v2")))
  X1 <- matrix(c(1, 0, 1, 1, 2, 1, 1, 2, 3, 1), 5, 2,
               dimnames = list(NULL, c("e1", "e2")))
  X2 <- matrix(c(1, 1, 0, 1, 2, 1, 1, 3, 2, 2), 5, 2,
               dimnames = list(NULL, c("e2", "e3")))
  data <- suppressWarnings(dkge_data(list(B1, B2), designs = list(X1, X2),
                                     subject_ids = c("s1", "s2")))
  K <- matrix(c(1.5, 0.4, 0.2,
                0.4, 1.2, 0.3,
                0.2, 0.3, 1.1), 3, 3,
              dimnames = list(data$effects, data$effects))
  fit <- dkge_fit(data, K = K, rank = 2, w_method = "none",
                  missingness = "none", keep_X = TRUE)

  masks <- dkge:::.dkge_obs_masks_from_provenance(
    fit$provenance, fit$subject_ids, nrow(K)
  )
  oracle_blocks <- lapply(seq_along(fit$Braw), function(s) {
    Bzero <- fit$Braw[[s]]
    Bzero[!masks[[s]], ] <- 0
    fit$Khalf %*% t(fit$R) %*% Bzero
  })
  oracle <- do.call(cbind, oracle_blocks)

  expect_equal(fit$X_concat, oracle, tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(fit$Chat, tcrossprod(oracle), tolerance = 1e-10,
               ignore_attr = TRUE)

  old_blocks <- lapply(seq_along(fit$Btil), function(s) {
    obs <- which(masks[[s]])
    out <- matrix(0, nrow(K), ncol(fit$Btil[[s]]))
    out[obs, ] <- fit$Khalf[obs, obs, drop = FALSE] %*%
      fit$Btil[[s]][obs, , drop = FALSE]
    out
  })
  expect_gt(norm(oracle - do.call(cbind, old_blocks), "F"), 1e-3)
})

test_that("pair-normalized moments are q-space-only and fail closed", {
  effects <- c("e1", "e2")
  make_subject <- function(B, counts, id) {
    rownames(B) <- effects
    X <- diag(2)
    colnames(X) <- effects
    suppressWarnings(dkge_subject(B, X, id = id, effect_n = counts))
  }
  subjects <- list(
    make_subject(matrix(c(1, 2, 3, 1), 2, 2), c(e1 = 2, e2 = 10), "s1"),
    make_subject(matrix(c(4, 1, 2, 5), 2, 2), c(e1 = 9, e2 = 3), "s2")
  )
  data <- dkge_data(subjects)
  fit <- dkge_fit(
    data, K = diag(2), rank = 1, w_method = "none",
    effect_scaling = "none", effect_weights = dkge_effect_weights("count")
  )

  expect_s3_class(fit, "dkge_qspace")
  expect_equal(fit$representation, "qspace_moment")
  expect_null(fit$v)
  expect_null(fit$X_concat)
  expect_false(inherits(fit, "multiblock_biprojector"))
  expect_equal(fit$Chat %*% fit$eig_vectors_full,
               sweep(fit$eig_vectors_full, 2, fit$eig_values_full, "*"),
               tolerance = 1e-10, ignore_attr = TRUE)
  Qr <- fit$eig_vectors_full[, seq_len(fit$rank), drop = FALSE]
  expect_equal(tcrossprod(fit$s),
               Qr %*% diag(fit$eig_values_full[seq_len(fit$rank)], fit$rank) %*% t(Qr),
               tolerance = 1e-10, ignore_attr = TRUE)
  expect_error(
    dkge_project_blocks(fit, lapply(subjects, `[[`, "beta")),
    "requires representation='block_biprojector'"
  )
  expect_error(
    dkge_fit(
      data, K = diag(2), rank = 1, w_method = "none", keep_X = TRUE,
      effect_scaling = "none", effect_weights = dkge_effect_weights("count")
    ),
    "has no exact subject-by-voxel block factor"
  )
})

test_that("indefinite analytic moments retain diagnostics but not negative components", {
  effects <- c("e1", "e2")
  X <- diag(2)
  colnames(X) <- effects
  B <- matrix(c(3, 0, 0, 0), 2, 2,
              dimnames = list(effects, c("v1", "v2")))
  make_subject <- function(id) {
    Lambda <- diag(c(0.1, 1), 2)
    dimnames(Lambda) <- list(effects, effects)
    suppressWarnings(dkge_subject(
      B, X, id = id, effect_noise_cov = Lambda,
      residual_variance = c(v1 = 1, v2 = 1), residual_df = 1
    ))
  }
  data <- dkge_data(list(make_subject("s1"), make_subject("s2")))
  fit <- suppressWarnings(dkge_fit(
    data, K = diag(2), rank = 2, w_method = "none",
    effect_scaling = "none", debias = "analytic"
  ))

  expect_equal(fit$representation, "qspace_moment")
  expect_null(fit$v)
  expect_equal(fit$rank, 1)
  expect_lt(min(fit$eig_values_full), 0)
  expect_equal(fit$moment_diagnostics$transformed$negative_count, 1)
  expect_true(all(fit$sdev > 0))
  expect_equal(tcrossprod(fit$s),
               fit$eig_vectors_full[, 1, drop = FALSE] %*%
                 fit$eig_values_full[1] %*%
                 t(fit$eig_vectors_full[, 1, drop = FALSE]),
               tolerance = 1e-12, ignore_attr = TRUE)
  expect_error(
    dkge_project_block(fit, 1, B),
    "requires representation='block_biprojector'"
  )
})

test_that("singular PSD kernels retain a K-orthonormal supported basis", {
  fixture <- make_algebra_fixture(q = 3, seed = 6207)
  K <- matrix(1, 3, 3,
              dimnames = list(fixture$effects, fixture$effects))
  expect_warning(
    dkge_fit(fixture$betas, fixture$designs, K = K, rank = 2,
             w_method = "none", keep_X = TRUE),
    "exceeds effective rank"
  )
  fit <- suppressWarnings(dkge_fit(
    fixture$betas, fixture$designs, K = K, rank = 2,
    w_method = "none", keep_X = TRUE
  ))

  expect_equal(fit$rank, 1)
  expect_equal(crossprod(fit$U, K %*% fit$U), matrix(1),
               tolerance = 1e-8, ignore_attr = TRUE)
  expect_equal(fit$Chat, tcrossprod(fit$X_concat), tolerance = 1e-9,
               ignore_attr = TRUE)
})

test_that("q-space moments are equivariant to effect permutation", {
  fixture <- make_algebra_fixture(q = 4, seed = 6211)
  fit <- dkge_fit(fixture$betas, fixture$designs, K = fixture$K, rank = 2,
                  w_method = "none", effect_scaling = "none")
  perm <- c(3, 1, 4, 2)
  fit_perm <- dkge_fit(
    lapply(fixture$betas, function(B) B[perm, , drop = FALSE]),
    lapply(fixture$designs, function(X) X[, perm, drop = FALSE]),
    K = fixture$K[perm, perm], rank = 2, w_method = "none",
    effect_scaling = "none"
  )

  expect_equal(fit_perm$Chat, fit$Chat[perm, perm], tolerance = 1e-9,
               ignore_attr = TRUE)
  expect_equal(fit_perm$eig_values_full, fit$eig_values_full,
               tolerance = 1e-9)
  Q <- fit$eig_vectors_full[, seq_len(fit$rank), drop = FALSE]
  Qp <- fit_perm$eig_vectors_full[, seq_len(fit_perm$rank), drop = FALSE]
  expect_equal(tcrossprod(Qp), tcrossprod(Q[perm, , drop = FALSE]),
               tolerance = 1e-8, ignore_attr = TRUE)
})
