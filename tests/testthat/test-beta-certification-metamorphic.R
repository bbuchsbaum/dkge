testthat::local_edition(3)

beta_certification_fixture <- function(seed = 82026L) {
  set.seed(seed)
  q <- 4L
  p <- 7L
  effects <- paste0("e", seq_len(q))
  subjects <- paste0("s", seq_len(5L))
  betas <- lapply(subjects, function(id) {
    out <- matrix(rnorm(q * p), q, p)
    dimnames(out) <- list(effects, paste0("v", seq_len(p)))
    out
  })
  names(betas) <- subjects
  designs <- lapply(subjects, function(id) {
    out <- diag(q)
    colnames(out) <- effects
    out
  })
  names(designs) <- subjects
  A <- matrix(rnorm(q * q), q, q)
  K <- crossprod(A) + 0.25 * diag(q)
  dimnames(K) <- list(effects, effects)
  list(betas = betas, designs = designs, K = K,
       effects = effects, subjects = subjects)
}

beta_certification_fit <- function(fx, betas = fx$betas,
                                   designs = fx$designs, K = fx$K,
                                   omega = NULL) {
  dkge_fit(
    betas, designs, K = K, Omega_list = omega, rank = 3,
    w_method = "none", effect_scaling = "none", keep_X = TRUE
  )
}

beta_projector <- function(fit) {
  Q <- fit$eig_vectors_full[, seq_len(fit$rank), drop = FALSE]
  tcrossprod(Q)
}

test_that("beta geometry is equivariant to labels, voxel bases, and scale", {
  fx <- beta_certification_fixture()
  base <- beta_certification_fit(fx)

  effect_perm <- c(3L, 1L, 4L, 2L)
  permuted <- beta_certification_fit(
    fx,
    betas = lapply(fx$betas, function(B) B[effect_perm, , drop = FALSE]),
    designs = lapply(
      fx$designs,
      function(X) X[, effect_perm, drop = FALSE]
    ),
    K = fx$K[effect_perm, effect_perm, drop = FALSE]
  )
  expect_equal(
    permuted$Chat,
    base$Chat[effect_perm, effect_perm, drop = FALSE],
    tolerance = 1e-9,
    ignore_attr = TRUE
  )
  expect_equal(
    beta_projector(permuted),
    beta_projector(base)[effect_perm, effect_perm, drop = FALSE],
    tolerance = 1e-8,
    ignore_attr = TRUE
  )

  set.seed(82027L)
  rotated_betas <- lapply(fx$betas, function(B) {
    Q <- qr.Q(qr(matrix(rnorm(ncol(B)^2), ncol(B), ncol(B))))
    B %*% Q
  })
  rotated <- beta_certification_fit(fx, betas = rotated_betas)
  expect_equal(rotated$Chat, base$Chat, tolerance = 1e-9,
               ignore_attr = TRUE)
  expect_equal(rotated$eig_values_full, base$eig_values_full,
               tolerance = 1e-9)
  expect_equal(beta_projector(rotated), beta_projector(base),
               tolerance = 1e-8, ignore_attr = TRUE)

  scale <- 7.25
  scaled <- beta_certification_fit(
    fx, betas = lapply(fx$betas, `*`, scale)
  )
  expect_identical(scaled$rank, base$rank)
  expect_equal(scaled$Chat, scale^2 * base$Chat, tolerance = 1e-8,
               ignore_attr = TRUE)
  expect_equal(scaled$eig_values_full, scale^2 * base$eig_values_full,
               tolerance = 1e-8)
  expect_equal(beta_projector(scaled), beta_projector(base),
               tolerance = 1e-8, ignore_attr = TRUE)
})

test_that("diagonal vector and matrix spatial weights are identical", {
  fx <- beta_certification_fixture(seed = 82028L)
  omega_vector <- lapply(fx$subjects, function(id) seq(0.5, 1.7, length.out = 7))
  names(omega_vector) <- fx$subjects
  omega_matrix <- lapply(omega_vector, diag)

  vector_fit <- beta_certification_fit(fx, omega = omega_vector)
  matrix_fit <- beta_certification_fit(fx, omega = omega_matrix)

  expect_equal(matrix_fit$Chat, vector_fit$Chat, tolerance = 1e-10,
               ignore_attr = TRUE)
  expect_equal(matrix_fit$eig_values_full, vector_fit$eig_values_full,
               tolerance = 1e-10)
  expect_equal(beta_projector(matrix_fit), beta_projector(vector_fit),
               tolerance = 1e-9, ignore_attr = TRUE)
})

test_that("streamed prediction agrees with the batch prediction path", {
  fx <- beta_certification_fixture(seed = 82029L)
  fit <- beta_certification_fit(fx)
  contrasts <- diag(length(fx$effects))[, 1:2, drop = FALSE]
  rownames(contrasts) <- fx$effects
  colnames(contrasts) <- c("e1", "e2")

  batch <- dkge_predict(fit, fx$betas, contrasts)
  loader <- list(
    n = function() length(fx$betas),
    B = function(s) fx$betas[[s]]
  )
  streamed <- dkge_predict_stream(fit, loader, contrasts)

  expect_equal(unname(streamed$values), unname(batch$values), tolerance = 0)
  expect_equal(unname(streamed$A_list), unname(batch$A_list), tolerance = 0)
})

test_that("direct moment and exact block-factor decompositions agree", {
  fx <- beta_certification_fixture(seed = 82030L)
  fit <- beta_certification_fit(fx)
  expect_identical(fit$representation, "block_biprojector")

  factorized <- svd(fit$X_concat, nu = nrow(fit$X_concat), nv = 0)
  direct_values <- fit$eig_values_full[seq_along(factorized$d)]
  expect_equal(factorized$d^2, direct_values, tolerance = 1e-9)
  expect_equal(
    tcrossprod(factorized$u[, seq_len(fit$rank), drop = FALSE]),
    beta_projector(fit),
    tolerance = 1e-8,
    ignore_attr = TRUE
  )
})

test_that("fit geometry is invariant to named subject order", {
  fx <- beta_certification_fixture(seed = 82031L)
  base <- beta_certification_fit(fx)
  subject_perm <- c(4L, 1L, 5L, 2L, 3L)
  reordered <- beta_certification_fit(
    fx,
    betas = fx$betas[subject_perm],
    designs = fx$designs[subject_perm]
  )

  expect_equal(reordered$Chat, base$Chat, tolerance = 1e-10,
               ignore_attr = TRUE)
  expect_equal(reordered$eig_values_full, base$eig_values_full,
               tolerance = 1e-10)
  expect_equal(beta_projector(reordered), beta_projector(base),
               tolerance = 1e-9, ignore_attr = TRUE)
})

test_that("between-subject targets preserve named feature-column meaning", {
  set.seed(82032L)
  n <- 18L
  subject_ids <- paste0("s", seq_len(n))
  X <- cbind(intercept = 1, x = rnorm(n), z = rnorm(n))
  rownames(X) <- subject_ids
  coef <- matrix(
    rnorm(ncol(X) * 6L), nrow = ncol(X),
    dimnames = list(colnames(X), paste0("feature", seq_len(6L)))
  )
  Y <- X %*% coef + matrix(rnorm(n * 6L, sd = 0.05), n, 6L)
  dimnames(Y) <- list(subject_ids, colnames(coef))
  base <- dkge_between_rrr(dkge_make_target(Y = Y), X, rank = 2)

  feature_perm <- c(5L, 2L, 6L, 1L, 4L, 3L)
  reordered <- dkge_between_rrr(
    dkge_make_target(Y = Y[, feature_perm, drop = FALSE]), X, rank = 2
  )

  expect_equal(
    reordered$coef[, colnames(base$coef), drop = FALSE],
    base$coef,
    tolerance = 1e-9
  )
  expect_equal(
    reordered$fitted[, colnames(base$fitted), drop = FALSE],
    base$fitted,
    tolerance = 1e-9
  )
})
