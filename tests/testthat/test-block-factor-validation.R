make_block_factor_validation_fixture <- function(relative_drift) {
  sdev <- c(20, 1e-5)
  Xstar <- diag(c(sdev[1], sdev[2] * (1 + relative_drift)))
  fit <- list(
    representation = "block_biprojector",
    Chat = diag(sdev^2),
    rank = 2L,
    v = diag(c(1, 1 + relative_drift)),
    s = diag(sdev),
    eig_vectors_full = diag(2),
    sdev = sdev
  )
  list(fit = fit, Xstar = Xstar)
}

test_that("near-rank-deficient block factors pass only within measured backward error", {
  for (relative_drift in c(5e-7, 2e-6)) {
    fixture <- make_block_factor_validation_fixture(relative_drift)
    fit <- fixture$fit
    Xstar <- fixture$Xstar
    diagnostics <- .dkge_block_factor_diagnostics(fit, Xstar)

    # The factorization residual is tiny at Chat's scale, but division by the
    # retained 1e-10 component variance measurably amplifies it in V'V.
    expect_gt(diagnostics$chat_scale, 100)
    expect_lt(diagnostics$chat_relative_error, 1e-14)
    expect_gt(diagnostics$condition_proxy, 1e10)
    expect_gt(diagnostics$orthogonality_error, 1e-7)
    expect_lte(diagnostics$orthogonality_error,
               diagnostics$orthogonality_tolerance)

    # Independent dense-SVD oracle: accepting the factor must not hide a wrong
    # rank-r block reconstruction.
    dense <- svd(Xstar, nu = fit$rank, nv = fit$rank)
    oracle <- dense$u[, seq_len(fit$rank), drop = FALSE] %*%
      diag(dense$d[seq_len(fit$rank)], nrow = fit$rank) %*%
      t(dense$v[, seq_len(fit$rank), drop = FALSE])
    represented <- fit$s %*% t(fit$v)
    relative_error <- norm(represented - oracle, "F") /
      max(1, norm(oracle, "F"))
    expect_lt(relative_error, 1e-12)
    factor_relative_error <- norm(fit$Chat - tcrossprod(Xstar), "F") /
      max(1, norm(fit$Chat, "F"))
    expect_lt(factor_relative_error, 1e-14)
  }
})

test_that("conditioning alone never excuses unrelated loading corruption", {
  Xstar <- diag(c(1, 1e-5))
  fit <- list(
    representation = "block_biprojector",
    Chat = tcrossprod(Xstar),
    rank = 2L,
    v = diag(2),
    eig_vectors_full = diag(2),
    sdev = c(1, 1e-5)
  )

  diagnostics <- .dkge_block_factor_diagnostics(fit, Xstar)
  expect_gt(diagnostics$condition_proxy, 1e9)
  expect_equal(diagnostics$orthogonality_tolerance, 1e-10)
  expect_invisible(.dkge_validate_block_factor(fit, Xstar))

  corrupted <- fit
  corrupted$v[2, 2] <- 1 + 1e-5
  expect_error(
    .dkge_validate_block_factor(corrupted, Xstar),
    "not orthonormal.*scale-aware bound"
  )
})

test_that("block-factor validation rejects malformed and invalid factors", {
  Xstar <- diag(2)
  fit <- list(
    representation = "block_biprojector",
    Chat = diag(2),
    rank = 2L,
    v = diag(2),
    eig_vectors_full = diag(2),
    sdev = c(1, 1)
  )
  expect_invisible(.dkge_validate_block_factor(fit, Xstar))

  nonsymmetric <- fit
  nonsymmetric$Chat[1, 2] <- 0.1
  expect_error(
    .dkge_validate_block_factor(nonsymmetric, Xstar),
    "Chat is not symmetric"
  )

  wrong_factor <- Xstar
  wrong_factor[1, 1] <- 1.1
  expect_error(
    .dkge_validate_block_factor(fit, wrong_factor),
    "Chat is not Xstar Xstar"
  )

  wrong_shape <- fit
  wrong_shape$v <- matrix(1, nrow = 1, ncol = 2)
  expect_error(
    .dkge_validate_block_factor(wrong_shape, Xstar),
    "ncol\\(Xstar\\)-by-rank"
  )

  zero_scale <- fit
  zero_scale$sdev[2] <- 0
  expect_error(
    .dkge_validate_block_factor(zero_scale, Xstar),
    "positive finite sdev"
  )

  underflow_scale <- fit
  underflow_scale$sdev[2] <- 1e-200
  expect_error(
    .dkge_validate_block_factor(underflow_scale, Xstar),
    "positive finite sdev"
  )

  wrong_eigenvectors <- fit
  wrong_eigenvectors$eig_vectors_full <- matrix(1, nrow = 1, ncol = 2)
  expect_error(
    .dkge_validate_block_factor(wrong_eigenvectors, Xstar),
    "eigenvectors do not match"
  )

  nonfinite <- fit
  nonfinite$v[1, 1] <- NA_real_
  expect_error(
    .dkge_validate_block_factor(nonfinite, Xstar),
    "finite ncol\\(Xstar\\)-by-rank"
  )
})
