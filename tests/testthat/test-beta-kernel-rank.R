beta_rank_fixture <- function(scale = 1, K = diag(2), rank = 2L) {
  effects <- c("e1", "e2")
  B <- scale * diag(c(1, 1e-5))
  rownames(B) <- effects
  X <- diag(2)
  colnames(X) <- effects
  dimnames(K) <- list(effects, effects)
  suppressWarnings(dkge_fit(
    list(s1 = B, s2 = B),
    list(s1 = X, s2 = X),
    K = K,
    rank = rank,
    w_method = "none",
    effect_scaling = "none"
  ))
}

test_that("singular PSD kernels use only their numerical range", {
  fit <- beta_rank_fixture(K = diag(c(1, 0)))

  expect_equal(fit$kernel_rank, 1L)
  expect_equal(fit$kernel_nullity, 1L)
  expect_equal(fit$rank, 1L)
  expect_equal(crossprod(fit$U, fit$K %*% fit$U), diag(1), tolerance = 1e-9)
  expect_equal(unname(fit$Khalf), diag(c(1, 0)), tolerance = 1e-12)
  expect_equal(unname(fit$Kihalf), diag(c(1, 0)), tolerance = 1e-12)
  expect_named(
    fit$spectral_diagnostics,
    c("absolute_tolerance", "relative_tolerance", "kernel_tolerance",
      "moment_tolerance", "kernel_rank", "kernel_nullity",
      "moment_rank", "effective_rank"),
    ignore.order = TRUE
  )
  expect_identical(dkge_diagnostics(fit)$spectral,
                   fit$spectral_diagnostics)
})

test_that("zero kernels fail before an arbitrary component is constructed", {
  expect_error(
    beta_rank_fixture(K = matrix(0, 2, 2)),
    "[Kk]ernel.*rank zero|zero-rank kernel",
    class = "dkge_zero_rank_error"
  )
})

test_that("effective rank and whitened subspace are scale metamorphic", {
  scales <- 10^seq(-10, 10, by = 5)
  fits <- lapply(scales, beta_rank_fixture)
  ranks <- vapply(fits, `[[`, integer(1), "rank")

  expect_length(unique(ranks), 1L)
  reference <- fits[[1]]$Khalf %*% fits[[1]]$U
  reference_projector <- tcrossprod(reference)
  for (fit in fits[-1]) {
    whitened <- fit$Khalf %*% fit$U
    expect_equal(tcrossprod(whitened), reference_projector, tolerance = 1e-8)
    expect_equal(
      crossprod(fit$U, fit$K %*% fit$U),
      diag(fit$rank),
      tolerance = 1e-8
    )
  }
})

test_that("kernel rescaling preserves range rank and geometric subspace", {
  effects <- c("e1", "e2")
  B <- diag(c(1, 0.4))
  rownames(B) <- effects
  X <- diag(2)
  colnames(X) <- effects
  K0 <- diag(c(2, 0.5))
  scales <- 10^seq(-10, 10, by = 5)
  fits <- lapply(scales, function(scale) {
    K <- scale * K0
    dimnames(K) <- list(effects, effects)
    dkge_fit(
      list(s1 = B, s2 = B), list(s1 = X, s2 = X),
      K = K, rank = 2, w_method = "none", effect_scaling = "none"
    )
  })

  expect_identical(vapply(fits, `[[`, integer(1), "kernel_rank"),
                   rep(2L, length(scales)))
  expect_identical(vapply(fits, `[[`, integer(1), "rank"),
                   rep(2L, length(scales)))
  reference <- tcrossprod(fits[[1]]$Khalf %*% fits[[1]]$U)
  for (fit in fits[-1]) {
    expect_equal(tcrossprod(fit$Khalf %*% fit$U), reference,
                 tolerance = 1e-8)
  }
})

test_that("ill-conditioned kernels discard unsupported metric directions", {
  fit <- beta_rank_fixture(K = diag(c(1, 1e-12)))
  expect_equal(fit$kernel_rank, 1L)
  expect_equal(fit$kernel_nullity, 1L)
  expect_equal(fit$rank, 1L)
  expect_equal(crossprod(fit$U, fit$K %*% fit$U), diag(1),
               tolerance = 1e-9)
})

test_that("fold bases cap rank by their training moment", {
  effects <- c("e1", "e2")
  B1 <- matrix(c(1, 0, 0, 0), nrow = 2,
               dimnames = list(effects, c("p1", "p2")))
  B2 <- matrix(c(2, 0, 0, 0), nrow = 2,
               dimnames = list(effects, c("p1", "p2")))
  B3 <- matrix(c(0, 1, 0, 0), nrow = 2,
               dimnames = list(effects, c("p1", "p2")))
  X <- diag(2)
  colnames(X) <- effects
  fit <- suppressWarnings(dkge_fit(
    list(s1 = B1, s2 = B2, s3 = B3),
    list(s1 = X, s2 = X, s3 = X),
    K = diag(2), rank = 2,
    w_method = "none", effect_scaling = "none"
  ))

  fold <- .dkge_build_fold_bases(
    fit, assignments = list(3L), ridge = 0,
    align = FALSE, loader_scope = "heldout"
  )$folds[[1]]

  expect_equal(fold$training_rank, 1L)
  expect_equal(ncol(fold$basis), 1L)
  expect_equal(crossprod(fold$basis, fit$K %*% fold$basis),
               diag(1), tolerance = 1e-9)
})
