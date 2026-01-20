# test-fit.R
# Direct diagnostics for dkge_fit helpers

library(testthat)

make_fit_fixture <- function(S = 3, q = 3, P = 4, T = 20, seed = 900) {
  set.seed(seed)
  betas <- replicate(S, matrix(rnorm(q * P), q, P), simplify = FALSE)
  designs <- replicate(S, {
    X <- matrix(rnorm(T * q), T, q)
    qr.Q(qr(X))
  }, simplify = FALSE)
  list(betas = betas, designs = designs, K = diag(q))
}

test_that("pooled design Cholesky matches Gram matrix", {
  fixture <- make_fit_fixture()
  ruler <- dkge:::.dkge_compute_shared_ruler(fixture$designs)
  expect_equal(ruler$G_pool, Reduce(`+`, lapply(fixture$designs, crossprod)))
  expect_equal(t(ruler$R) %*% ruler$R, ruler$G_pool, tolerance = 1e-10)
})

test_that("row standardisation matches manual computation", {
  fixture <- make_fit_fixture()
  ruler <- dkge:::.dkge_compute_shared_ruler(fixture$designs)
  Btil <- dkge:::.dkge_row_standardize(fixture$betas, ruler$R)
  expect_equal(Btil[[1]], t(ruler$R) %*% fixture$betas[[1]])
})

test_that("kernel roots reconstruct original kernel", {
  K <- matrix(c(2, 0.5, 0.5, 1.5), 2, 2)
  roots <- dkge:::.dkge_kernel_roots(K)
  expect_equal(roots$Khalf %*% roots$Khalf, K, tolerance = 1e-10)
  expect_equal(roots$Khalf %*% roots$Kihalf, diag(2), tolerance = 1e-10)
})

test_that("subject weights handle MFA and Omega matrices", {
  Btil <- list(matrix(1:6, 3, 2), matrix(2:7, 3, 2))
  Omega_list <- list(diag(c(1, 2)), matrix(c(2, 1, 1, 3), 2))
  Khalf <- diag(3)

  w_mfa <- dkge:::.dkge_subject_weights(Btil, list(NULL, NULL), Khalf,
                                        w_method = "mfa_sigma1", w_tau = 0)
  expect_true(all(is.finite(w_mfa)))

  w_mat <- dkge:::.dkge_subject_weights(Btil, Omega_list, Khalf,
                                        w_method = "energy", w_tau = 0.5)
  expect_equal(length(w_mat), length(Btil))
  expect_true(all(w_mat > 0))
})

test_that("dkge_fit accepts raw lists and dkge_data", {
  fixture <- make_fit_fixture()
  fit1 <- dkge_fit(fixture$betas, fixture$designs, K = fixture$K, rank = 2,
                   keep_X = TRUE)
  data_bundle <- dkge_data(fixture$betas, designs = fixture$designs)
  fit2 <- dkge_fit(data_bundle, K = fixture$K, rank = 2, keep_X = FALSE)

  expect_equal(fit1$rank, 2)
  expect_equal(fit2$rank, 2)
  expect_true(!is.null(fit1$X_concat))
  expect_null(fit2$X_concat)
  expect_equal(fit1$K, fixture$K)
})

test_that("dkge_fit honours ridge and Omega weighting", {
  fixture <- make_fit_fixture()
  Omega <- lapply(fixture$betas, function(B) runif(ncol(B)))
  fit <- dkge_fit(fixture$betas, fixture$designs, K = fixture$K, rank = 2,
                  ridge = 0.5, Omega_list = Omega)
  expect_true(all(eigen(fit$Chat, symmetric = TRUE)$values >= -1e-8))
  expect_equal(length(fit$weights), length(fixture$betas))
})

test_that("cpca arguments are respected", {
  fixture <- make_fit_fixture(q = 5)
  fit <- dkge_fit(fixture$betas, fixture$designs, K = fixture$K, rank = 3,
                  cpca_blocks = 1:2, cpca_part = "design", cpca_ridge = 0.1)
  expect_true(!is.null(fit$cpca))
  expect_equal(fit$cpca$part, "design")
  expect_equal(fit$cpca$blocks, 1:2)
})

test_that("dkge_fit detects dimension mismatches", {
  fixture <- make_fit_fixture()
  bad_designs <- fixture$designs
  bad_designs[[1]] <- bad_designs[[1]][, -1]
  expect_error(dkge_fit(fixture$betas, bad_designs, K = fixture$K, rank = 2))
})

# -------------------------------------------------------------------------
# K-orthonormality property tests
# -------------------------------------------------------------------------

test_that("K-orthonormality holds with identity kernel", {
  withr::local_seed(1001)
  fixture <- make_fit_fixture(S = 3, q = 4, P = 5, T = 20)
  fit <- dkge_fit(fixture$betas, fixture$designs, K = fixture$K, rank = 3,
                  w_method = "none")

  # U^T K U = I_r (identity kernel K=I so this is simply orthonormality)
  UtKU <- t(fit$U) %*% fit$K %*% fit$U
  expect_lt(max(abs(UtKU - diag(fit$rank))), 1e-8)
})

test_that("K-orthonormality holds with RBF (ordinal) kernel", {
  withr::local_seed(1002)
  q <- 5
  # Construct an RBF kernel for ordinal data
  dists <- as.matrix(dist(1:q))
  sigma <- 1.5
  K <- exp(-dists^2 / (2 * sigma^2))
  K <- (K + t(K)) / 2  # ensure symmetric

  fixture <- make_fit_fixture(S = 4, q = q, P = 6, T = 25)
  fit <- dkge_fit(fixture$betas, fixture$designs, K = K, rank = 3,
                  w_method = "none")

  UtKU <- t(fit$U) %*% fit$K %*% fit$U
  expect_lt(max(abs(UtKU - diag(fit$rank))), 1e-8)
})

test_that("K-orthonormality holds with multi-factor kernel", {
  withr::local_seed(1003)
  q <- 6
  # Multi-factor kernel: block diagonal with 2 factors
  K <- matrix(0, q, q)
  K[1:3, 1:3] <- 1
  K[4:6, 4:6] <- 1
  diag(K) <- 1
  K <- K + 0.1 * diag(q)  # ensure PSD

  fixture <- make_fit_fixture(S = 3, q = q, P = 5, T = 20)
  fit <- dkge_fit(fixture$betas, fixture$designs, K = K, rank = 4,
                  w_method = "none")

  UtKU <- t(fit$U) %*% fit$K %*% fit$U
  expect_lt(max(abs(UtKU - diag(fit$rank))), 1e-8)
})

test_that("K-orthonormality holds for rank=1 edge case", {
  withr::local_seed(1004)
  fixture <- make_fit_fixture(S = 3, q = 4, P = 5, T = 20)
  fit <- dkge_fit(fixture$betas, fixture$designs, K = fixture$K, rank = 1,
                  w_method = "none")

  expect_equal(fit$rank, 1)
  UtKU <- t(fit$U) %*% fit$K %*% fit$U
  expect_lt(max(abs(UtKU - diag(fit$rank))), 1e-8)
})

test_that("K-orthonormality holds for full rank (rank=q)", {
  withr::local_seed(1005)
  q <- 4
  fixture <- make_fit_fixture(S = 3, q = q, P = 5, T = 20)
  fit <- dkge_fit(fixture$betas, fixture$designs, K = fixture$K, rank = q,
                  w_method = "none")

  # May be less than q if some eigenvalues are near zero
  expect_lte(fit$rank, q)
  UtKU <- t(fit$U) %*% fit$K %*% fit$U
  expect_lt(max(abs(UtKU - diag(fit$rank))), 1e-8)
})

test_that("K-orthonormality holds for near-full rank (rank=q-1)", {
  withr::local_seed(1006)
  q <- 5
  fixture <- make_fit_fixture(S = 3, q = q, P = 6, T = 20)
  fit <- dkge_fit(fixture$betas, fixture$designs, K = fixture$K, rank = q - 1,
                  w_method = "none")

  expect_lte(fit$rank, q - 1)
  UtKU <- t(fit$U) %*% fit$K %*% fit$U
  expect_lt(max(abs(UtKU - diag(fit$rank))), 1e-8)
})

test_that("K-orthonormality holds with ridge > 0", {
  withr::local_seed(1007)
  fixture <- make_fit_fixture(S = 3, q = 4, P = 5, T = 20)
  fit <- dkge_fit(fixture$betas, fixture$designs, K = fixture$K, rank = 3,
                  ridge = 0.5, w_method = "none")

  UtKU <- t(fit$U) %*% fit$K %*% fit$U
  expect_lt(max(abs(UtKU - diag(fit$rank))), 1e-8)
})

test_that("K-orthonormality holds with MFA subject weights", {
  withr::local_seed(1008)
  fixture <- make_fit_fixture(S = 4, q = 4, P = 5, T = 20)
  fit <- dkge_fit(fixture$betas, fixture$designs, K = fixture$K, rank = 3,
                  w_method = "mfa_sigma1", w_tau = 0.3)

  UtKU <- t(fit$U) %*% fit$K %*% fit$U
  expect_lt(max(abs(UtKU - diag(fit$rank))), 1e-8)
})

test_that("K-orthonormality holds with energy subject weights", {
  withr::local_seed(1009)
  fixture <- make_fit_fixture(S = 4, q = 4, P = 5, T = 20)
  fit <- dkge_fit(fixture$betas, fixture$designs, K = fixture$K, rank = 3,
                  w_method = "energy", w_tau = 0.5)

  UtKU <- t(fit$U) %*% fit$K %*% fit$U
  expect_lt(max(abs(UtKU - diag(fit$rank))), 1e-8)
})
