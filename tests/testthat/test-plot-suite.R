library(testthat)
library(dkge)

make_plot_fit <- function(S = 3, P = NULL, seed = 1) {
  set.seed(seed)
  factors <- list(A = list(L = 2), B = list(L = 2), time = list(L = 4))
  dk <- design_kernel(factors, basis = "effect")
  q <- nrow(dk$K)
  labels <- rownames(dk$K)
  P <- P %||% (q + 5)
  betas <- replicate(S, {
    B <- matrix(rnorm(q * P), q, P)
    rownames(B) <- labels
    B
  }, simplify = FALSE)
  designs <- replicate(S, {
    X <- diag(q)
    dimnames(X) <- list(labels, labels)
    X
  }, simplify = FALSE)
  fit <- dkge_fit(betas, designs, K = dk, rank = min(3, q))
  list(fit = fit, kernel = dk)
}

test_that("plot suite falls back gracefully when patchwork missing", {
  skip_if_not_installed("ggplot2")
  fixture <- make_plot_fit()
  fit <- fixture$fit
  orig_requireNamespace <- base::requireNamespace
  plot_obj <- with_mocked_bindings(
    dkge_plot_suite(fit, save_path = NULL),
    .package = "base",
    requireNamespace = function(pkg, ...) {
      if (pkg == "patchwork") return(FALSE)
      orig_requireNamespace(pkg, ...)
    }
  )
  expect_s3_class(plot_obj, "ggplot")
})

test_that("dkge_principal_angles_K computes K-metric basis angles", {
  K <- diag(3)
  U1 <- diag(3)[, 1:2]
  U2 <- U1

  angles <- dkge_principal_angles_K(U1, U2, K)

  expect_equal(length(angles), 2L)
  expect_equal(angles, c(0, 0), tolerance = 1e-10)
})

test_that("principal angles depend on the subspace, not on the basis scaling", {
  K <- diag(c(4, 1, 1))
  U1 <- diag(3)[, 1:2]

  # Any other basis of the same subspace must give zero angles.
  expect_equal(dkge_principal_angles_K(U1, 0.5 * U1, K), c(0, 0), tolerance = 1e-6)
  mixer <- matrix(c(1, 2, -1, 3), 2, 2)
  expect_equal(dkge_principal_angles_K(U1, U1 %*% mixer, K), c(0, 0), tolerance = 1e-6)

  # K-orthogonal directions are 90 degrees apart.
  e1 <- matrix(c(1, 0, 0), ncol = 1)
  e2 <- matrix(c(0, 1, 0), ncol = 1)
  expect_equal(dkge_principal_angles_K(e1, e2, K), pi / 2, tolerance = 1e-10)
})

test_that("principal angles respect the numerical rank of each basis", {
  K <- diag(c(2, 1, 1))
  e12 <- diag(3)[, 1:2]
  # A rank-one subspace handed over with a duplicated column. Padding it out to
  # two columns with an arbitrary complement (as a plain QR does) would invent a
  # direction shared with span{e1, e2} and report angles (0, 90).
  e3_dup <- cbind(c(0, 0, 1), c(0, 0, 1))

  angles <- dkge_principal_angles_K(e12, e3_dup, K)
  expect_equal(length(angles), 1L)
  expect_equal(angles, pi / 2, tolerance = 1e-8)
  expect_equal(dkge_principal_angles_K(e3_dup, e12, K), pi / 2, tolerance = 1e-8)

  # A rank-deficient basis against itself: zeros, one per dimension it spans.
  self <- dkge_principal_angles_K(e3_dup, e3_dup, K)
  expect_equal(length(self), 1L)
  expect_equal(self, 0, tolerance = 1e-6)

  rank2_dup <- cbind(e12, e12[, 1])
  self2 <- dkge_principal_angles_K(rank2_dup, rank2_dup, K)
  expect_equal(length(self2), 2L)
  expect_equal(self2, c(0, 0), tolerance = 1e-6)

  # Full-rank behaviour is unchanged.
  expect_equal(length(dkge_principal_angles_K(e12, e12, K)), 2L)

  # A basis spanning nothing is an error rather than a fabricated subspace.
  expect_error(dkge_principal_angles_K(matrix(0, 3, 2), e12, K),
               "spans nothing")
})

test_that("subspace stability plots survive rank-deficient bases", {
  skip_if_not_installed("ggplot2")
  K <- diag(3)
  full <- diag(3)[, 1:2]
  deficient <- cbind(c(1, 0, 0), c(1, 0, 0))
  p <- dkge_plot_subspace_stability(list(full, deficient), K, consensus = full)
  expect_s3_class(p, "ggplot")
  # Two angles for the full-rank basis, one for the rank-one basis.
  expect_equal(nrow(p$data), 3L)
})

test_that("principal angles under a non-identity K match a hand computation", {
  K <- diag(c(4, 1, 1))
  # Both bases are K-orthonormal: u^T K u = 1.
  u1 <- matrix(c(0.5, 0, 0), ncol = 1)
  u2 <- matrix(c(1 / (2 * sqrt(2)), 1 / sqrt(2), 0), ncol = 1)
  expect_equal(as.numeric(t(u1) %*% K %*% u1), 1)
  expect_equal(as.numeric(t(u2) %*% K %*% u2), 1)

  # cos(theta) = u1' K u2 = 1/sqrt(2)  =>  theta = 45 degrees.
  expect_equal(dkge_principal_angles_K(u1, u2, K), pi / 4, tolerance = 1e-8)
  expect_equal(dkge_principal_angles_K(u1, u2, K, orthonormalize = FALSE),
               pi / 4, tolerance = 1e-8)

  # A non-K-orthonormal representation of the same subspace gives the same angle.
  expect_equal(dkge_principal_angles_K(u1, 7 * u2, K), pi / 4, tolerance = 1e-8)
})
