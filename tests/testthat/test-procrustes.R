# test-procrustes.R
# Focused diagnostics for K-Procrustes utilities

library(testthat)

make_procrustes_toy <- function(q = 4, r = 2, seed = 101) {
  set.seed(seed)
  K <- crossprod(matrix(rnorm(q * q), q, q)) + diag(q) * 0.1
  U_ref <- dkge_k_orthonormalize(matrix(rnorm(q * r), q, r), K)
  list(K = K, U_ref = U_ref)
}

test_that("dkge_k_orthonormalize produces K-orthonormal basis", {
  toy <- make_procrustes_toy()
  W <- matrix(rnorm(nrow(toy$K) * 3), nrow(toy$K), 3)
  U <- dkge_k_orthonormalize(W, toy$K)
  gram <- crossprod(U, toy$K %*% U)
  expect_equal(gram, diag(ncol(U)), tolerance = 1e-8)
})

test_that("dkge_procrustes_K aligns basis to reference", {
  toy <- make_procrustes_toy()
  noise <- matrix(rnorm(length(toy$U_ref)) * 0.005, nrow(toy$U_ref))
  U_target <- dkge_k_orthonormalize(toy$U_ref + noise, toy$K)
  pr <- dkge_procrustes_K(toy$U_ref, U_target, toy$K, allow_reflection = FALSE)
  expect_equal(t(pr$U_aligned) %*% toy$K %*% pr$U_aligned, diag(ncol(toy$U_ref)), tolerance = 1e-6)
  M <- t(toy$U_ref) %*% toy$K %*% pr$U_aligned
  expect_gt(min(svd(M)$d), 0.99)
})

test_that("dkge_align_bases_K returns aligned set and scores", {
  toy <- make_procrustes_toy()
  bases <- lapply(1:3, function(i) {
    noise <- matrix(rnorm(length(toy$U_ref)) * 0.02, nrow(toy$U_ref))
    dkge_k_orthonormalize(toy$U_ref + noise, toy$K)
  })
  aligned <- dkge_align_bases_K(bases, toy$K, ref = toy$U_ref, allow_reflection = FALSE)
  expect_length(aligned$U_aligned, length(bases))
  expect_true(all(sapply(aligned$U_aligned, function(U) {
    P <- t(toy$U_ref) %*% toy$K %*% U
    min(svd(P)$d) > 0.99
  })))
  expect_equal(length(aligned$score), length(bases))
})

test_that("dkge_procrustes_K recovers the rotation and maximises overlap for large angles", {
  toy <- make_procrustes_toy(q = 6, r = 3, seed = 11)
  Uref <- dkge_k_orthonormalize(toy$U_ref, toy$K)
  set.seed(202)
  Rot <- qr.Q(qr(matrix(rnorm(9), 3, 3)))   # large orthogonal rotation in component space
  U <- Uref %*% Rot                          # still K-orthonormal
  pr <- dkge_procrustes_K(Uref, U, toy$K, allow_reflection = TRUE)

  # Aligned basis matches the reference in the K-metric (the transposed rotation
  # would anti-align and leave t(Uref) K U_aligned != I).
  expect_equal(t(Uref) %*% toy$K %*% pr$U_aligned, diag(3), tolerance = 1e-8)

  # Overlap trace attains the optimum = sum of singular values of C.
  C <- t(Uref) %*% toy$K %*% U
  expect_equal(sum(diag(C %*% pr$R)), sum(svd(C)$d), tolerance = 1e-8)

  # Alignment is not worse than the raw (un-aligned) input in the K-metric.
  d_raw <- sum((toy$K %*% (U - Uref))^2)
  d_aln <- sum((toy$K %*% (pr$U_aligned - Uref))^2)
  expect_lte(d_aln, d_raw + 1e-8)
})

test_that("dkge_procrustes_K without reflection returns a proper rotation", {
  toy <- make_procrustes_toy(q = 5, r = 3, seed = 21)
  Uref <- dkge_k_orthonormalize(toy$U_ref, toy$K)
  set.seed(7)
  Refl <- qr.Q(qr(matrix(rnorm(9), 3, 3)))
  if (det(Refl) > 0) Refl[, 1] <- -Refl[, 1]   # force a reflection (det = -1)
  U <- Uref %*% Refl
  pr <- dkge_procrustes_K(Uref, U, toy$K, allow_reflection = FALSE)
  expect_equal(det(pr$R), 1, tolerance = 1e-8)
  expect_equal(t(pr$R) %*% pr$R, diag(3), tolerance = 1e-8)
  C <- t(Uref) %*% toy$K %*% U
  expect_equal(pr$d, sum(diag(C %*% pr$R)), tolerance = 1e-10)
  expect_equal(pr$d, 1, tolerance = 1e-8)
  expect_equal(pr$unconstrained_d, 3, tolerance = 1e-8)
})

test_that("K-orthonormalization rejects indefinite and metric-rank-deficient inputs", {
  # A previously validated matrix cannot lend its verdict to a modified copy.
  valid <- diag(2)
  dkge_k_orthonormalize(diag(2), valid)
  valid[2, 2] <- -1
  expect_error(
    dkge_k_orthonormalize(diag(2), valid),
    "positive semidefinite"
  )

  K_psd <- diag(c(1, 1, 0))
  U <- dkge_k_orthonormalize(diag(3)[, 1:2, drop = FALSE], K_psd)
  expect_equal(crossprod(U, K_psd %*% U), diag(2), tolerance = 1e-12)
  expect_error(dkge_k_orthonormalize(diag(3), K_psd), "K metric")
})

test_that("Procrustes fails closed for bases outside its K-orthonormal contract", {
  K <- diag(3)
  expect_error(dkge_procrustes_K(diag(3) * 2, diag(3), K),
               "Uref.*K-orthonormal")
  expect_error(dkge_procrustes_K(diag(3), diag(3)[, 1:2], K),
               "conformable")
})

test_that("consensus is K-orthonormal and invariant to basis-list order", {
  K <- diag(5)
  Uref <- diag(5)[, 1:3, drop = FALSE]
  R1 <- diag(c(-1, 1, 1))
  theta <- pi / 3
  R2 <- matrix(c(cos(theta), -sin(theta), 0,
                 sin(theta),  cos(theta), 0,
                 0,           0,          1), 3, 3, byrow = TRUE)
  bases <- list(Uref, Uref %*% R1, Uref %*% R2)

  forward <- dkge_consensus_basis_K(bases, K)
  reverse <- dkge_consensus_basis_K(rev(bases), K)
  expect_true(forward$converged)
  expect_true(reverse$converged)
  expect_equal(crossprod(forward$U, K %*% forward$U), diag(3), tolerance = 1e-10)
  overlap <- svd(crossprod(forward$U, K %*% reverse$U), nu = 0, nv = 0)$d
  expect_equal(overlap, rep(1, 3), tolerance = 1e-10)
})

test_that("alignment and consensus validate references, weights, and controls", {
  K <- diag(3)
  bases <- list(diag(3), diag(3))
  expect_error(dkge_align_bases_K(bases, K, ref = 3), "valid index")
  expect_error(dkge_align_bases_K(bases, K, weights = 1), "weights")
  expect_error(dkge_consensus_basis_K(bases, K, weights = c(1, -1)), "weights")
  expect_error(dkge_consensus_basis_K(bases, K, max_iter = 0), "max_iter")
  expect_error(dkge_consensus_basis_K(bases, K, tol = Inf), "tol")
})
