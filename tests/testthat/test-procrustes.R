# test-procrustes.R
# Focused diagnostics for K-Procrustes utilities

library(testthat)

make_procrustes_toy <- function(q = 4, r = 2, seed = 101) {
  set.seed(seed)
  K <- crossprod(matrix(rnorm(q * q), q, q)) + diag(q) * 0.1
  U_ref <- qr.Q(qr(matrix(rnorm(q * r), q, r)))
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
