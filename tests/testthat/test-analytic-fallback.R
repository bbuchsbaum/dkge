# tests/testthat/test-analytic-fallback.R
# Tests for analytic LOSO fallback paths - complete coverage of all safety conditions
context("Analytic LOSO fallback paths")

library(testthat)

# ---------- helpers (local to tests) ----------

#' Symmetric eigenvalue square root
sym_eig_sqrt <- function(M, inv = FALSE, jitter = 1e-10) {
  M <- (M + t(M)) / 2
  ee <- eigen(M, symmetric = TRUE)
  vals <- pmax(ee$values, jitter)
  V <- ee$vectors
  if (!inv) {
    V %*% (diag(sqrt(vals), length(vals))) %*% t(V)
  } else {
    V %*% (diag(1 / sqrt(vals), length(vals))) %*% t(V)
  }
}

#' Build a small, consistent dkge "fit" object for testing
#' @param seed Random seed
#' @param q Number of design effects
#' @param r Rank
#' @param S Number of subjects
#' @param P Number of clusters
#' @param low_leverage_subject Subject index to have low leverage
#' @param low_scale Scale factor for low-leverage subject
.make_fit <- function(seed = 1L,
                      q = 8L, r = 3L, S = 10L, P = 12L,
                      low_leverage_subject = 1L,
                      low_scale = 1e-3) {
  set.seed(seed)

  # SPD design kernel
  A <- matrix(rnorm(q * q), q, q)
  K <- crossprod(A) + diag(q) * 0.5  # well-conditioned
  Khalf  <- sym_eig_sqrt(K, inv = FALSE)
  Kihalf <- sym_eig_sqrt(K, inv = TRUE)

  # Subject betas (already "row-standardized" in tests)
  Btil <- lapply(seq_len(S), function(s) {
    B <- matrix(rnorm(q * P), q, P)
    if (s == low_leverage_subject) B <- B * low_scale
    B
  })

  # Subject weights: make the low-leverage subject small
  weights <- rep(1, S)
  weights[low_leverage_subject] <- 0.05

  # Contributions and Chat (consistent with how dkge_fit builds them)
  contribs <- vector("list", S)
  Chat <- matrix(0, q, q)

  for (s in seq_len(S)) {
    right <- Btil[[s]] %*% t(Btil[[s]])            # q x q
    S_s <- Khalf %*% right %*% Khalf               # q x q (PSD)
    contribs[[s]] <- (S_s + t(S_s)) / 2
    Chat <- Chat + weights[s] * contribs[[s]]
  }
  Chat <- (Chat + t(Chat)) / 2

  # Full eigen-decomposition of Chat (V in K^{1/2} domain), U = K^{-1/2} V
  eg <- eigen(Chat, symmetric = TRUE)
  V_full <- eg$vectors
  lambda_full <- eg$values
  stopifnot(ncol(V_full) == q, length(lambda_full) == q)

  U_hat <- V_full[, seq_len(r), drop = FALSE]      # q x r (Euclidean-orthonormal)
  U <- Kihalf %*% U_hat                             # q x r (K-orthonormal)

  # Shared ruler R (identity is fine for contrasts in tests)
  R <- diag(q)

  fit <- list(
    class            = "dkge",
    U                = U,
    sdev             = sqrt(pmax(lambda_full[seq_len(r)], 0)),
    evals            = lambda_full,
    eig_vectors_full = V_full,
    eig_values_full  = lambda_full,
    K                = K,
    Khalf            = Khalf,
    Kihalf           = Kihalf,
    R                = R,
    Chat             = Chat,
    contribs         = contribs,
    weights          = weights,
    Btil             = Btil,
    subject_ids      = paste0("sub", seq_len(S)),
    rank             = r,
    solver           = "pooled"  # Default to pooled solver
  )
  class(fit) <- c("dkge", "list")
  fit
}

#' Recompute subject contribution S_s from Btil_s and Khalf
.recompute_contrib <- function(Btil_s, Khalf) {
  right <- Btil_s %*% t(Btil_s)
  S_s <- Khalf %*% right %*% Khalf
  (S_s + t(S_s)) / 2
}

# ---------- Test 1: Fallback when solver is not pooled ----------
test_that("Analytic LOSO falls back when solver is not pooled", {
  skip_on_cran()
  fit <- .make_fit(seed = 101, q = 8, r = 3, S = 8, P = 10)
  fit$solver <- "jd"  # Non-pooled solver

  cvec <- rnorm(nrow(fit$U))
  result <- dkge_analytic_loso(fit, s = 1, c = cvec)

  expect_equal(result$method, "fallback")
  expect_equal(result$diagnostic$reason, "solver_not_pooled")

  # Verify fallback produces same result as direct LOSO
  exact <- dkge_loso_contrast(fit, s = 1, c = cvec)
  expect_lt(max(abs(result$v - exact$v)), 1e-12)
})

# ---------- Test 2: Fallback when voxel_weights are nonuniform ----------
test_that("Analytic LOSO falls back when voxel_weights are nonuniform", {
  skip_on_cran()
  fit <- .make_fit(seed = 102, q = 8, r = 3, S = 8, P = 10)
  fit$voxel_weights <- c(1.0, 0.5, 1.0, 0.8, 1.2, 0.9, 1.1, 1.0)  # Nonuniform

  cvec <- rnorm(nrow(fit$U))
  result <- dkge_analytic_loso(fit, s = 1, c = cvec)

  expect_equal(result$method, "fallback")
  expect_equal(result$diagnostic$reason, "nonuniform_voxel_weights")
})

# ---------- Test 3: No fallback when voxel_weights are uniform ----------
test_that("Analytic LOSO does NOT fall back when voxel_weights are uniform", {
  skip_on_cran()
  fit <- .make_fit(seed = 103, q = 8, r = 3, S = 10, P = 10,
                   low_leverage_subject = 1, low_scale = 1e-3)
  fit$voxel_weights <- rep(1, 8)  # Uniform weights

  cvec <- rnorm(nrow(fit$U))
  result <- dkge_analytic_loso(fit, s = 1, c = cvec, tol = 1e-8)

  expect_equal(result$method, "analytic")
})
