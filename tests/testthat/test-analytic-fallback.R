# tests/testthat/test-analytic-fallback.R
# Tests for analytic LOSO fallback paths - complete coverage of all safety conditions

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
  result <- dkge_analytic_loso(fit, s = 1, contrasts = cvec)

  expect_equal(result$method, "fallback")
  expect_equal(result$diagnostic$reason, "solver_not_pooled")

  # Verify fallback produces same result as direct LOSO
  exact <- dkge_loso_contrast(fit, s = 1, contrasts = cvec)
  expect_lt(max(abs(result$v - exact$v)), 1e-12)
})

# ---------- Test 2: Fallback when voxel_weights are nonuniform ----------
test_that("Analytic LOSO falls back when voxel_weights are nonuniform", {
  skip_on_cran()
  fit <- .make_fit(seed = 102, q = 8, r = 3, S = 8, P = 10)
  fit$voxel_weights <- c(1.0, 0.5, 1.0, 0.8, 1.2, 0.9, 1.1, 1.0)  # Nonuniform

  cvec <- rnorm(nrow(fit$U))
  result <- dkge_analytic_loso(fit, s = 1, contrasts = cvec)

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
  result <- dkge_analytic_loso(fit, s = 1, contrasts = cvec, tol = 1e-8)

  expect_equal(result$method, "analytic")
})

# ---------- Test 4: Fallback when eig_vectors_full is NULL ----------
test_that("Analytic LOSO falls back when eig_vectors_full is NULL", {
  skip_on_cran()
  fit <- .make_fit(seed = 104, q = 8, r = 3, S = 8, P = 10)
  fit$eig_vectors_full <- NULL

  cvec <- rnorm(nrow(fit$U))
  result <- dkge_analytic_loso(fit, s = 1, contrasts = cvec)

  expect_equal(result$method, "fallback")
  expect_equal(result$diagnostic$reason, "missing_full_decomposition")
})

# ---------- Test 5: Fallback when eig_values_full is NULL ----------
test_that("Analytic LOSO falls back when eig_values_full is NULL", {
  skip_on_cran()
  fit <- .make_fit(seed = 105, q = 8, r = 3, S = 8, P = 10)
  fit$eig_values_full <- NULL

  cvec <- rnorm(nrow(fit$U))
  result <- dkge_analytic_loso(fit, s = 1, contrasts = cvec)

  expect_equal(result$method, "fallback")
  expect_equal(result$diagnostic$reason, "missing_full_decomposition")
})

# ---------- Test 6: Fallback when perturbation magnitude is large ----------
test_that("Analytic LOSO falls back when perturbation magnitude is large", {
  skip_on_cran()

  # This test triggers the perturbation_magnitude fallback condition.
  # The fallback triggers when: |coeffs[-j]| = |w_s * H[-j, j] / gaps[-j]| > 0.1
  # where H = V^T S_s V is the coupling matrix of subject contribution S_s.
  #
  # Strategy:
  # 1. Create a fit object with well-separated eigenvalues
  # 2. Manipulate stored eigenvalues to have very small gaps (but > eigengap tol)
  # 3. Construct subject 1 contribution S_s that has large off-diagonal coupling
  #    with respect to the Chat eigenvectors

  fit <- .make_fit(seed = 206, q = 6, r = 3, S = 4, P = 8)

  # Step 1: Get the V_full basis and create nearly-degenerate eigenvalues
  V_full <- fit$eig_vectors_full
  q <- nrow(V_full)

  # Create eigenvalues with tiny gaps (above eigengap tol but small)
  # Gap must be > 1e-6 (default eigengap tol) but small enough that w_s * H / gap > 0.1
  # Gap = 0.001, H = 1, w_s = 0.5 => coeff = 0.5 * 1 / 0.001 = 500 >> 0.1
  new_lambda <- c(100.0, 99.999, 99.998, 10, 5, 1)  # Top 3 very close
  fit$eig_values_full <- new_lambda
  fit$evals <- new_lambda

  # Step 2: Construct subject 1 contribution with large off-diagonal coupling
  # S_s = w w^T where w mixes eigenvectors v_1 and v_2
  w_vec <- V_full[, 1] + 2 * V_full[, 2]  # Linear combo of v_1 and v_2
  S_new <- w_vec %*% t(w_vec)
  S_new <- (S_new + t(S_new)) / 2  # Symmetrize

  fit$contribs[[1]] <- S_new
  fit$weights[1] <- 0.5  # Moderate weight

  # Verify setup: check that H has large off-diagonal
  H <- t(V_full) %*% S_new %*% V_full
  # H[1, 2] should be ~ 2 (from the mix of v_1 + 2*v_2)
  expect_gt(abs(H[1, 2]), 1.5)

  cvec <- rnorm(nrow(fit$U))
  result <- dkge_analytic_loso(fit, s = 1, contrasts = cvec, tol = 1e-8)

  # Should fall back due to perturbation magnitude (coeff >> 0.1)
  expect_equal(result$method, "fallback")
  expect_equal(result$diagnostic$reason, "perturbation_magnitude")

  # Verify fallback produces valid result (same as direct LOSO)
  exact <- dkge_loso_contrast(fit, s = 1, contrasts = cvec)
  expect_lt(max(abs(result$v - exact$v)), 1e-12)
})

# ---------- Test 7: Analytic matches iterative within tolerance ----------
test_that("Analytic matches iterative LOSO within tolerance for stable input", {
  skip_on_cran()

  # Low-leverage subject = stable regime for perturbation
  # This tests the canonical accuracy requirement: cosine(v, v_exact) > 0.98
  fit <- .make_fit(seed = 107, q = 10, r = 3, S = 12, P = 12,
                   low_leverage_subject = 1, low_scale = 1e-4)

  cvec <- rnorm(nrow(fit$U))

  exact <- dkge_loso_contrast(fit, s = 1, contrasts = cvec)
  approx <- dkge_analytic_loso(fit, s = 1, contrasts = cvec, tol = 1e-8)

  expect_equal(approx$method, "analytic")

  # Value comparison: cosine similarity > 0.98 (canonical threshold from Phase 3 criteria)
  cosine_v <- sum(approx$v * exact$v) / (sqrt(sum(approx$v^2)) * sqrt(sum(exact$v^2)))
  expect_gt(abs(cosine_v), 0.98)

  # Relative error on values (complementary check)
  rel_err_v <- sqrt(sum((approx$v - exact$v)^2)) / sqrt(sum(exact$v^2))
  expect_lt(rel_err_v, 0.01)  # < 1% relative error

  # Eigenvalue comparison: relative error < 1%
  approx_evals <- approx$evals
  exact_evals <- exact$evals[seq_len(length(approx_evals))]
  rel_err_lambda <- max(abs(approx_evals - exact_evals)) / max(abs(exact_evals))
  expect_lt(rel_err_lambda, 0.01)

  # Basis K-orthonormality check (basis is valid even if rotated)
  G <- t(approx$basis) %*% fit$K %*% approx$basis
  expect_lt(max(abs(G - diag(ncol(approx$basis)))), 1e-6)
})

# ---------- Test 8: Diagnostic information is populated correctly ----------
test_that("Analytic LOSO populates diagnostic information", {
  skip_on_cran()

  fit <- .make_fit(seed = 108, q = 8, r = 3, S = 10, P = 10,
                   low_leverage_subject = 1, low_scale = 1e-3)
  cvec <- rnorm(nrow(fit$U))

  result <- dkge_analytic_loso(fit, s = 1, contrasts = cvec, tol = 1e-8)

  expect_equal(result$method, "analytic")
  expect_true(!is.null(result$diagnostic))
  expect_equal(result$diagnostic$reason, "analytic")
  expect_true(is.numeric(result$diagnostic$min_eigengap))
  expect_true(is.numeric(result$diagnostic$max_abs_coeff))
  expect_true(is.numeric(result$diagnostic$threshold_eigengap))
  expect_true(is.numeric(result$diagnostic$threshold_coeff))
})
