# tests/testthat/test_dkge_analytic.R
context("Analytic LOSO DKGE: accuracy, fallback, invariance")

library(testthat)

# ---------- helpers (local to tests) ----------
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

# Build a small, consistent dkge "fit" object sufficient for dkge_analytic_loso() & dkge_loso_contrast()
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
    right <- Btil[[s]] %*% t(Btil[[s]])            # q×q
    S_s <- Khalf %*% right %*% Khalf               # q×q (PSD)
    contribs[[s]] <- (S_s + t(S_s)) / 2
    Chat <- Chat + weights[s] * contribs[[s]]
  }
  Chat <- (Chat + t(Chat)) / 2

  # Full eigen-decomposition of Chat (V in K^{1/2} domain), U = K^{-1/2} V
  eg <- eigen(Chat, symmetric = TRUE)
  V_full <- eg$vectors
  lambda_full <- eg$values
  stopifnot(ncol(V_full) == q, length(lambda_full) == q)

  U_hat <- V_full[, seq_len(r), drop = FALSE]      # q×r (Euclidean-orthonormal)
  U <- Kihalf %*% U_hat                             # q×r (K-orthonormal)

  # Shared ruler R (identity is fine for contrasts in tests)
  R <- diag(q)

  fit <- list(
    class            = "dkge",  # cosmetic; we'll set class below
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
    rank             = r
  )
  class(fit) <- c("dkge", "list")
  fit
}

# Norm helpers
rel_err <- function(x, y) {
  dx <- sqrt(sum((x - y)^2))
  dy <- sqrt(sum(y^2))
  if (dy < 1e-12) dx else dx / dy
}

max_abs <- function(M) max(abs(M))

# ---------- 1) Accuracy under low-leverage perturbation ----------
test_that("Analytic LOSO matches exact LOSO for low-leverage subject", {
  skip_on_cran()

  fit <- .make_fit(seed = 42, q = 9, r = 3, S = 12, P = 10,
                   low_leverage_subject = 2, low_scale = 5e-4)

  # Random contrast in design space
  set.seed(99)
  cvec <- rnorm(nrow(fit$U))

  # Exact LOSO (ground truth) vs Analytic
  exact   <- dkge_loso_contrast(fit, s = 2, c = cvec, ridge = 0)
  approx  <- dkge_analytic_loso(fit, s = 2, c = cvec, tol = 1e-8, fallback = TRUE)

  # Should use analytic path (not fallback) in this benign regime
  expect_identical(approx$method, "analytic")

  # Value agreement (subject contrast vector)
  expect_lt(rel_err(approx$v, exact$v), 2e-3)

  # Eigenvalue first-order update vs exact LOSO eigenvalues (top r)
  lam_exact  <- exact$evals[seq_len(fit$rank)]
  lam_approx <- approx$evals
  expect_lt(max(abs(lam_approx - lam_exact)), 1e-3)

  # K-orthonormality of returned held-out basis
  G <- t(approx$basis) %*% fit$K %*% approx$basis
  expect_lt(max_abs(G - diag(fit$rank)), 5e-8)
})

# ---------- 2) Fallback under near-degeneracy ----------
test_that("Analytic LOSO falls back when eigengap is tiny (near-degenerate spectrum)", {
  skip_on_cran()

  fit <- .make_fit(seed = 7, q = 8, r = 3, S = 8, P = 9,
                   low_leverage_subject = 1, low_scale = 1e-3)

  # Force a tiny eigengap in the stored full spectrum to trigger fallback
  lam <- fit$eig_values_full
  lam[2] <- lam[1] - 1e-10            # |gap| < tol -> fallback
  fit$eig_values_full <- lam

  # Keep an orthonormal V (arbitrary; only used to compute H before fallback)
  # No need to make V consistent with Chat since we expect immediate fallback.
  Q <- qr.Q(qr(matrix(rnorm(nrow(fit$U) * nrow(fit$U)), nrow(fit$U))))
  fit$eig_vectors_full <- Q

  cvec <- rnorm(nrow(fit$U))

  exact  <- dkge_loso_contrast(fit, s = 1, c = cvec, ridge = 0)
  approx <- dkge_analytic_loso(fit, s = 1, c = cvec, tol = 1e-6, fallback = TRUE)

  # Must take fallback path
  expect_identical(approx$method, "fallback")

  # And produce the same numbers as exact LOSO (since it calls dkge_loso_contrast)
  expect_lt(rel_err(approx$v, exact$v), 1e-12)     # essentially identical
  expect_lt(max_abs(approx$basis - exact$basis), 1e-10)
})

# ---------- 3) Basis permutation/rotation invariance (values) ----------
test_that("Analytic LOSO values are invariant to permutation of stored eigenbasis", {
  skip_on_cran()

  fit <- .make_fit(seed = 123, q = 10, r = 4, S = 10, P = 11,
                   low_leverage_subject = 3, low_scale = 1e-3)

  cvec <- rnorm(nrow(fit$U))

  # Original analytic LOSO
  a1 <- dkge_analytic_loso(fit, s = 3, c = cvec, tol = 1e-8, fallback = TRUE)
  expect_identical(a1$method, "analytic")

  # Permute the full eigenbasis (columns) and eigenvalues consistently
  perm <- c(2,1,3,4,5,6,7,8,9,10)      # simple swap of first two
  fit2 <- fit
  fit2$eig_vectors_full <- fit$eig_vectors_full[, perm, drop = FALSE]
  fit2$eig_values_full  <- fit$eig_values_full[perm]

  a2 <- dkge_analytic_loso(fit2, s = 3, c = cvec, tol = 1e-8, fallback = TRUE)
  expect_identical(a2$method, "analytic")

  # The *projected values* v must be invariant to any orthonormal relabeling
  expect_lt(rel_err(a2$v, a1$v), 1e-10)

  # And the returned basis remains K-orthonormal
  G2 <- t(a2$basis) %*% fit$K %*% a2$basis
  expect_lt(max_abs(G2 - diag(fit$rank)), 5e-8)
})