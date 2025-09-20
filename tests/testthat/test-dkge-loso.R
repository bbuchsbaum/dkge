# test-dkge-loso.R
context("dkge_loso_contrast: mathematical correctness, ridge behavior, permutation equivariance")

library(testthat)

# ------------------------- local helpers -------------------------
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

rel_err <- function(x, y) {
  dx <- sqrt(sum((x - y)^2))
  dy <- sqrt(sum(y^2))
  if (dy < 1e-12) dx else dx / dy
}

max_abs <- function(M) max(abs(M))

# Build a consistent synthetic dkge "fit" sufficient for dkge_loso_contrast()
.make_fit <- function(seed = 11L, q = 8L, r = 3L, S = 9L, P = 12L) {
  set.seed(seed)

  # SPD design kernel K
  A <- matrix(rnorm(q * q), q, q)
  K <- crossprod(A) + diag(q) * 0.5
  Khalf  <- sym_eig_sqrt(K, inv = FALSE)
  Kihalf <- sym_eig_sqrt(K, inv = TRUE)

  # Shared ruler R (Cholesky of a pooled Gram)
  Gpool <- crossprod(matrix(rnorm(5 * q), 5, q)) + diag(q) * 0.2
  R <- chol(Gpool)

  # Subject betas (already "row-standardized" for the test)
  Btil <- lapply(seq_len(S), function(s) matrix(rnorm(q * P), q, P))

  # Positive weights
  weights <- rexp(S) + 0.1

  # Per-subject contributions and pooled Chat
  contribs <- vector("list", S)
  Chat <- matrix(0, q, q)
  for (s in seq_len(S)) {
    right <- Btil[[s]] %*% t(Btil[[s]])     # q×q
    S_s <- Khalf %*% right %*% Khalf        # q×q (PSD)
    S_s <- (S_s + t(S_s)) / 2
    contribs[[s]] <- S_s
    Chat <- Chat + weights[s] * S_s
  }
  Chat <- (Chat + t(Chat)) / 2

  # Full eigendecomposition (for completeness)
  eg <- eigen(Chat, symmetric = TRUE)
  Vfull <- eg$vectors
  lamfull <- eg$values

  # Set a reference U with exactly r columns (only its ncol() is used by LOSO)
  Uref <- Kihalf %*% Vfull[, seq_len(r), drop = FALSE]

  fit <- list(
    U                = Uref,
    evals            = lamfull,
    eig_vectors_full = Vfull,
    eig_values_full  = lamfull,
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

# ------------------------- 1) Core math, K-orthonormality, R dependence -------------------------
test_that("LOSO: matches manual pipeline; basis is K-orthonormal; depends on R (ruler)", {
  skip_on_cran()

  fit <- .make_fit(seed = 3, q = 9, r = 3, S = 8, P = 10)
  q <- nrow(fit$U); r <- ncol(fit$U)
  s <- 2L

  # Non-trivial contrast and non-identity R
  set.seed(99)
  cvec <- rnorm(q)

  # Function output
  out <- dkge_loso_contrast(fit, s = s, c = cvec, ridge = 0)

  # ---- Manual recomputation (should match function) ----
  Chat_minus <- fit$Chat - fit$weights[s] * fit$contribs[[s]]
  Chat_minus <- (Chat_minus + t(Chat_minus)) / 2
  egm <- eigen(Chat_minus, symmetric = TRUE)
  Uminus <- fit$Kihalf %*% egm$vectors[, seq_len(r), drop = FALSE]

  c_tilde <- backsolve(fit$R, cvec, transpose = FALSE)
  alpha   <- t(Uminus) %*% fit$K %*% c_tilde

  Bts <- fit$Btil[[s]]
  A_s <- t(Bts) %*% fit$K %*% Uminus
  v_manual <- as.numeric(A_s %*% alpha)

  # Numerical agreement
  expect_lt(rel_err(out$v, v_manual), 1e-12)
  expect_lt(max_abs(out$basis - Uminus), 1e-10)
  expect_lt(rel_err(out$alpha, alpha), 1e-12)

  # K-orthonormality: Uminus^T K Uminus = I_r
  G <- t(out$basis) %*% fit$K %*% out$basis
  expect_lt(max_abs(G - diag(r)), 5e-10)

  # Sanity check that R matters: if we *ignore* R in the manual step, alpha shifts
  alpha_noR <- t(Uminus) %*% fit$K %*% cvec
  expect_gt(max_abs(alpha_noR - alpha), 1e-6)  # using the ruler does change coordinates
})

# ------------------------- 2) Ridge behavior: eigenvalue shift -------------------------
test_that("LOSO: ridge shifts eigenvalues by +ridge and preserves K-orthonormality", {
  skip_on_cran()

  fit <- .make_fit(seed = 10, q = 8, r = 4, S = 7, P = 9)
  q <- nrow(fit$U); r <- ncol(fit$U)
  s <- 3L
  cvec <- rnorm(q)

  out0 <- dkge_loso_contrast(fit, s = s, c = cvec, ridge = 0)
  ridge <- 0.125
  out1 <- dkge_loso_contrast(fit, s = s, c = cvec, ridge = ridge)

  # Independent check: eigenvalues of (Chat_minus + ridge * I)
  Chat_minus <- fit$Chat - fit$weights[s] * fit$contribs[[s]]
  Chat_minus <- (Chat_minus + t(Chat_minus)) / 2
  lam0 <- eigen(Chat_minus, symmetric = TRUE)$values
  lam1 <- eigen(Chat_minus + ridge * diag(q), symmetric = TRUE)$values

  # Function returns full spectrum in $evals; check top r agree with lam1
  expect_lt(max_abs(out1$evals[seq_len(q)] - lam1), 1e-10)

  # And lam1 ≈ lam0 + ridge
  expect_lt(max_abs((lam0 + ridge) - lam1), 1e-10)

  # K-orthonormality still holds
  G1 <- t(out1$basis) %*% fit$K %*% out1$basis
  expect_lt(max_abs(G1 - diag(r)), 5e-10)
})

# ------------------------- 3) Equivariance to permutation of subject's cluster order -------------------------
test_that("LOSO: permuting columns of B_s permutes v_s identically; basis & alpha unchanged", {
  skip_on_cran()

  fit <- .make_fit(seed = 21, q = 10, r = 3, S = 6, P = 11)
  q <- nrow(fit$U); r <- ncol(fit$U)
  s <- 4L
  cvec <- rnorm(q)

  # Original result
  out <- dkge_loso_contrast(fit, s = s, c = cvec, ridge = 0)

  # Permute columns (clusters) for subject s only
  P <- ncol(fit$Btil[[s]])
  perm <- sample(P)
  Pi   <- diag(P)[, perm, drop = FALSE]   # column-permutation matrix
  fit2 <- fit
  fit2$Btil[[s]] <- fit$Btil[[s]] %*% Pi

  # Note: contrib S_s = Khalf B B^T Khalf is invariant to column permutations,
  # so Chat_minus is unchanged; the basis and alpha must be identical.
  out2 <- dkge_loso_contrast(fit2, s = s, c = cvec, ridge = 0)

  # Basis & alpha equality
  expect_lt(max_abs(out2$basis - out$basis), 1e-10)
  expect_lt(rel_err(out2$alpha, out$alpha), 1e-12)

  # v is permuted in the same way (right-multiplication by Pi permutes columns => v permutes)
  expect_lt(max_abs(out2$v - as.numeric(t(Pi) %*% out$v)), 1e-12)
})