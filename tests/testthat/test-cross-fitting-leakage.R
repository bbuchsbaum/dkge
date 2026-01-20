# test-cross-fitting-leakage.R
# Tests for data leakage prevention and cross-fitting mechanics

library(testthat)

# ========================== Local Helpers ==========================

#' Compute symmetric matrix square root (or inverse)
.sym_eig_sqrt <- function(M, inv = FALSE, jitter = 1e-10) {
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

#' Relative error between two vectors/matrices
.rel_err <- function(x, y) {
  dx <- sqrt(sum((x - y)^2))
  dy <- sqrt(sum(y^2))
  if (dy < 1e-12) dx else dx / dy
}

#' Maximum absolute value
.max_abs <- function(M) max(abs(M))

#' Build a synthetic dkge "fit" for testing LOSO leakage prevention
.make_fit <- function(seed = 11L, q = 8L, r = 3L, S = 9L, P = 12L) {
  set.seed(seed)

  # SPD design kernel K
  A <- matrix(rnorm(q * q), q, q)
  K <- crossprod(A) + diag(q) * 0.5
  Khalf  <- .sym_eig_sqrt(K, inv = FALSE)
  Kihalf <- .sym_eig_sqrt(K, inv = TRUE)

  # Shared ruler R (Cholesky of a pooled Gram)
  Gpool <- crossprod(matrix(rnorm(5 * q), 5, q)) + diag(q) * 0.2
  R <- chol(Gpool)

  # Subject betas
  Btil <- lapply(seq_len(S), function(s) matrix(rnorm(q * P), q, P))

  # Positive weights
  weights <- rexp(S) + 0.1

  # Per-subject contributions and pooled Chat
  contribs <- vector("list", S)
  Chat <- matrix(0, q, q)
  for (s in seq_len(S)) {
    right <- Btil[[s]] %*% t(Btil[[s]])     # q x q
    S_s <- Khalf %*% right %*% Khalf        # q x q (PSD)
    S_s <- (S_s + t(S_s)) / 2
    contribs[[s]] <- S_s
    Chat <- Chat + weights[s] * S_s
  }
  Chat <- (Chat + t(Chat)) / 2

  # Full eigendecomposition
  eg <- eigen(Chat, symmetric = TRUE)
  Vfull <- eg$vectors
  lamfull <- eg$values

  # Set a reference U with exactly r columns
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

#' Recompute subject contribution for extreme injection test
.recompute_contrib <- function(Btil_s, Khalf) {
  right <- Btil_s %*% t(Btil_s)
  S_s <- Khalf %*% right %*% Khalf
  (S_s + t(S_s)) / 2
}

# ========================== Test 1: LOSO basis differs from full basis ==========================

test_that("LOSO: held-out basis U_minus differs from full basis U", {
  skip_on_cran()

  fit <- .make_fit(seed = 42, q = 10, r = 4, S = 7, P = 15)
  q <- nrow(fit$U)
  r <- ncol(fit$U)

  # Test multiple subjects, not just s=1
  for (s in c(1L, 3L, 5L)) {
    cvec <- rnorm(q)

    out <- dkge_loso_contrast(fit, s = s, c = cvec, ridge = 0)
    U_minus <- out$basis
    U_full <- fit$U

    # U_minus should be different from U_full when subject has non-trivial contribution
    frob_diff <- norm(U_minus - U_full, "F")

    # Subject's contribution to Chat should be non-zero (non-trivial)
    contrib_frob <- norm(fit$contribs[[s]], "F")
    expect_gt(contrib_frob, 1e-6,
              label = sprintf("Subject %d should have non-trivial contribution", s))

    # Therefore U_minus should differ from U_full
    # Use a tolerance that accounts for subjects with smaller contribution
    expect_gt(frob_diff, 1e-8,
              label = sprintf("LOSO basis for subject %d should differ from full basis", s))
  }
})

# ========================== Test 2: Extreme value injection does not affect U_minus ==========================

test_that("LOSO: extreme values in held-out subject do not affect U_minus", {
  skip_on_cran()

  fit <- .make_fit(seed = 77, q = 9, r = 3, S = 8, P = 12)
  q <- nrow(fit$U)
  s <- 4L
  cvec <- rnorm(q)

  # Get baseline LOSO basis for subject s
  out1 <- dkge_loso_contrast(fit, s = s, c = cvec, ridge = 0)

  # Create modified fit with extreme values in subject s only
  fit_extreme <- fit
  fit_extreme$Btil[[s]] <- fit$Btil[[s]] * 1000  # Scale up by 1000x

  # Recompute contribution for subject s (Chat will also need updating)
  fit_extreme$contribs[[s]] <- .recompute_contrib(fit_extreme$Btil[[s]], fit$Khalf)

  # Update Chat to reflect the extreme contribution
  Chat_new <- matrix(0, q, q)
  for (i in seq_along(fit_extreme$Btil)) {
    Chat_new <- Chat_new + fit_extreme$weights[i] * fit_extreme$contribs[[i]]
  }
  fit_extreme$Chat <- (Chat_new + t(Chat_new)) / 2

  # Get LOSO basis for same subject with extreme values
  out2 <- dkge_loso_contrast(fit_extreme, s = s, c = cvec, ridge = 0)

  # Key test: out1$basis should equal out2$basis because subject s is excluded
  # from both held-out basis computations
  expect_lt(.max_abs(out1$basis - out2$basis), 1e-10,
            label = "LOSO basis should be identical regardless of held-out subject's values")

  # Also check eigenvalues are the same
  expect_lt(.max_abs(out1$evals - out2$evals), 1e-10,
            label = "LOSO eigenvalues should be identical regardless of held-out subject's values")

  # Verify Chat_minus is the same (computed without subject s)
  # Compute manually for verification
  # Note: use 1e-7 tolerance due to floating-point accumulation in Chat computation
  Chat_minus_manual1 <- fit$Chat - fit$weights[s] * fit$contribs[[s]]
  Chat_minus_manual2 <- fit_extreme$Chat - fit_extreme$weights[s] * fit_extreme$contribs[[s]]

  expect_lt(.max_abs(Chat_minus_manual1 - Chat_minus_manual2), 1e-7,
            label = "Chat_minus should be identical when held-out subject is excluded")
})

# ========================== Test 3: Different subjects get different held-out bases ==========================

test_that("LOSO: different subjects get different held-out bases", {
  skip_on_cran()

  fit <- .make_fit(seed = 99, q = 10, r = 4, S = 6, P = 14)
  q <- nrow(fit$U)
  cvec <- rnorm(q)

  # Get LOSO bases for subjects 1, 2, 3
  out1 <- dkge_loso_contrast(fit, s = 1L, c = cvec, ridge = 0)
  out2 <- dkge_loso_contrast(fit, s = 2L, c = cvec, ridge = 0)
  out3 <- dkge_loso_contrast(fit, s = 3L, c = cvec, ridge = 0)

  U1 <- out1$basis
  U2 <- out2$basis
  U3 <- out3$basis

  # All three bases should be different from each other
  # (unless subjects have nearly identical contributions, which is unlikely with random data)
  diff_12 <- norm(U1 - U2, "F")
  diff_13 <- norm(U1 - U3, "F")
  diff_23 <- norm(U2 - U3, "F")

  expect_gt(diff_12, 1e-8,
            label = "U_minus for subjects 1 and 2 should differ")
  expect_gt(diff_13, 1e-8,
            label = "U_minus for subjects 1 and 3 should differ")
  expect_gt(diff_23, 1e-8,
            label = "U_minus for subjects 2 and 3 should differ")
})
