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

  # Positive raw subject scores and their full-cohort normalization. LOSO must
  # re-normalize these scores over the training cohort rather than subtract a
  # full-cohort weighted contribution.
  weight_scores_raw <- rexp(S) + 0.1
  weights <- weight_scores_raw / mean(weight_scores_raw)

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
    subject_weight_scores_raw = weight_scores_raw,
    subject_weight_usable = rep(TRUE, S),
    w_method         = "fixture_raw",
    w_tau            = 0,
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

    out <- dkge_loso_contrast(fit, s = s, contrasts = cvec, ridge = 0)
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
  out1 <- dkge_loso_contrast(fit, s = s, contrasts = cvec, ridge = 0)

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
  out2 <- dkge_loso_contrast(fit_extreme, s = s, contrasts = cvec, ridge = 0)

  # Key test: out1$basis should equal out2$basis because subject s is excluded
  # from both held-out basis computations
  expect_lt(.max_abs(out1$basis - out2$basis), 1e-10,
            label = "LOSO basis should be identical regardless of held-out subject's values")

  # Also check eigenvalues are the same
  expect_lt(.max_abs(out1$evals - out2$evals), 1e-10,
            label = "LOSO eigenvalues should be identical regardless of held-out subject's values")

  # Verify the fold-repooled matrices directly. Subtracting a full-cohort
  # contribution is no longer the LOSO contract because raw subject scores are
  # re-normalized over the training cohort.
  train <- setdiff(seq_along(fit$Btil), s)
  Chat_minus1 <- dkge:::.dkge_fold_weight_context(fit, train)$Chat
  Chat_minus2 <- dkge:::.dkge_fold_weight_context(fit_extreme, train)$Chat
  expect_equal(Chat_minus1, Chat_minus2, tolerance = 1e-12,
               label = "Chat_minus should be identical when held-out subject is excluded")
})

# ========================== Test 3: Different subjects get different held-out bases ==========================

test_that("LOSO: different subjects get different held-out bases", {
  skip_on_cran()

  fit <- .make_fit(seed = 99, q = 10, r = 4, S = 6, P = 14)
  q <- nrow(fit$U)
  cvec <- rnorm(q)

  # Get LOSO bases for subjects 1, 2, 3
  out1 <- dkge_loso_contrast(fit, s = 1L, contrasts = cvec, ridge = 0)
  out2 <- dkge_loso_contrast(fit, s = 2L, contrasts = cvec, ridge = 0)
  out3 <- dkge_loso_contrast(fit, s = 3L, contrasts = cvec, ridge = 0)

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

# ========================== Test 4: K-fold with k=S equals LOSO exactly ==========================

test_that("KFOLD: k=S produces identical per-subject values to LOSO", {
  skip_on_cran()

  set.seed(101)

  # Use dkge_sim_toy to create data with known ground truth
  # S=5 subjects (more than existing S=3 test)
  factors <- list(A = list(L = 2), B = list(L = 2))
  sim <- dkge_sim_toy(factors, active_terms = c("A", "B"), S = 5, P = 12, snr = 9)

  fit <- dkge_fit(sim$B_list,
                  designs = sim$X_list,
                  K = sim$K,
                  w_method = "none",
                  rank = 2)

  # Use the true contrast
  c_true <- as.numeric(sim$U_true[, 1])

  # Run LOSO
  loso <- dkge_contrast(fit, c_true, method = "loso", align = FALSE)

  # Run K-fold with k=S (should equal LOSO)
  S <- length(sim$B_list)
  kf <- dkge_contrast(fit, c_true, method = "kfold", folds = S, align = FALSE)

  # Compare per-subject values with tighter tolerance (1e-10)
  for (s in seq_len(S)) {
    subj_id <- fit$subject_ids[s]
    v_loso <- loso$values[[1]][[subj_id]]
    v_kf <- kf$values[[1]][[subj_id]]

    rmse <- sqrt(mean((v_loso - v_kf)^2))
    expect_lt(rmse, 1e-10,
              label = sprintf("Subject %s: K-fold values should match LOSO", subj_id))
  }
})

# ========================== Test 5: K-fold vs LOSO basis comparison ==========================

test_that("KFOLD: held-out bases match between k=S K-fold and LOSO", {
  skip_on_cran()

  set.seed(202)

  # Smaller test to examine bases directly
  factors <- list(A = list(L = 2))
  sim <- dkge_sim_toy(factors, active_terms = "A", S = 4, P = 10, snr = 8)

  fit <- dkge_fit(sim$B_list,
                  designs = sim$X_list,
                  K = sim$K,
                  w_method = "none",
                  rank = 1)

  c_true <- as.numeric(sim$U_true[, 1])

  # Run both methods without alignment
  loso <- dkge_contrast(fit, c_true, method = "loso", align = FALSE)
  S <- length(sim$B_list)
  kf <- dkge_contrast(fit, c_true, method = "kfold", folds = S, align = FALSE)

  # Extract bases from metadata
  loso_bases <- loso$metadata$bases
  kf_bases <- kf$metadata$fold_bases

  # Verify bases match (accounting for possible sign flips)
  for (s in seq_len(S)) {
    basis_loso <- loso_bases[[s]]
    basis_kf <- kf_bases[[s]]

    # Check dimensions match
    expect_equal(dim(basis_loso), dim(basis_kf),
                 label = sprintf("Subject %d: basis dimensions should match", s))

    # Use cosine similarity to account for sign flips
    # For rank-1, |cosine| should be exactly 1
    # For higher rank, use Frobenius norm after sign alignment
    if (ncol(basis_loso) == 1) {
      cosine <- abs(sum(basis_loso * basis_kf))
      expect_gt(cosine, 1 - 1e-10,
                label = sprintf("Subject %d: basis cosine should be ~1", s))
    } else {
      # For multi-column basis, check subspace agreement via SVD
      gram <- t(basis_loso) %*% fit$K %*% basis_kf
      svals <- svd(gram, nu = 0, nv = 0)$d
      # All singular values should be ~1 for matching subspaces
      expect_lt(max(abs(svals - 1)), 1e-8,
                label = sprintf("Subject %d: subspace singular values should be ~1", s))
    }
  }
})

# ========================== Test 6: LOSO includes all-but-one subject ==========================

test_that("LOSO: Chat_minus renormalizes raw weights over all-but-one subject", {
  skip_on_cran()

  fit <- .make_fit(seed = 303, q = 8, r = 3, S = 6, P = 10)
  q <- nrow(fit$U)

  # Test for multiple subjects
  for (s in c(1L, 3L, 6L)) {
    train <- setdiff(seq_along(fit$Btil), s)
    fold_weights <- fit$subject_weight_scores_raw[train]
    fold_weights <- fold_weights / mean(fold_weights)
    expect_equal(mean(fold_weights), 1, tolerance = 1e-14)

    # Verify Chat_minus can be reconstructed from remaining subjects
    Chat_from_remaining <- matrix(0, q, q)
    for (j in seq_along(train)) {
      Chat_from_remaining <- Chat_from_remaining +
        fold_weights[j] * fit$contribs[[train[j]]]
    }
    Chat_from_remaining <- (Chat_from_remaining + t(Chat_from_remaining)) / 2

    # Verify LOSO function uses correct Chat_minus by comparing eigenvalues
    cvec <- rnorm(q)
    out <- dkge_loso_contrast(fit, s = s, contrasts = cvec, ridge = 0)

    # Eigenvalues of manual Chat_minus should match returned evals
    manual_eig <- eigen(Chat_from_remaining, symmetric = TRUE)$values
    expect_lt(.max_abs(out$evals - manual_eig), 1e-10,
              label = sprintf("Subject %d: LOSO eigenvalues should match manual Chat_minus", s))
  }
})
