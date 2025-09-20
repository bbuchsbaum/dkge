# test-dkge-kfold.R
context("dkge_kfold: folds, equivalence to LOSO at k=S, permutation equivariance")

library(testthat)

# ------------------------- shared local helpers -------------------------
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

# Minimal synthetic dkge "fit" with all fields used by .dkge_contrast_kfold()
.make_fit <- function(seed = 1L, q = 9L, r = 3L, S = 10L, P = 12L) {
  set.seed(seed)

  # SPD design kernel K and its roots
  A <- matrix(rnorm(q * q), q, q)
  K <- crossprod(A) + diag(q) * 0.5
  Khalf  <- sym_eig_sqrt(K, inv = FALSE)
  Kihalf <- sym_eig_sqrt(K, inv = TRUE)

  # pooled design Cholesky factor (R'R = pooled design Gram)
  Gpool <- crossprod(matrix(rnorm(5 * q), 5, q)) + diag(q) * 0.25
  R <- chol(Gpool)

  # subject betas (already row-standardized for the test)
  Btil <- lapply(seq_len(S), function(s) matrix(rnorm(q * P), q, P))

  # positive subject weights
  weights <- rexp(S) + 0.1

  # per-subject contributions and pooled Chat
  contribs <- vector("list", S)
  Chat <- matrix(0, q, q)
  for (s in seq_len(S)) {
    right <- Btil[[s]] %*% t(Btil[[s]])    # q×q
    S_s <- Khalf %*% right %*% Khalf       # q×q PSD
    S_s <- (S_s + t(S_s)) / 2
    contribs[[s]] <- S_s
    Chat <- Chat + weights[s] * S_s
  }
  Chat <- (Chat + t(Chat)) / 2

  # reference U only to define r (its columns aren't used by kfold directly)
  eg <- eigen(Chat, symmetric = TRUE)
  Uref <- Kihalf %*% eg$vectors[, seq_len(r), drop = FALSE]

  fit <- list(
    U        = Uref,           # q×r
    K        = K,
    Khalf    = Khalf,
    Kihalf   = Kihalf,
    R        = R,
    Chat     = Chat,
    contribs = contribs,       # list of q×q
    weights  = weights,        # length S
    Btil     = Btil,           # list of q×P
    subject_ids = paste0("sub", seq_len(S))
  )
  class(fit) <- c("dkge", "list")
  fit
}

# helper: find which fold contains subject s
.which_fold <- function(folds, s) {
  which(vapply(folds$assignments, function(ix) s %in% ix, logical(1)))
}

# ------------------------- 1) folds: reproducible, balanced, printable -------------------------
test_that("dkge_define_folds(subject): reproducible seeds, full coverage, near-balanced", {
  skip_on_cran()

  fit <- .make_fit(S = 11)
  k   <- 4

  f1 <- dkge_define_folds(fit, type = "subject", k = k, seed = 2024)
  f2 <- dkge_define_folds(fit, type = "subject", k = k, seed = 2024)
  f3 <- dkge_define_folds(fit, type = "subject", k = k, seed = 7)

  # same seed -> identical assignments
  expect_identical(f1$assignments, f2$assignments)
  # different seed -> likely different permutation
  expect_false(identical(f1$assignments, f3$assignments))

  # full coverage, no out-of-range indices
  all_idx <- sort(unique(unlist(f1$assignments)))
  expect_identical(all_idx, seq_len(length(fit$Btil)))

  # near-balanced fold sizes (differ by at most 1)
  sizes <- vapply(f1$assignments, length, integer(1))
  expect_lte(max(sizes) - min(sizes), 1L)

  # align flag stored; class and metadata present
  expect_true(inherits(f1, "dkge_folds"))
  expect_true(is.list(f1$metadata))
  expect_false(isTRUE(f1$align))   # default align=FALSE unless specified

  # printing shows header and fold count
  expect_output(print(f1), "DKGE Fold Definition")
  expect_output(print(f1), paste0("Folds: ", k))
})

# ------------------------- 2) k = S reproduces LOSO: values match; bases are K-orthonormal & subspace-equal -------------------------
test_that(".dkge_contrast_kfold: with k=S equals LOSO and bases are valid", {
  skip_on_cran()

  q <- 10; r <- 3; S <- 7; P <- 8
  fit <- .make_fit(seed = 42, q = q, r = r, S = S, P = P)

  # one random contrast with full length q
  set.seed(9)
  cvec <- rnorm(q)
  clist <- list(c1 = cvec)

  # make deterministic one-subject-per-fold assignment
  folds <- dkge_define_folds(fit, type = "subject", k = S, seed = 101)

  # K-fold cross-fitting
  kres <- .dkge_contrast_kfold(fit, clist, folds = folds, ridge = 0, parallel = FALSE, verbose = FALSE)

  # LOSO baseline (per subject)
  loso_vs <- vector("list", S)
  loso_Us <- vector("list", S)
  for (s in seq_len(S)) {
    out <- dkge_loso_contrast(fit, s, cvec, ridge = 0)
    loso_vs[[s]] <- out$v
    loso_Us[[s]] <- out$basis
  }

  # 2.a values match per subject
  for (s in seq_len(S)) {
    v_kfold <- kres$values$c1[[s]]
    expect_lt(rel_err(v_kfold, loso_vs[[s]]), 1e-12)
  }

  # 2.b bases: K-orthonormal and subspace equality with LOSO
  for (s in seq_len(S)) {
    fidx <- .which_fold(folds, s)
    U_fold <- kres$metadata$fold_bases[[fidx]]
    # K-orthonormality
    G <- t(U_fold) %*% fit$K %*% U_fold
    expect_lt(max_abs(G - diag(r)), 5e-10)

    # Subspace equality t(U_loso) K U_fold has singular values ~1
    U_loso <- loso_Us[[s]]
    X <- t(U_loso) %*% fit$K %*% U_fold
    sv <- svd(X)$d
    expect_gt(min(sv), 1 - 1e-10)
    expect_lt(max_abs((sv - 1)), 1e-10)
  }

  # 2.c alpha in metadata has expected dimensions and matches formula
  a_fold <- kres$metadata$fold_alphas[[1]]  # k × r
  expect_true(is.matrix(a_fold))
  expect_equal(ncol(a_fold), r)
  expect_equal(nrow(a_fold), S)

  # check first fold alpha explicitly
  c_tilde <- backsolve(fit$R, cvec, transpose = FALSE)
  U1 <- kres$metadata$fold_bases[[1]]
  alpha1 <- t(U1) %*% fit$K %*% c_tilde
  expect_lt(rel_err(a_fold[1, ], as.numeric(alpha1)), 1e-12)
})

# ------------------------- 3) Equivariance: permuting a held-out subject's cluster order -------------------------
test_that(".dkge_contrast_kfold: permuting B_s columns only permutes that subject's v_s; bases unchanged", {
  skip_on_cran()

  q <- 9; r <- 3; S <- 6; P <- 10
  fit <- .make_fit(seed = 7, q = q, r = r, S = S, P = P)
  cvec <- rnorm(q)
  clist <- list(c1 = cvec)

  # k = S folds with fixed seed
  folds <- dkge_define_folds(fit, type = "subject", k = S, seed = 777)

  # baseline
  base <- .dkge_contrast_kfold(fit, clist, folds = folds, ridge = 0, parallel = FALSE, verbose = FALSE)

  # pick a subject and permute its cluster order (columns of Btil)
  s0 <- 3L
  Pmat <- diag(P)[, sample(P), drop = FALSE]
  fit2 <- fit
  fit2$Btil[[s0]] <- fit$Btil[[s0]] %*% Pmat

  # run again
  alt <- .dkge_contrast_kfold(fit2, clist, folds = folds, ridge = 0, parallel = FALSE, verbose = FALSE)

  # find the fold index that holds out s0
  fidx <- .which_fold(folds, s0)

  # Bases are computed on training set (thus invariant to reordering columns in held-out B_s)
  expect_lt(max_abs(alt$metadata$fold_bases[[fidx]] - base$metadata$fold_bases[[fidx]]), 1e-12)

  # The held-out value vector for s0 is permuted accordingly
  v_base <- base$values$c1[[s0]]
  v_alt  <- alt$values$c1[[s0]]
  expect_lt(max_abs(v_alt - as.numeric(t(Pmat) %*% v_base)), 1e-12)

  # Other subjects unaffected
  others <- setdiff(seq_len(S), s0)
  for (s in others) {
    expect_lt(rel_err(alt$values$c1[[s]], base$values$c1[[s]]), 1e-12)
  }
})
