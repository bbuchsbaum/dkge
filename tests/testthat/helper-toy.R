# Helper utilities for weight-related tests

testthat::local_edition(3)

# ---- Toy kernel info (2 x 2 x 3) with identity mapping ----
toy_kernel_info <- function() {
  levels <- list(A = 2L, B = 2L, time = 3L)
  Qcell  <- prod(unlist(levels))
  map    <- diag(Qcell)
  K      <- diag(Qcell)
  list(
    K = K,
    map = map,
    info = list(levels = levels)
  )
}

# ---- Toy betas (subject list of Q x V matrices) ----
toy_betas <- function(nsub = 3L, Q = 12L, V = 8L, seed = 123) {
  set.seed(seed)
  lapply(seq_len(nsub), function(s) {
    M <- matrix(rnorm(Q * V, mean = 0, sd = 1), nrow = Q, ncol = V)
    if (s == 1) M[1:2, ] <- M[1:2, ] + 0.7
    if (s == 2) M[3:4, ] <- M[3:4, ] + 0.5
    M
  })
}

# ---- Minimal dkge-like fit for fold-weight tests ----
toy_fold_fit <- function(nsub = 3L, Q = 12L, V = 8L, seed = 123) {
  B_list <- toy_betas(nsub = nsub, Q = Q, V = V, seed = seed)
  r <- Q
  fit <- list(
    Btil = B_list,
    Omega = vector("list", nsub),
    weights = rep(1, nsub),
    K = diag(Q),
    Khalf = diag(Q),
    Kihalf = diag(Q),
    U = diag(Q)[, seq_len(r), drop = FALSE],
    subject_ids = paste0("sub", seq_len(nsub)),
    kernel_info = toy_kernel_info(),
    weight_spec = dkge_weights(adapt = "none"),
    folds_index = list(nsub),
    B_list = B_list
  )
  class(fit) <- c("dkge", "list")
  fit
}

# ---- Small real dkge fit (uses dkge()) for update tests ----
toy_real_fit <- function(nsub = 3L, Q = 4L, V = 6L, weights = NULL, seed = 321) {
  set.seed(seed)
  effects <- paste0("eff", seq_len(Q))
  betas <- lapply(seq_len(nsub), function(s) {
    matrix(rnorm(Q * V), nrow = Q, ncol = V,
           dimnames = list(effects, paste0("v", seq_len(V))))
  })
  designs <- lapply(seq_len(nsub), function(s) {
    diag(Q)
  })
  designs <- lapply(designs, function(X) {
    colnames(X) <- effects
    X
  })
  kernel <- diag(Q)
  dkge(betas, designs = designs, kernel = kernel, keep_inputs = TRUE, weights = weights)
}

# ---- Reference weighted covariance builder (column scaling) ----
reference_weighted_G <- function(B_list, w, ridge = 0, block = 5000L) {
  stopifnot(length(B_list) > 0)
  Q <- nrow(B_list[[1L]])
  G <- matrix(0, Q, Q)
  sw <- sum(w)
  for (B in B_list) {
    V <- ncol(B)
    for (j in seq(1L, V, by = block)) {
      idx <- j:min(V, j + block - 1L)
      Bw  <- sweep(B[, idx, drop = FALSE], 2L, sqrt(w[idx]), "*")
      G   <- G + Bw %*% t(Bw)
    }
  }
  G <- G / (sw * length(B_list))
  if (ridge > 0) {
    G <- G + (ridge * sum(diag(G)) / nrow(G)) * diag(nrow(G))
  }
  G
}
