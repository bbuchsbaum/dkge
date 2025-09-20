
# dkge-stream.R (v0.5)
# Memory-light, two-pass streaming fit over subjects/runs.
# Pass 1: accumulate G_pool = Σ X'X to build the pooled design Cholesky factor R (R'R = G_pool).
# Pass 2: stream subjects again to update the tiny q×q matrix  Ĉ = K^{1/2} C K^{1/2}.
# Optionally cache per-subject contributions S_s to disk for LOSO without recomputation.

# A lightweight "loader" interface:
# - n(): number of subjects
# - X(s): returns T_s×q design
# - B(s): returns q×P_s GLM betas for subject s
# - Omega(s): returns length-P_s or P_s×P_s (optional; may be NULL)
# Helper to build a loader from in-memory lists is provided below.

#' Build a simple loader from in-memory lists
#' @export
dkge_make_loader_from_lists <- function(B_list, X_list, Omega_list = NULL) {
  stopifnot(length(B_list) == length(X_list))
  if (is.null(Omega_list)) Omega_list <- vector("list", length(B_list))
  list(
    n = function() length(B_list),
    X = function(s) X_list[[s]],
    B = function(s) B_list[[s]],
    Omega = function(s) Omega_list[[s]]
  )
}

#' Streaming DKGE fit (two-pass, memory-light)
#'
#' @param loader list with functions n(), X(s), B(s), Omega(s)
#' @param K q×q design kernel
#' @param rank target rank
#' @param ridge small ridge in K-metric
#' @param w_method "mfa_sigma1"|"energy"|"none"
#' @param w_tau shrinkage toward equal weights (0..1)
#' @param cache_dir optional directory to save per-subject contributions S_s (as RDS)
#' @param verbose logical
#' @return a dkge object with S_s stored either in-memory or as file paths in $contribs_paths
#' @export
dkge_fit_streamed <- function(loader, K, rank = 6, ridge = 0, w_method = "mfa_sigma1", w_tau = 0.3,
                              cache_dir = NULL, verbose = TRUE) {
  S <- loader$n()
  X1 <- loader$X(1)
  q  <- ncol(X1)

  # ---- Pass 1: pooled Gram G_pool ----
  if (verbose) message("Pass 1/2: accumulating pooled Gram...")
  G_pool <- matrix(0, q, q)
  for (s in seq_len(S)) {
    Xs <- loader$X(s)
    G_pool <- G_pool + crossprod(Xs)
  }
  diag(G_pool) <- diag(G_pool) + 1e-10
  R <- chol(G_pool)  # upper, R'R = G_pool

  # K roots
  Ksym <- (K + t(K))/2
  ee <- eigen(Ksym, symmetric = TRUE)
  vals <- pmax(ee$values, 1e-10); V <- ee$vectors
  Khalf  <- V %*% (diag(sqrt(vals), length(vals))) %*% t(V)
  Kihalf <- V %*% (diag(1/sqrt(vals), length(vals))) %*% t(V)

  # ---- Pass 2: accumulate tiny Ĉ and compute subject weights ----
  if (verbose) message("Pass 2/2: accumulating compressed covariance...")
  Chat <- matrix(0, q, q)
  weights <- rep(1, S)
  contribs <- vector("list", S)
  contribs_paths <- NULL
  if (!is.null(cache_dir)) {
    dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
    contribs_paths <- character(S)
  }

  # helper: subject weight
  compute_weight <- function(Btil, Ω, method) {
    if (method == "none") return(1)
    Xs <- if (is.null(Ω)) Khalf %*% Btil else {
      if (is.vector(Ω)) Khalf %*% (Btil * rep(sqrt(Ω), each = q))
      else Khalf %*% Btil %*% sqrt(Ω)
    }
    d <- svd(Xs, nu = 0, nv = 0)$d
    if (method == "mfa_sigma1") w <- 1 / (d[1]^2 + 1e-12)
    else if (method == "energy") w <- 1 / (sum(Xs*Xs) + 1e-12)
    else w <- 1
    w
  }

  for (s in seq_len(S)) {
    Bs <- loader$B(s)           # q×P_s
    Ωs <- try(loader$Omega(s), silent = TRUE); if (inherits(Ωs, "try-error")) Ωs <- NULL
    Btil <- t(R) %*% Bs

    weights[s] <- compute_weight(Btil, Ωs, w_method)
    # contribution S_s = Khalf * (Btil Ω Btil') * Khalf
    if (is.null(Ωs)) {
      right <- Btil %*% t(Btil)
    } else if (is.vector(Ωs)) {
      right <- (Btil * rep(Ωs, each = q)) %*% t(Btil)
    } else {
      right <- Btil %*% Ωs %*% t(Btil)
    }
    S_s <- Khalf %*% right %*% Khalf

    if (is.null(cache_dir)) {
      contribs[[s]] <- S_s
    } else {
      fp <- file.path(cache_dir, sprintf("dkge_contrib_s%04d.rds", s))
      saveRDS(S_s, fp)
      contribs_paths[s] <- fp
    }
    Chat <- Chat + weights[s] * S_s
  }

  # shrink toward equal weights, normalize mean=1
  if (w_method != "none") {
    weights <- (1 - w_tau) * 1 + w_tau * (weights / mean(weights))
  }

  if (ridge > 0) Chat <- Chat + ridge * diag(q)
  eg <- eigen((Chat + t(Chat))/2, symmetric = TRUE)
  U <- Kihalf %*% eg$vectors[, seq_len(rank), drop = FALSE]

  out <- list(U = U, evals = eg$values,
              R = R, K = K, Khalf = Khalf, Kihalf = Kihalf,
              Chat = Chat, weights = weights,
              contribs = if (is.null(cache_dir)) contribs else NULL,
              contribs_paths = contribs_paths,
              loader = loader)     # keep loader to fetch B on demand
  class(out) <- "dkge_stream"
  out
}

#' LOSO contrast using a streamed fit (reads per-subject contribution from disk if needed)
#' @export
dkge_loso_contrast_stream <- function(fit, s, c, ridge = 0) {
  stopifnot(inherits(fit, "dkge_stream"))
  q <- nrow(fit$U)
  Cs <- if (!is.null(fit$contribs)) fit$contribs[[s]] else readRDS(fit$contribs_paths[s])
  Chat_minus <- fit$Chat - fit$weights[s] * Cs
  if (ridge > 0) Chat_minus <- Chat_minus + ridge * diag(q)
  eg <- eigen((Chat_minus + t(Chat_minus))/2, symmetric = TRUE)
  r <- ncol(fit$U)
  Uminus <- fit$Kihalf %*% eg$vectors[, seq_len(r), drop = FALSE]

  # Load subject's B and project on the fly
  Bs <- fit$loader$B(s)
  Btil_s <- t(fit$R) %*% Bs
  A_s <- t(Btil_s) %*% fit$K %*% Uminus

  c_tilde <- backsolve(fit$R, c, transpose = FALSE)
  alpha   <- t(Uminus) %*% fit$K %*% c_tilde
  v_s     <- as.numeric(A_s %*% alpha)

  list(v = v_s, alpha = alpha, basis = Uminus)
}
