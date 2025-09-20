
# dkge-sinkhorn.R
# Entropic OT (Sinkhorn) transport from subject clusters to a medoid parcellation.

# Internal: pairwise squared distances between two point sets (rows), returns n×m.
.pairwise_sqdist <- function(A, B) {
  # A: n×d, B: m×d
  n <- nrow(A); m <- nrow(B)
  An <- rowSums(A*A); Bn <- rowSums(B*B)
  # ||a-b||^2 = ||a||^2 + ||b||^2 - 2 a b'
  outer(An, rep(1, m)) + outer(rep(1, n), Bn) - 2 * (A %*% t(B))
}

# Internal: build cost matrix between clusters using embedding + spatial + (optional) size
# Aemb_s, Aemb_ref: P_s×r and Q×r (unit rows)
# X_s, X_ref: P_s×3 and Q×3 centroids (mm)
.dkge_cost_matrix <- function(Aemb_s, Aemb_ref, X_s, X_ref,
                              lambda_emb = 1, lambda_spa = 0.5, sigma_mm = 15,
                              sizes_s = NULL, sizes_ref = NULL, lambda_size = 0) {
  C_emb <- .pairwise_sqdist(Aemb_s, Aemb_ref)                    # P×Q
  C_spa <- .pairwise_sqdist(X_s / sigma_mm, X_ref / sigma_mm)
  C <- lambda_emb * C_emb + lambda_spa * C_spa
  if (!is.null(sizes_s) && !is.null(sizes_ref) && lambda_size > 0) {
    C <- C + lambda_size * (.outer_log2(sizes_s, sizes_ref))
  }
  C
}

.outer_log2 <- function(a, b) {
  # (log a - log b)^2
  la <- log(pmax(a, 1)); lb <- log(pmax(b, 1))
  outer(la, lb, function(x,y) (x - y)^2)
}

# Sinkhorn-Knopp in log-stabilized form
# Returns transport plan T (P×Q) coupling mu and nu.
sinkhorn_plan <- function(C, mu, nu, epsilon = 0.05, max_iter = 200, tol = 1e-6) {
  stopifnot(all(mu >= 0), all(nu >= 0), abs(sum(mu) - sum(nu)) < 1e-6)
  K <- exp(-C / epsilon)                           # P×Q
  u <- rep(1 / nrow(K), nrow(K))
  v <- rep(1 / ncol(K), ncol(K))

  # iterate in standard scaling (fast in practice for well-scaled costs)
  for (it in seq_len(max_iter)) {
    u_old <- u; v_old <- v
    # Avoid division by zero by flooring denominators
    Ku  <- K %*% v;  Ku[Ku < 1e-300] <- 1e-300
    u   <- mu / Ku
    KtU <- t(K) %*% u; KtU[KtU < 1e-300] <- 1e-300
    v   <- nu / KtU
    # check marginal error
    if (it %% 5 == 0) {
      err1 <- max(abs(u * (K %*% v) - mu))
      err2 <- max(abs(v * (t(K) %*% u) - nu))
      if (max(err1, err2) < tol) break
    }
    # damping to improve stability
    if (it %% 50 == 0) {
      u <- 0.5 * (u + u_old)
      v <- 0.5 * (v + v_old)
    }
  }
  T <- (u * K) %*% diag(v, nrow = length(v))
  T
}

#' Transport cluster values to a medoid via Sinkhorn OT (mass-preserving)
#'
#' @param v_list list of vectors v_s (cluster values) length P_s
#' @param A_list list of P_s×r embeddings (cluster loadings) per subject
#' @param centroids list of P_s×3 matrices of MNI centroids
#' @param sizes list of length-P_s vectors (cluster sizes/masses); if NULL uses 1's
#' @param medoid index of medoid subject
#' @param lambda_emb, lambda_spa, sigma_mm cost weights and spatial scale
#' @param epsilon entropic regularization
#' @param max_iter, tol Sinkhorn parameters
#' @return list(value = Q-vector on medoid clusters, subj_values = S×Q matrix, plans = list of T_s)
#' @export
dkge_transport_to_medoid_sinkhorn <- function(v_list, A_list, centroids, sizes = NULL,
                                              medoid,
                                              lambda_emb = 1, lambda_spa = 0.5, sigma_mm = 15,
                                              epsilon = 0.05, max_iter = 200, tol = 1e-6) {
  S <- length(v_list)
  if (is.null(sizes)) sizes <- lapply(v_list, function(v) rep(1, length(v)))
  Ahat <- lapply(A_list, function(A) A / pmax(sqrt(rowSums(A^2)), 1e-8))
  Q <- nrow(A_list[[medoid]])
  Y <- matrix(NA_real_, S, Q)
  plans <- vector("list", S)

  sizes_ref <- sizes[[medoid]]; sizes_ref <- sizes_ref / sum(sizes_ref)

  for (s in seq_len(S)) {
    if (s == medoid) {
      Y[s,] <- v_list[[s]]
      plans[[s]] <- diag(sizes_ref, nrow = Q)       # identity coupling
      next
    }
    As <- Ahat[[s]]; Xs <- centroids[[s]]; szs <- sizes[[s]]; szs <- szs / sum(szs)
    Aref <- Ahat[[medoid]]; Xref <- centroids[[medoid]]

    C <- .dkge_cost_matrix(As, Aref, Xs, Xref, lambda_emb, lambda_spa, sigma_mm)
    Tplan <- sinkhorn_plan(C, mu = szs, nu = sizes_ref, epsilon = epsilon, max_iter = max_iter, tol = tol)
    plans[[s]] <- Tplan
    # push values to medoid clusters with mass-preserving plan
    num <- as.numeric(t(v_list[[s]]) %*% Tplan)      # length Q
    den <- pmax(colSums(Tplan), 1e-12)
    Y[s,] <- num / den
  }
  list(value = apply(Y, 2, stats::median, na.rm = TRUE),
       subj_values = Y,
       plans = plans)
}
