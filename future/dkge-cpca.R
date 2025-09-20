
# dkge-cpca.R
# Optional CPCA-style row-side "inside-span" filtering implemented entirely in q-space.
# The idea: given the compressed matrix  Ĉ = Σ_s (K^{1/2} B̃_s Ω_s^{1/2}) (⋯)^T (q×q),
# and a row subspace T (q×q0) that encodes "design-explainable" effects,
# build the (Euclidean) projector in the K^{1/2}-space and split:
#   P̂ = K^{1/2} T (T' K T)^{-1} T' K^{1/2}
#   Ĉ_design = P̂ Ĉ P̂,   Ĉ_resid = (I − P̂) Ĉ (I − P̂)
# Both parts stay inside span(Ĉ) and are column-orthogonal in the stacked data (CPCA property).

#' K-orthogonal projector onto span(T) in effect space
#' @param T q×q0 matrix whose columns span the row subspace of interest
#' @param K q×q PSD metric (design kernel)
#' @return list(P_K = T (T' K T)^{-1} T' K, P_hat = K^{1/2} T (T' K T)^{-1} T' K^{1/2})
#' @export
dkge_projector_K <- function(T, K) {
  stopifnot(is.matrix(T), is.matrix(K), nrow(T) == nrow(K))
  M <- crossprod(T, K %*% T)                 # q0×q0
  Minv <- solve((M + t(M))/2)                # symmetric solve
  P_K <- T %*% Minv %*% crossprod(T, K)      # K-selfadjoint projector: P_K' K = K P_K
  # To get a Euclidean-orthogonal projector in the K^{1/2}-space:
  # P̂ = K^{1/2} T (T' K T)^{-1} T' K^{1/2}
  eig <- eigen((K + t(K))/2, symmetric = TRUE)
  vals <- pmax(eig$values, 1e-12); V <- eig$vectors
  Khalf <- V %*% (diag(sqrt(vals), length(vals))) %*% t(V)
  Phat <- Khalf %*% T %*% Minv %*% t(T) %*% Khalf
  list(P_K = P_K, P_hat = Phat)
}

#' Split a compressed matrix Ĉ into design vs residual parts (CPCA inside-span, q-space)
#' @param Chat q×q compressed matrix in K^{1/2}-metric (i.e., Chat = K^{1/2} C K^{1/2})
#' @param T q×q0 row subspace basis
#' @param K q×q
#' @return list(Chat_design, Chat_resid, P_hat)
#' @export
dkge_cpca_split_chat <- function(Chat, T, K) {
  pr <- dkge_projector_K(T, K)
  P <- pr$P_hat
  Iq <- diag(nrow(Chat))
  Chat <- (Chat + t(Chat))/2
  C1 <- P %*% Chat %*% P
  C2 <- (Iq - P) %*% Chat %*% (Iq - P)
  list(Chat_design = (C1 + t(C1))/2, Chat_resid = (C2 + t(C2))/2, P_hat = P)
}

#' Fit DKGE bases on CPCA-filtered parts (design/residual)
#' @param fit a dkge object (from dkge_fit or dkge_fit_fast)
#' @param blocks integer vector of row indices to define T via selection, or
#' @param T optional q×q0 explicit basis (if provided, blocks ignored)
#' @param part "design" | "resid" | "both"
#' @param rank target rank (default = ncol(fit$U))
#' @param ridge small ridge in K-metric
#' @return list(U_design, U_resid, evals_design, evals_resid, P_hat)
#' @export
dkge_fit_cpca <- function(fit, blocks = NULL, T = NULL,
                          part = c("design","resid","both"),
                          rank = NULL, ridge = 0) {
  part <- match.arg(part)
  q <- nrow(fit$U)
  if (is.null(rank)) rank <- ncol(fit$U)
  if (is.null(T)) {
    stopifnot(!is.null(blocks), all(blocks >= 1), all(blocks <= q))
    T <- diag(1, q)[, unique(blocks), drop = FALSE]
  }
  sp <- dkge_cpca_split_chat(fit$Chat, T, fit$K)
  out <- list(P_hat = sp$P_hat)
  if (part %in% c("design","both")) {
    C1 <- sp$Chat_design; if (ridge > 0) C1 <- C1 + ridge * diag(q)
    eg1 <- eigen(C1, symmetric = TRUE)
    out$U_design <- fit$Kihalf %*% eg1$vectors[, seq_len(rank), drop = FALSE]
    out$evals_design <- eg1$values
  }
  if (part %in% c("resid","both")) {
    C2 <- sp$Chat_resid; if (ridge > 0) C2 <- C2 + ridge * diag(q)
    eg2 <- eigen(C2, symmetric = TRUE)
    out$U_resid <- fit$Kihalf %*% eg2$vectors[, seq_len(rank), drop = FALSE]
    out$evals_resid <- eg2$values
  }
  out
}
