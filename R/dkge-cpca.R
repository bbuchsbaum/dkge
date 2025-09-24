# dkge-cpca.R
# Optional CPCA-style row-side filtering implemented entirely in q-space.

#' K-orthogonal projector onto span(T) in effect space
#'
#' Builds both the K-selfadjoint projector operating in effect coordinates and
#' the Euclidean projector in the K^{1/2} metric. These projectors are useful
#' for CPCA-style splits of the compressed covariance.
#'
#' All computations take place in the metric induced by the design kernel
#' `K`. Choosing a non-identity kernel therefore changes the projector: smooth
#' kernels diffuse energy across correlated effects, whereas a diagonal kernel
#' keeps the split tied to the raw effect indices. When the design-aligned
#' subspace is not well represented by coordinate axes, supply `cpca_T` with
#' columns that already capture the kernel geometry so the CPCA basis respects
#' those relationships.
#'
#' @param T qxq0 matrix whose columns span the subspace of interest.
#' @param K qxq positive semi-definite design kernel.
#' @return List with entries `P_K` (effect-space projector) and `P_hat`
#'   (projector in the K^{1/2} metric).
#' @export
dkge_projector_K <- function(T, K) {
  stopifnot(is.matrix(T), is.matrix(K), nrow(T) == nrow(K))
  M <- crossprod(T, K %*% T)
  Msym <- (M + t(M)) / 2
  chol_M <- tryCatch(chol(Msym), error = function(e) NULL)
  Minv <- if (is.null(chol_M)) {
    solve(Msym)
  } else {
    chol2inv(chol_M)
  }
  P_K <- T %*% Minv %*% crossprod(T, K)

  eig <- eigen((K + t(K)) / 2, symmetric = TRUE)
  vals <- pmax(eig$values, 1e-12)
  V <- eig$vectors
  Khalf <- V %*% diag(sqrt(vals), length(vals)) %*% t(V)
  P_hat <- Khalf %*% T %*% Minv %*% t(T) %*% Khalf

  list(P_K = P_K, P_hat = P_hat)
}

#' Split compressed covariance into design/residual parts (CPCA inside-span)
#'
#' @param Chat qxq compressed covariance expressed in the K^{1/2} metric.
#' @param T qxq0 matrix describing the design-aligned subspace.
#' @param K qxq design kernel.
#' @return List with `Chat_design`, `Chat_resid`, and the projector `P_hat`.
#' @export
dkge_cpca_split_chat <- function(Chat, T, K) {
  pr <- dkge_projector_K(T, K)
  P <- pr$P_hat
  Iq <- diag(nrow(Chat))

  Chat <- (Chat + t(Chat)) / 2
  C1 <- P %*% Chat %*% P
  C2 <- (Iq - P) %*% Chat %*% (Iq - P)

  list(
    Chat_design = (C1 + t(C1)) / 2,
    Chat_resid = (C2 + t(C2)) / 2,
    P_hat = P
  )
}

#' Fit DKGE bases on CPCA-filtered components
#'
#' Computes separate DKGE bases on the CPCA design-aligned or residual parts of
#' the compressed covariance. Helpful for decomposing the latent space into
#' interpretable sub-structures.
#'
#' @param fit A `dkge` object from [dkge_fit()] or [dkge()].
#' @param blocks Optional integer vector of effect indices defining the subspace
#'   (used when `T` is not supplied).
#' @param T Optional explicit qxq0 basis matrix; overrides `blocks` when given.
#' @param part Which portion to return: "design", "resid", or "both".
#' @param rank Target rank for the returned bases (defaults to `ncol(fit$U)`).
#' @param ridge Optional ridge term applied before the eigen decompositions.
#' @return List containing the requested bases/eigenvalues and the projector.
#' @export
dkge_fit_cpca <- function(fit, blocks = NULL, T = NULL,
                          part = c("design", "resid", "both"),
                          rank = NULL, ridge = 0) {
  stopifnot(inherits(fit, "dkge"))
  part <- match.arg(part)

  q <- nrow(fit$U)
  if (is.null(rank)) rank <- ncol(fit$U)

  if (is.null(T)) {
    stopifnot(!is.null(blocks), length(blocks) >= 1, all(blocks >= 1), all(blocks <= q))
    T <- diag(1, q)[, unique(blocks), drop = FALSE]
  }

  split <- dkge_cpca_split_chat(fit$Chat, T, fit$K)
  out <- list(P_hat = split$P_hat)

  if (part %in% c("design", "both")) {
    C1 <- split$Chat_design
    if (ridge > 0) C1 <- C1 + ridge * diag(q)
    eg1 <- eigen((C1 + t(C1)) / 2, symmetric = TRUE)
    out$U_design <- fit$Kihalf %*% eg1$vectors[, seq_len(rank), drop = FALSE]
    out$evals_design <- eg1$values
  }

  if (part %in% c("resid", "both")) {
    C2 <- split$Chat_resid
    if (ridge > 0) C2 <- C2 + ridge * diag(q)
    eg2 <- eigen((C2 + t(C2)) / 2, symmetric = TRUE)
    out$U_resid <- fit$Kihalf %*% eg2$vectors[, seq_len(rank), drop = FALSE]
    out$evals_resid <- eg2$values
  }

  out
}

#' Fit DKGE with CPCA filtering
#'
#' Convenience wrapper around [dkge()] that enables CPCA filtering with a
#' slightly terser interface. Either `cpca_blocks` or `cpca_T` must be provided.
#' Additional arguments are forwarded to [dkge()].
#'
#' @inheritParams dkge
#' @inheritParams dkge_fit
#' @param ... Additional arguments passed to [dkge()].
#' @export
dkge_cpca_fit <- function(..., cpca_blocks = NULL, cpca_T = NULL,
                           cpca_part = c("design", "resid", "both"),
                           cpca_ridge = 0) {
  if (is.null(cpca_blocks) && is.null(cpca_T)) {
    stop("Provide either cpca_blocks or cpca_T for CPCA filtering.")
  }
  cpca_part <- match.arg(cpca_part)
  dkge(..., cpca_blocks = cpca_blocks, cpca_T = cpca_T,
       cpca_part = cpca_part, cpca_ridge = cpca_ridge)
}
