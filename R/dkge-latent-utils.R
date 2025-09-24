#' DKGE latent-space utilities
#'
#' Helpers for moving between subject-level cluster representations and the
#' shared DKGE latent basis.
#'
#' @name dkge-latent-utils
NULL

#' Project subject clusters into the DKGE latent space
#'
#' For each subject, computes \eqn{Z_s = B_s^\top K U} so that every row
#' represents a subject cluster embedded in the \eqn{r}-dimensional DKGE
#' latent space. These projections are commonly used for training classifiers
#' or computing Haufe-style encoders.
#'
#' @param fit Fitted object of class `dkge` containing `Btil`, `K`, and `U`.
#' @return A list of length `S` (number of subjects). Element `s` is a
#'   `P_s x r` matrix holding the latent representation of subject `s`'s
#'   clusters.
#' @export
#' @examples
#' \dontrun{
#' Z_list <- dkge_project_clusters_to_latent(fit)
#' str(Z_list[[1]])  # P_1 x r matrix
#' }
dkge_project_clusters_to_latent <- function(fit) {
  stopifnot(inherits(fit, "dkge"))
  Bs_list <- fit$Btil
  U <- fit$U
  K <- fit$K
  stopifnot(is.list(Bs_list), !is.null(U), !is.null(K))
  stopifnot(ncol(U) >= 1L)

  KU <- K %*% U
  lapply(Bs_list, function(Bs) {
    stopifnot(is.matrix(Bs))
    proj <- crossprod(Bs, KU)  # P_s x r
    as.matrix(proj)
  })
}

#' Cluster-to-latent loadings for DKGE subjects
#'
#' Computes \eqn{A_s = B_s^\top K U}, the linear map that pulls latent-space
#' vectors back to subject cluster space. Multiplying \eqn{A_s} by a latent
#' coefficient vector produces a Haufe-style decoder at the subject level.
#'
#' @inheritParams dkge_project_clusters_to_latent
#' @return A list of `P_s x r` matrices of cluster loadings.
#' @export
dkge_cluster_loadings <- function(fit) {
  stopifnot(inherits(fit, "dkge"))
  Bs_list <- fit$Btil
  U <- fit$U
  K <- fit$K
  stopifnot(is.list(Bs_list), !is.null(U), !is.null(K))

  KU <- K %*% U
  lapply(Bs_list, function(Bs) {
    stopifnot(is.matrix(Bs))
    loadings <- crossprod(Bs, KU)
    as.matrix(loadings)
  })
}

#' @keywords internal
#' @noRd
.dkge_standardize_matrix <- function(X, center = TRUE, scale = TRUE) {
  mu <- if (center) colMeans(X) else rep(0, ncol(X))
  sdv <- if (scale) pmax(apply(X, 2, stats::sd), 1e-8) else rep(1, ncol(X))
  Xs <- sweep(sweep(X, 2, mu, FUN = "-"), 2, sdv, FUN = "/")
  list(X = Xs, mean = mu, sd = sdv)
}
