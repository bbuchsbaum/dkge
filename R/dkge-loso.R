# dkge-loso.R
# Leave-one-subject-out DKGE contrasts.


#' Leave-one-subject-out DKGE contrast
#'
#' @param fit `dkge` object
#' @param s Subject index (1-based)
#' @param contrasts Contrast vector in the original design basis
#' @param ridge Optional ridge when recomputing the held-out compressed matrix
#' @return List with fields `v`, `alpha`, and `basis`
#' @keywords internal
#' @export
#' @examples
#' \donttest{
#' toy <- dkge_sim_toy(
#'   factors = list(A = list(L = 2), B = list(L = 3)),
#'   active_terms = c("A", "B"), S = 4, P = 20, snr = 5
#' )
#' fit <- dkge_fit(toy$B_list, toy$X_list, toy$K, rank = 2)
#' c_vec <- c(1, -1, rep(0, 3))
#' result <- dkge_loso_contrast(fit, s = 1, contrasts = c_vec)
#' }
dkge_loso_contrast <- function(fit, s, contrasts, ridge = 0) {
  stopifnot(inherits(fit, "dkge"), s >= 1L, s <= length(fit$Btil))
  q <- nrow(fit$U)
  stopifnot(length(contrasts) == q)

  train_ids <- setdiff(seq_len(length(fit$Btil)), s)
  ctx <- .dkge_fold_weight_context(fit, train_ids, ridge = ridge)
  Chat_minus <- ctx$Chat
  weight_eval <- ctx$weights

  eig_minus <- eigen(Chat_minus, symmetric = TRUE)
  r <- ncol(fit$U)
  Uminus <- fit$Kihalf %*% eig_minus$vectors[, seq_len(r), drop = FALSE]

  c_tilde <- backsolve(fit$R, contrasts, transpose = FALSE)
  alpha <- t(Uminus) %*% fit$K %*% c_tilde

  loader_weights <- weight_eval$total
  Bts <- fit$Btil[[s]]
  Bw <- if (is.null(loader_weights)) Bts else sweep(Bts, 2L, sqrt(pmax(loader_weights, 0)), "*")
  A_s <- t(Bw) %*% fit$K %*% Uminus
  v_s <- as.numeric(A_s %*% alpha)

  list(v = v_s, alpha = alpha, basis = Uminus, evals = eig_minus$values)
}
