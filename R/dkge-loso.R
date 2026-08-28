# dkge-loso.R
# Leave-one-subject-out DKGE contrasts.


#' Leave-one-subject-out DKGE contrast
#'
#' @param fit `dkge` object
#' @param s Subject index (1-based)
#' @param contrasts Contrast vector in the original design basis
#' @param ridge Optional ridge when recomputing the held-out compressed matrix
#' @return List with fields `v`, `alpha`, `basis`, `loadings`, and a typed
#'   `alignment_receipt` binding those values to their held-out training and
#'   preprocessing provenance.
#' @details Exact estimator replay currently supports pooled, non-CPCA fits.
#'   JD and CPCA fits fail closed rather than substituting an ordinary pooled
#'   eigensolve and mislabelling it as the fitted estimator.
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
  .dkge_assert_crossfit_estimator_supported(fit, "LOSO contrast estimation")
  q <- nrow(fit$U)
  stopifnot(length(contrasts) == q)
  .dkge_validate_kernel_contrasts(list(contrast1 = as.numeric(contrasts)), fit)

  train_ids <- setdiff(seq_len(length(fit$Btil)), s)
  ctx <- .dkge_fold_weight_context(fit, train_ids, ridge = ridge)
  Chat_minus <- ctx$Chat
  weight_eval <- ctx$weights

  eig_minus <- eigen(Chat_minus, symmetric = TRUE)
  r <- ncol(fit$U)
  fold_contract <- .dkge_spectral_contract(eig_minus$values)
  fold_rank <- min(fit$kernel_rank %||% qr(fit$K)$rank,
                   fold_contract$rank)
  if (fold_rank < r) {
    .dkge_abort(
      sprintf(
        "LOSO training data have effective rank %d, below fitted rank %d; refit at rank <= %d.",
        fold_rank, r, fold_rank
      ),
      "dkge_fold_rank_error"
    )
  }
  Uminus <- fit$Kihalf %*% eig_minus$vectors[, seq_len(r), drop = FALSE]
  Uminus <- dkge_k_orthonormalize(Uminus, fit$K)

  c_tilde <- backsolve(fit$R, contrasts, transpose = FALSE)
  alpha <- t(Uminus) %*% fit$K %*% c_tilde

  Bts <- fit$Btil[[s]]
  loader_weights <- .dkge_subject_loader_weights(weight_eval$total, Bts)
  Bw <- if (is.null(loader_weights)) Bts else sweep(Bts, 2L, sqrt(pmax(loader_weights, 0)), "*")
  Bmodel <- .dkge_apply_fit_spatial(fit, Bw, subject = s)
  A_s <- t(Bmodel) %*% fit$K %*% Uminus
  v_s <- as.numeric(A_s %*% alpha)

  preprocessing <- .dkge_alignment_preprocessing_receipt(
    fit, s, Bts, voxel_weights = loader_weights
  )
  alignment_receipt <- .dkge_make_alignment_receipt(
    fit = fit,
    subject = s,
    train_ids = train_ids,
    basis = Uminus,
    evals = eig_minus$values,
    loadings = A_s,
    preprocessing = preprocessing,
    alphas = list(contrast1 = as.numeric(alpha)),
    fold_index = s,
    holdout = s,
    subject_weights = ctx$subject_weights,
    subject_weight_source = ctx$subject_weight_source,
    method = "loso",
    eligible = TRUE,
    eligibility_reason = "exact_heldout_basis"
  )

  list(
    v = v_s,
    alpha = alpha,
    basis = Uminus,
    evals = eig_minus$values,
    loadings = A_s,
    training_subject_indices = train_ids,
    training_subject_ids = alignment_receipt$training_subject_ids,
    basis_id = alignment_receipt$basis_id,
    basis_hash = alignment_receipt$basis_hash,
    preprocessing = preprocessing,
    alignment_receipt = alignment_receipt
  )
}
