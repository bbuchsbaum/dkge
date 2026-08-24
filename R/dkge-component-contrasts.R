# dkge-component-contrasts.R

#' Construct contrasts that isolate fitted DKGE components
#'
#' Returns the effect-space contrast matrix `C = R %*% U`, where `R` is the
#' pooled-design ruler and `U` is the fitted K-orthonormal basis. These are the
#' contrasts that isolate fitted component coordinates because
#' `crossprod(U, K %*% backsolve(R, C))` equals the identity.
#'
#' This is deliberately different from [dkge_component_saliences()], which
#' returns the dual read-out basis `K %*% U`. Saliences describe how effects load
#' on components; passing `K %*% U` back to [dkge_contrast()] applies the kernel
#' a second time and generally does not isolate components. With LOSO or K-fold
#' contrasts, genuine fold-to-fold basis variation can still mix coordinates
#' relative to the full-data basis.
#'
#' @param fit A fitted `dkge` object.
#' @param comps Components to include. Defaults to all fitted components and
#'   otherwise accepts numeric component indices.
#' @return A numeric effects-by-components matrix suitable for the `contrasts`
#'   argument of [dkge_contrast()].
#' @examples
#' toy <- dkge_sim_toy(
#'   factors = list(cond = list(L = 3)), active_terms = "cond",
#'   S = 3, P = 10, snr = 4
#' )
#' fit <- dkge(toy$B_list, toy$X_list, K = toy$K, rank = 2)
#' C <- dkge_component_contrasts(fit)
#' round(crossprod(fit$U, fit$K %*% backsolve(fit$R, C)), 10)
#' @seealso [dkge_component_saliences()], [dkge_contrast()]
#' @export
dkge_component_contrasts <- function(fit, comps = NULL) {
  stopifnot(inherits(fit, "dkge"))
  if (is.null(comps)) comps <- seq_len(ncol(fit$U))
  comps <- .dkge_component_indexer(fit, comps)
  idx <- .dkge_effect_indexer(fit)

  U_sub <- fit$U[, comps, drop = FALSE]
  contrasts <- fit$R %*% U_sub
  rownames(contrasts) <- idx$names
  colnames(contrasts) <- .dkge_component_labels(comps)
  contrasts
}
