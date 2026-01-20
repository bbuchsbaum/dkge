# dkge-fit-from-kernels.R
# Helper to enter the DKGE pipeline from aligned subject effect kernels.

#' Fit DKGE from precomputed subject effect kernels
#'
#' Converts a list of subject-level effect kernels \eqn{K_s \in \mathbb{R}^{q \times q}} into
#' synthetic GLM inputs that reuse [dkge_fit()] without modifying the core
#' implementation. Each kernel is factorised into a symmetric square root,
#' scaled to keep the pooled design metric unchanged, and paired with an
#' identity design matrix so the resulting DKGE fit matches the supplied
#' kernels.
#'
#' @param K_list List of symmetric positive semi-definite matrices sharing the
#'   same effect ordering.
#' @param effect_ids Character vector of length \eqn{q} naming the shared effect
#'   (anchor) indices.
#' @param subject_ids Optional character vector naming subjects. Defaults to the
#'   names of `K_list` or sequential identifiers.
#' @param design_kernel Optional design kernel passed to [dkge_fit()]. Defaults
#'   to the \eqn{q \times q} identity matrix, which matches the whitened anchor
#'   setup.
#' @param sqrt_tol Eigenvalue tolerance used when extracting square roots.
#' @param ... Additional arguments forwarded to [dkge_fit()].
#'
#' @return A `dkge` object identical to one obtained from [dkge_fit()], with
#'   provenance annotated to record the kernel-driven construction.
#' @export
#'
#' @examples
#' \dontrun{
#' q <- 5
#' Ks <- replicate(3, {
#'   X <- matrix(rnorm(q * q), q)
#'   S <- crossprod(X)
#'   S / sqrt(sum(diag(S)))
#' }, simplify = FALSE)
#' fit <- dkge_fit_from_kernels(Ks, effect_ids = paste0("z", seq_len(q)))
#' }
dkge_fit_from_kernels <- function(K_list,
                                  effect_ids,
                                  subject_ids = NULL,
                                  design_kernel = NULL,
                                  sqrt_tol = 1e-10,
                                  ...) {
  stopifnot(is.list(K_list), length(K_list) >= 1L)
  q <- length(effect_ids)
  stopifnot(is.numeric(q), q >= 1L)

  if (is.null(subject_ids)) {
    subject_ids <- names(K_list)
    if (is.null(subject_ids)) {
      subject_ids <- paste0("subject", seq_len(length(K_list)))
    }
  } else {
    stopifnot(length(subject_ids) == length(K_list))
  }

  effect_ids <- as.character(effect_ids)
  stopifnot(length(effect_ids) == q)

  symmetrize <- function(M) {
    M <- as.matrix(M)
    stopifnot(nrow(M) == q, ncol(M) == q)
    0.5 * (M + t(M))
  }

  sqrt_psd <- function(M) {
    eg <- eigen(M, symmetric = TRUE)
    vals <- pmax(eg$values, 0)
    keep <- vals > sqrt_tol
    if (!any(keep)) {
      matrix(0, q, 0)
    } else {
      eg$vectors[, keep, drop = FALSE] %*%
        (diag(sqrt(vals[keep]), nrow = sum(keep)))
    }
  }

  K_list <- lapply(K_list, symmetrize)
  S <- length(K_list)
  scale_factor <- sqrt(S)

  designs <- vector("list", S)
  betas <- vector("list", S)
  design_template <- diag(1, q)
  rownames(design_template) <- paste0("anchor_trial_", seq_len(q))
  colnames(design_template) <- effect_ids

  for (s in seq_len(S)) {
    sqrtKs <- sqrt_psd(K_list[[s]])
    beta <- sqrtKs / scale_factor
    if (nrow(beta) != q) {
      # zero-rank subject already handled -> q x 0 with proper row count
      beta <- matrix(beta, nrow = q)
    }
    rownames(beta) <- effect_ids
    if (ncol(beta) > 0 && is.null(colnames(beta))) {
      colnames(beta) <- paste0("kernel_factor_", seq_len(ncol(beta)))
    }
    betas[[s]] <- beta
    designs[[s]] <- design_template
  }

  rank_cols <- vapply(betas, ncol, integer(1))
  if (!any(rank_cols > 0)) {
    warning("All subject kernels are numerically rank-0 after tolerance; returning zero factors.",
            call. = FALSE)
  }

  data_bundle <- dkge_data(betas = betas,
                           designs = designs,
                           subject_ids = subject_ids)

  if (is.null(design_kernel)) {
    design_kernel <- diag(1, q)
    dimnames(design_kernel) <- list(effect_ids, effect_ids)
  } else if (is.list(design_kernel) && !is.null(design_kernel$K)) {
    design_kernel <- design_kernel$K
  }

  stopifnot(is.matrix(design_kernel), nrow(design_kernel) == q, ncol(design_kernel) == q)

  fit <- dkge_fit(data_bundle, K = design_kernel, ...)

  provenance <- fit$provenance %||% list()
  provenance$kernel_fit <- list(
    subjects = subject_ids,
    effect_ids = effect_ids,
    sqrt_scale = scale_factor,
    source = "dkge_fit_from_kernels"
  )
  fit$provenance <- provenance
  fit
}
