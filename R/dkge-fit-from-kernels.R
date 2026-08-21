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
#'   provenance annotated to record the kernel-driven construction. Effect
#'   matrices use `effect_ids` as canonical dimnames, subject-indexed lists use
#'   `subject_ids`, and component axes are labelled `LV1`, `LV2`, and so on.
#' @export
#'
#' @examples
#' \donttest{
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
  effect_ids <- as.character(effect_ids)
  q <- length(effect_ids)
  if (q < 1L || anyNA(effect_ids) || any(!nzchar(effect_ids)) ||
      anyDuplicated(effect_ids)) {
    stop("`effect_ids` must contain unique, non-empty effect labels.",
         call. = FALSE)
  }

  if (is.null(subject_ids)) {
    subject_ids <- names(K_list)
    if (is.null(subject_ids)) {
      subject_ids <- paste0("subject", seq_len(length(K_list)))
    }
  } else {
    if (length(subject_ids) != length(K_list)) {
      stop("`subject_ids` must have one entry per subject kernel.",
           call. = FALSE)
    }
  }
  subject_ids <- as.character(subject_ids)
  if (anyNA(subject_ids) || any(!nzchar(subject_ids)) ||
      anyDuplicated(subject_ids)) {
    stop("`subject_ids` must contain unique, non-empty subject labels.",
         call. = FALSE)
  }

  symmetrize <- function(M) {
    M <- as.matrix(M)
    if (!is.numeric(M) || !identical(dim(M), c(q, q)) ||
        any(!is.finite(M))) {
      stop("Each subject kernel must be a finite numeric q x q matrix.",
           call. = FALSE)
    }

    rn <- rownames(M)
    cn <- colnames(M)
    if (is.null(rn) && is.null(cn)) {
      dimnames(M) <- list(effect_ids, effect_ids)
    } else {
      # A single named axis is sufficient for a symmetric kernel: the unnamed
      # axis necessarily shares its positional ordering. When both axes are
      # named, reorder them independently before any numeric operation so a
      # differently ordered but correctly labelled matrix retains its
      # name-to-value mapping.
      if (is.null(rn)) rn <- cn
      if (is.null(cn)) cn <- rn
      valid_axis <- function(x) {
        length(x) == q && !anyNA(x) && all(nzchar(x)) &&
          !anyDuplicated(x) && setequal(x, effect_ids)
      }
      if (!valid_axis(rn) || !valid_axis(cn)) {
        stop(
          "Subject-kernel dimnames must be unique permutations of `effect_ids`.",
          call. = FALSE
        )
      }
      M <- M[match(effect_ids, rn), match(effect_ids, cn), drop = FALSE]
      dimnames(M) <- list(effect_ids, effect_ids)
    }

    out <- 0.5 * (M + t(M))
    dimnames(out) <- list(effect_ids, effect_ids)
    out
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
  names(K_list) <- subject_ids
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

  component_ids <- paste0("LV", seq_len(ncol(fit$U)))
  full_component_ids <- paste0("LV", seq_len(ncol(fit$eig_vectors_full)))
  label_effect_components <- function(x, components = component_ids) {
    if (is.matrix(x) && nrow(x) == q && ncol(x) == length(components)) {
      dimnames(x) <- list(effect_ids, components)
    }
    x
  }
  label_effect_square <- function(x) {
    if (is.matrix(x) && identical(dim(x), c(q, q))) {
      dimnames(x) <- list(effect_ids, effect_ids)
    }
    x
  }
  label_subject_matrices <- function(x, square = FALSE) {
    if (!is.list(x) || length(x) != length(subject_ids)) return(x)
    x <- lapply(x, function(M) {
      if (!is.matrix(M)) return(M)
      if (square) {
        label_effect_square(M)
      } else if (nrow(M) == q) {
        rownames(M) <- effect_ids
        M
      } else {
        M
      }
    })
    names(x) <- subject_ids
    x
  }

  for (nm in c("R", "K", "Khalf", "Kihalf", "Chat", "Chat_sym",
               "effect_moment", "pair_counts", "pair_weight", "pair_ess")) {
    fit[[nm]] <- label_effect_square(fit[[nm]])
  }
  fit$U <- label_effect_components(fit$U)
  fit$s <- label_effect_components(fit$s)
  fit$KU <- label_effect_components(fit$KU)
  fit$scores_matrix <- label_effect_components(fit$scores_matrix)
  fit$eig_vectors_full <- label_effect_components(
    fit$eig_vectors_full, full_component_ids
  )
  if (is.matrix(fit$v) && ncol(fit$v) == length(component_ids)) {
    colnames(fit$v) <- component_ids
  }
  if (length(fit$sdev) == length(component_ids)) names(fit$sdev) <- component_ids
  if (length(fit$evals) == length(full_component_ids)) {
    names(fit$evals) <- full_component_ids
  }
  if (length(fit$eig_values_full) == length(full_component_ids)) {
    names(fit$eig_values_full) <- full_component_ids
  }
  for (nm in c("contribs", "effect_moments", "effect_moments_raw",
               "noise_moments")) {
    fit[[nm]] <- label_subject_matrices(fit[[nm]], square = TRUE)
  }
  for (nm in c("Braw", "Btil", "Omega", "subjects", "effect_precision")) {
    fit[[nm]] <- label_subject_matrices(fit[[nm]], square = FALSE)
  }
  if (is.list(fit$block_indices) &&
      length(fit$block_indices) == length(subject_ids)) {
    names(fit$block_indices) <- subject_ids
  }

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
