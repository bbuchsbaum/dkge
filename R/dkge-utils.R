# dkge-utils.R
# Shared helper utilities for DKGE

#' Null-coalescing helper
#'
#' Returns `b` when `a` is `NULL`, otherwise returns `a`.
#'
#' @name grapes-or-or-grapes
#' @keywords internal
NULL

#' @rdname grapes-or-or-grapes
#' @param a Primary value tested for `NULL`.
#' @param b Fallback value returned when `a` is `NULL`.
#' @usage a \%||\% b
#' @keywords internal
`%||%` <- function(a, b) if (is.null(a)) b else a

#' Apply helper with optional parallelism
#'
#' Wraps `lapply()` with an optional future.apply backend so callers can enable
#' `parallel = TRUE` without repeating boilerplate dependency checks.
#'
#' @param X Vector or list to iterate over.
#' @param FUN Function to apply.
#' @param parallel Logical; if `TRUE`, uses `future.apply::future_lapply()`.
#' @param ... Additional arguments passed to the apply backend.
#' @return List of results matching `lapply()` semantics.
#' @keywords internal
.dkge_apply <- function(X, FUN, parallel = FALSE, ...) {
  if (parallel) {
    if (!requireNamespace("future.apply", quietly = TRUE)) {
      stop("parallel=TRUE requires the future.apply package; install it or set parallel=FALSE.",
           call. = FALSE)
    }
    future.apply::future_lapply(X, FUN, ...)
  } else {
    lapply(X, FUN, ...)
  }
}

#' Check whether verbose output should be emitted
#'
#' Uses the per-call `verbose` flag combined with the global
#' `options(dkge.verbose = TRUE)` toggle.
#'
#' @keywords internal
.dkge_verbose <- function(verbose) {

  isTRUE(verbose) && isTRUE(getOption("dkge.verbose", TRUE))
}

# -------------------------------------------------------------------------
# Numerical robustness utilities ------------------------------------------
# -------------------------------------------------------------------------

#' Check matrix rank for design and/or beta matrices
#'
#' Detects rank deficiency and emits informative warnings identifying the
#' culprit subject and the nature of the problem.
#'
#' @param design Design matrix to check (T x q).
#' @param beta Optional beta matrix to check (q x P).
#' @param subject_id Optional subject identifier for warning messages.
#' @return List with `design_rank` and `beta_rank` (if beta provided).
#' @keywords internal
#' @noRd
.dkge_check_rank <- function(design, beta = NULL, subject_id = NULL) {
  subject_label <- subject_id %||% "(unnamed)"
  result <- list(design_rank = NULL, beta_rank = NULL)


  if (!is.null(design) && is.matrix(design)) {
    qr_design <- qr(design)
    design_rank <- qr_design$rank
    expected_rank <- ncol(design)
    result$design_rank <- design_rank

    if (design_rank < expected_rank) {
      warning(sprintf(
        "Subject '%s': design matrix is rank-deficient (rank %d < %d columns). Effects may be aliased.",
        subject_label, design_rank, expected_rank
      ), call. = FALSE)
    }
  }

  if (!is.null(beta) && is.matrix(beta)) {
    beta_rank <- qr(beta)$rank
    result$beta_rank <- beta_rank

    if (beta_rank < nrow(beta)) {
      warning(sprintf(
        "Subject '%s': beta matrix has reduced rank (%d < %d effects).",
        subject_label, beta_rank, nrow(beta)
      ), call. = FALSE)
    }
  }

  result
}

#' Check matrix condition number against threshold
#'
#' Warns if the condition number exceeds the specified threshold, indicating
#' potential numerical instability.
#'
#' @param M Symmetric matrix to check.
#' @param threshold Condition number threshold (default 1e8).
#' @param name Descriptive name for the matrix (used in warning message).
#' @return The computed condition number.
#' @keywords internal
#' @noRd
.dkge_check_condition <- function(M, threshold = 1e8, name = "matrix") {
  cond <- kappa(M, exact = FALSE)
  if (cond > threshold) {
    warning(sprintf(
      "%s is ill-conditioned (condition number: %.2e > %.2e threshold). Results may be numerically unstable.",
      name, cond, threshold
    ), call. = FALSE)
  }
  cond
}

#' Identify and track voxels with non-finite values
#'
#' Scans each beta matrix for columns containing NA, NaN, or Inf values,
#' emits per-subject warnings, and returns metadata about exclusions.
#'
#' @param B_list List of beta matrices (q x P_s each).
#' @param subject_ids Optional character vector of subject identifiers.
#' @return List with `excluded_voxels` (list of integer vectors per subject),
#'   `excluded_counts` (integer vector), and `total_excluded` (integer).
#' @keywords internal
#' @noRd
.dkge_voxel_exclusion_mask <- function(B_list, subject_ids = NULL) {
  S <- length(B_list)
  excluded_voxels <- vector("list", S)
  excluded_counts <- integer(S)

  for (s in seq_len(S)) {
    B <- B_list[[s]]
    if (is.null(B) || !is.matrix(B) || ncol(B) == 0) {
      excluded_voxels[[s]] <- integer(0)
      excluded_counts[s] <- 0L
      next
    }

    bad_cols <- which(colSums(!is.finite(B)) > 0)
    excluded_voxels[[s]] <- bad_cols
    excluded_counts[s] <- length(bad_cols)

    if (length(bad_cols) > 0) {
      pct <- 100 * length(bad_cols) / ncol(B)
      subject_label <- if (!is.null(subject_ids) && length(subject_ids) >= s) {
        subject_ids[s]
      } else {
        as.character(s)
      }
      warning(sprintf(
        "Subject '%s': %d voxels (%.1f%%) excluded due to NA/NaN/Inf values.",
        subject_label, length(bad_cols), pct
      ), call. = FALSE)
    }
  }

  list(
    excluded_voxels = excluded_voxels,
    excluded_counts = excluded_counts,
    total_excluded = sum(excluded_counts)
  )
}
