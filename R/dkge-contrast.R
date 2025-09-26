# dkge-contrast.R
# Unified contrast engine for DKGE with multiple cross-fitting strategies

#' Compute DKGE contrasts with cross-fitting
#'
#' Main entry point for computing design contrasts after DKGE fit, with support
#' for multiple cross-fitting strategies to ensure unbiased estimation.
#'
#' @param fit A `dkge` object from [dkge_fit()] or [dkge()]
#' @param contrasts Either a q-length numeric vector for a single contrast,
#'   a named list of contrasts, or a qxk matrix where columns are contrasts
#' @param method Cross-fitting strategy: "loso" (leave-one-subject-out),
#'   "kfold" (K-fold cross-validation), or "analytic" (first-order approximation)
#' @param folds For method="kfold", either an integer K for random folds,
#'   or a list defining custom fold assignments (see [dkge_define_folds()])
#' @param ridge Optional ridge parameter added when recomputing held-out basis
#' @param parallel Logical; if TRUE uses `future.apply::future_lapply()` for
#'   per-subject work (requires the future.apply package)
#' @param verbose Logical; print progress messages
#' @param align Logical; if TRUE (default) align LOSO/K-fold bases to a common reference and compute a consensus basis for reporting.
#' @param transport Optional list describing how to transport subject-level
#'   contrasts to a shared reference parcellation. Supply either explicit
#'   transport matrices via `transforms`/`matrices`, or configuration for the
#'   medoid/atlas transport helpers (e.g., `method`, `centroids`, `medoid`). When
#'   provided, the resulting transport bundle is stored under
#'   `metadata$transport` for downstream reuse.
#' @param ... Additional arguments passed to method-specific functions
#'
#' @return A list with class `dkge_contrasts` containing:
#'   - `values`: Named list of contrast values (one P_s vector per subject per contrast)
#'   - `method`: Cross-fitting method used
#'   - `contrasts`: Input contrast specifications
#'   - `metadata`: Method-specific metadata (fold assignments, bases, etc.)
#'
#' @details
#' This function provides a unified interface to three cross-fitting strategies:
#'
#' 1. **LOSO** (`method = "loso"`): Recomputes the basis excluding each subject,
#'    then projects that subject's data. This is the gold standard for unbiased
#'    estimation but requires S eigen-decompositions.
#'
#' 2. **K-fold** (`method = "kfold"`): Splits data into K folds, recomputes basis
#'    excluding each fold, projects held-out data. More efficient than LOSO while
#'    maintaining good bias properties. Supports time-based, run-based, or custom
#'    fold definitions.
#'
#' 3. **Analytic** (`method = "analytic"`): Uses first-order eigenvalue perturbation
#'    theory to approximate the LOSO solution without full recomputation. Fast but
#'    may be less accurate when subjects have high leverage.
#'
#' All methods work entirely in the qxq design space and respect the K-metric
#' throughout. Multiple contrasts can be evaluated simultaneously for efficiency.
#'
#' @examples
#' \dontrun{
#' # Single contrast with LOSO
#' c1 <- c(1, -1, 0, 0, 0)  # Main effect difference
#' result <- dkge_contrast(fit, c1, method = "loso")
#'
#' # Multiple contrasts with K-fold
#' contrasts <- list(
#'   main1 = c(1, -1, 0, 0, 0),
#'   main2 = c(0, 0, 1, -1, 0),
#'   interaction = c(1, -1, -1, 1, 0)
#' )
#' result <- dkge_contrast(fit, contrasts, method = "kfold", folds = 5)
#'
#' # Fast analytic approximation
#' result <- dkge_contrast(fit, c1, method = "analytic")
#' }
#'
#' @seealso [dkge_loso_contrast()], [dkge_define_folds()], [dkge_infer()]
#' @export
dkge_contrast <- function(fit, contrasts,
                         method = c("loso", "kfold", "analytic"),
                         folds = NULL,
                         ridge = 0,
                         parallel = FALSE,
                         verbose = FALSE,
                         align = TRUE,
                         transport = NULL,
                         ...) {
  stopifnot(inherits(fit, "dkge"))
  method <- match.arg(method)

  # Normalize contrast input
  contrast_list <- .normalize_contrasts(contrasts, fit)

  # Dispatch to method
  result <- switch(method,
    loso = .dkge_contrast_loso(fit, contrast_list, ridge, parallel, verbose, align = align, ...),
    kfold = .dkge_contrast_kfold(fit, contrast_list, folds, ridge, parallel, verbose, align = align, ...),
    analytic = .dkge_contrast_analytic(fit, contrast_list, ridge, parallel, verbose, align = align, ...)
  )

  if (is.null(result$metadata)) {
    result$metadata <- list()
  }
  if (length(result$values) > 0) {
    first_values <- result$values[[1]]
    subject_ids <- names(first_values)
    if (is.null(subject_ids)) {
      subject_ids <- as.character(seq_along(first_values))
    }
    cluster_dims <- vapply(first_values, length, integer(1))
    names(cluster_dims) <- subject_ids
    result$metadata$cluster_dims <- cluster_dims
  } else {
    result$metadata$cluster_dims <- integer(0)
  }

  if (!is.null(transport)) {
    warning("`transport` argument to dkge_contrast() is deprecated; use `dkge_transport_contrasts_to_medoid()`.",
            call. = FALSE)
  }

  structure(result, class = "dkge_contrasts")
}

#' Normalize contrast specifications
#'
#' @param contrasts Various input formats
#' @param fit dkge object for validation
#' @return Named list of q-length numeric vectors
#' @keywords internal
#' @noRd
.normalize_contrasts <- function(contrasts, fit) {
  q <- nrow(fit$U)

  if (is.numeric(contrasts) && is.null(dim(contrasts))) {
    # Single vector
    stopifnot(length(contrasts) == q)
    return(list(contrast1 = as.numeric(contrasts)))
  }

  if (is.matrix(contrasts)) {
    # Matrix: columns are contrasts
    stopifnot(nrow(contrasts) == q)
    cn <- colnames(contrasts)
    if (is.null(cn)) cn <- paste0("contrast", seq_len(ncol(contrasts)))
    contrast_list <- lapply(seq_len(ncol(contrasts)), function(j) contrasts[, j])
    names(contrast_list) <- cn
    return(contrast_list)
  }

  if (is.list(contrasts)) {
    # Named list
    stopifnot(all(vapply(contrasts, length, integer(1)) == q))
    return(lapply(contrasts, as.numeric))
  }

  stop("contrasts must be a numeric vector, matrix, or named list")
}

#' LOSO contrast implementation
#'
#' @inheritParams dkge_contrast
#' @param contrast_list Normalized list of contrasts
#' @keywords internal
#' @noRd
.dkge_contrast_loso <- function(fit, contrast_list, ridge, parallel, verbose, align, ...) {
  S <- length(fit$Btil)
  n_contrasts <- length(contrast_list)
  verbose_flag <- .dkge_verbose(verbose)

  if (verbose_flag) {
    message(sprintf("Computing %d contrast(s) via LOSO for %d subjects", n_contrasts, S))
  }

  subject_labels <- fit$subject_ids %||% paste0("subject", seq_len(S))
  assignments <- lapply(seq_len(S), function(s) s)

  fold_info <- .dkge_build_fold_bases(
    fit,
    assignments = assignments,
    ridge = ridge,
    align = align,
    loader_scope = "heldout",
    verbose = verbose
  )
  folds <- fold_info$folds

  c_tilde_list <- lapply(contrast_list, function(ct) backsolve(fit$R, ct, transpose = FALSE))

  values <- vector("list", n_contrasts)
  names(values) <- names(contrast_list)
  alphas <- vector("list", n_contrasts)

  r <- ncol(fit$U)
  fold_row_names <- vapply(folds, function(fold) paste(subject_labels[fold$subjects], collapse = ","), character(1))

  for (i in seq_along(contrast_list)) {
    values[[i]] <- vector("list", S)
    names(values[[i]]) <- subject_labels
    alpha_mat <- matrix(NA_real_, nrow = length(folds), ncol = r)
    rownames(alpha_mat) <- fold_row_names

    for (fold in folds) {
      U_fold <- fold$basis
      alpha_vec <- as.numeric(t(U_fold) %*% fit$K %*% c_tilde_list[[i]])
      alpha_mat[fold$index, seq_along(alpha_vec)] <- alpha_vec

      holdout <- fold$subjects
      loaders <- fold$loaders
      subject_scores <- .dkge_apply(
        holdout,
        function(s) {
          loader <- loaders[[as.character(s)]]
          as.numeric(loader$A %*% alpha_vec)
        },
        parallel = parallel
      )

      for (pos in seq_along(holdout)) {
        s <- holdout[[pos]]
        values[[i]][[subject_labels[[s]]]] <- subject_scores[[pos]]
      }
    }

    alphas[[i]] <- alpha_mat
  }

  for (i in seq_along(values)) {
    missing <- which(vapply(values[[i]], is.null, logical(1)))
    if (length(missing) > 0) {
      warning(sprintf("Subjects %s not processed for contrast %s",
                      paste(subject_labels[missing], collapse = ","),
                      names(values)[i] %||% paste0("contrast", i)))
    }
  }

  metadata <- list(
    bases = lapply(folds, `[[`, "basis"),
    aligned_bases = lapply(folds, `[[`, "basis_aligned"),
    rotations = lapply(folds, `[[`, "rotation"),
    alphas = alphas,
    ridge = ridge,
    procrustes = if (align) list(alignment = fold_info$alignment, consensus = fold_info$consensus) else NULL
  )

  list(
    values = values,
    method = "loso",
    contrasts = contrast_list,
    metadata = metadata
  )
}

#' Analytic LOSO contrast implementation
#'
#' @inheritParams .dkge_contrast_loso
#' @keywords internal
#' @noRd
.dkge_contrast_analytic <- function(fit, contrast_list, ridge, parallel, verbose, align = TRUE, ...) {
  .dkge_contrast_analytic_impl(fit, contrast_list, ridge, parallel, verbose, align = align, ...)
}

#' Print method for dkge_contrasts
#'
#' @param x A dkge_contrasts object
#' @param ... Additional arguments (unused)
#' @export
print.dkge_contrasts <- function(x, ...) {
  n_contrasts <- length(x$contrasts)
  n_subjects <- length(x$values[[1]])

  cat("DKGE Contrasts\n")
  cat("--------------\n")
  cat(sprintf("Method: %s\n", x$method))
  cat(sprintf("Contrasts: %d\n", n_contrasts))
  cat(sprintf("Subjects: %d\n", n_subjects))

  if (n_contrasts <= 5) {
    cat("Contrast names:", paste(names(x$contrasts), collapse = ", "), "\n")
  } else {
    cat("Contrast names:", paste(names(x$contrasts)[1:5], collapse = ", "), "...\n")
  }

  if (!is.null(x$metadata$ridge) && x$metadata$ridge > 0) {
    cat(sprintf("Ridge: %g\n", x$metadata$ridge))
  }

  detail <- x$metadata$fallback_detail
  if (!is.null(detail) && nrow(detail)) {
    fallback_rows <- detail[is.na(detail$reason) | detail$reason != "analytic", , drop = FALSE]
    fallback_count <- nrow(fallback_rows)
    if (fallback_count > 0) {
      cat(sprintf("Fallback triggered for %d/%d subject-contrast pairs (see metadata$fallback_detail).\n",
                  fallback_count, nrow(detail)))
      top_rows <- head(fallback_rows, 5)
      for (k in seq_len(nrow(top_rows))) {
        row <- top_rows[k, , drop = FALSE]
        cat(sprintf("    subject=%s, contrast=%s, reason=%s\n",
                    row$subject, row$contrast, row$reason))
      }
      if (fallback_count > nrow(top_rows)) {
        cat("    ...\n")
      }
    }
  }

  invisible(x)
}

#' Extract contrast values as matrix
#'
#' @param x A dkge_contrasts object
#' @param contrast Name or index of contrast to extract
#' @param ... Additional arguments (not used)
#' @return SxP matrix of contrast values
#' @export
as.matrix.dkge_contrasts <- function(x, contrast = 1, ...) {
  if (is.character(contrast)) {
    contrast <- match(contrast, names(x$contrasts))
    if (is.na(contrast)) stop("Contrast not found")
  }

  value_list <- x$values[[contrast]]
  dims <- vapply(value_list, length, integer(1))
  unique_dims <- unique(dims)

  if (length(unique_dims) != 1) {
    msg <- "Subject cluster counts differ; use dkge_transport_contrasts_to_medoid() before stacking."
    cond <- structure(list(message = msg, call = sys.call()),
                     class = c("dkge_transport_needed", "error", "condition"))
    stop(cond)
  }

  do.call(rbind, value_list)
}

#' @export
as.data.frame.dkge_contrasts <- function(x, row.names = NULL, optional = FALSE, ...,
                                         stringsAsFactors = FALSE) {
  contrast_names <- names(x$contrasts)
  if (is.null(contrast_names) || any(!nzchar(contrast_names))) {
    contrast_names <- paste0("contrast", seq_along(x$contrasts))
  }

  rows <- vector("list", length(x$values))
  for (i in seq_along(x$values)) {
    subj_values <- x$values[[i]]
    if (!length(subj_values)) {
      next
    }
    subj_names <- names(subj_values)
    if (is.null(subj_names) || any(!nzchar(subj_names))) {
      subj_names <- paste0("subject", seq_along(subj_values))
    }

    entries <- Map(function(vals, subj) {
      cluster_ids <- names(vals)
      if (is.null(cluster_ids) || any(!nzchar(cluster_ids))) {
        cluster_ids <- paste0("cluster", seq_along(vals))
      }
      data.frame(
        contrast = contrast_names[[i]],
        subject = subj,
        component = cluster_ids,
        value = as.numeric(vals),
        method = x$method,
        stringsAsFactors = stringsAsFactors
      )
    }, subj_values, subj_names)
    rows[[i]] <- do.call(rbind, entries)
  }

  rows <- Filter(Negate(is.null), rows)
  result <- if (length(rows)) do.call(rbind, rows) else NULL
  if (is.null(result)) {
    result <- data.frame(
      contrast = character(0),
      subject = character(0),
      component = character(0),
      value = numeric(0),
      method = character(0),
      stringsAsFactors = stringsAsFactors
    )
  }
  if (!is.null(row.names)) {
    rownames(result) <- row.names
  } else {
    rownames(result) <- NULL
  }
  result
}
