# dkge-kfold.R
# K-fold cross-fitting for DKGE contrasts

#' Define folds for K-fold cross-validation
#'
#' Create fold assignments for subjects or time points in DKGE analysis.
#' Supports multiple strategies including random subject-level folds, time-based
#' splits, run-based partitions, and custom user-defined assignments.
#'
#' @param fit A `dkge` object or `dkge_data` bundle
#' @param type Type of fold definition:
#'   - `"subject"`: Random assignment of subjects to folds (default)
#'   - `"time"`: Split each subject's time series into temporal blocks
#'   - `"run"`: Use experimental runs as natural folds
#'   - `"custom"`: User provides fold assignments directly
#' @param k Number of folds (ignored for type="custom")
#' @param runs For type="run", a list of run indicators per subject
#' @param assignments For type="custom", a list of fold assignments
#' @param seed Random seed for reproducible fold assignment
#' @param align Logical; if TRUE (default) compute Procrustes alignment/consensus when folds are evaluated.
#' @param ... Additional arguments for specific fold types
#'
#' @return A `dkge_folds` object containing:
#'   - `type`: The fold type used
#'   - `k`: Number of folds
#'   - `assignments`: List specifying which data belongs to each fold
#'   - `metadata`: Additional information about fold creation
#'
#' @details
#' Different fold types serve different purposes:
#'
#' **Subject-level folds** (`type = "subject"`): Assigns entire subjects to folds.
#' This maintains subject independence and is appropriate when subjects are
#' exchangeable. Each fold will have approximately S/k subjects.
#'
#' **Time-based folds** (`type = "time"`): Splits each subject's time series into
#' k temporal blocks. Useful for assessing temporal stability or when early vs
#' late responses differ. Requires access to original time series dimensions.
#'
#' **Run-based folds** (`type = "run"`): Uses experimental runs as natural folds.
#' Common in fMRI where runs provide natural breaks. Requires run indicators.
#'
#' **Custom folds** (`type = "custom"`): Full control over fold assignments.
#' Supply a list where each element specifies indices for that fold.
#'
#' @examples
#' \dontrun{
#' # Random subject-level 5-fold
#' folds <- dkge_define_folds(fit, type = "subject", k = 5)
#'
#' # Time-based 3-fold (early, middle, late)
#' folds <- dkge_define_folds(fit, type = "time", k = 3)
#'
#' # Run-based using experimental structure
#' run_list <- list(sub1 = c(1,1,2,2,3,3), sub2 = c(1,1,2,2,3,3))
#' folds <- dkge_define_folds(fit, type = "run", runs = run_list)
#'
#' # Custom fold specification
#' folds <- dkge_define_folds(fit, type = "custom",
#'   assignments = list(fold1 = 1:10, fold2 = 11:20, fold3 = 21:30))
#' }
#'
#' @export
dkge_define_folds <- function(fit, type = c("subject", "time", "run", "custom"),
                             k = 5, runs = NULL, assignments = NULL,
                             seed = NULL, align = FALSE, ...) {
  type <- match.arg(type)

  # Extract data info
  if (inherits(fit, "dkge")) {
    n_subjects <- length(fit$Btil)
    subject_ids <- fit$subject_ids
  } else if (inherits(fit, "dkge_data")) {
    n_subjects <- fit$n_subjects
    subject_ids <- fit$subject_ids
  } else {
    stop("fit must be a dkge or dkge_data object")
  }

  if (!is.null(seed)) set.seed(seed)

  folds <- switch(type,
    subject = .define_subject_folds(n_subjects, k, subject_ids, seed = seed),
    time = .define_time_folds(fit, k, ...),
    run = .define_run_folds(fit, runs, ...),
    custom = .validate_custom_folds(assignments, n_subjects)
  )
  if (is.null(folds$metadata)) folds$metadata <- list()
  folds$metadata$seed <- seed
  folds$align <- align

  structure(folds, class = "dkge_folds")
}

#' Random subject-level fold assignment
#' @keywords internal
#' @noRd
.define_subject_folds <- function(n_subjects, k, subject_ids = NULL, seed = NULL) {
  stopifnot(k >= 2, k <= n_subjects)

  # Random permutation then split
  perm <- sample(n_subjects)
  fold_sizes <- rep(n_subjects %/% k, k)
  remainder <- n_subjects %% k
  if (remainder > 0) {
    fold_sizes[seq_len(remainder)] <- fold_sizes[seq_len(remainder)] + 1
  }

  assignments <- vector("list", k)
  start <- 1
  for (i in seq_len(k)) {
    end <- start + fold_sizes[i] - 1
    assignments[[i]] <- sort(perm[start:end])
    start <- end + 1
  }

  list(
    type = "subject",
    k = k,
    assignments = assignments,
    subject_ids = subject_ids,
    metadata = list(
      n_subjects = n_subjects,
      fold_sizes = fold_sizes,
      seed = seed,
      permutation = perm
    )
  )
}

#' Time-based fold definition (placeholder)
#' @keywords internal
#' @noRd
.define_time_folds <- function(fit, k, ...) {
  # This requires access to time series length per subject
  # Will be fully implemented when we have streaming/time-series integration
  stop("Time-based folds require access to original time series. Coming soon.")
}

#' Run-based fold definition
#' @keywords internal
#' @noRd
.define_run_folds <- function(fit, runs, ...) {
  stop("Run-level folds require per-run DKGE contributions; not yet supported.")
}

#' Validate custom fold assignments
#' @keywords internal
#' @noRd
.validate_custom_folds <- function(assignments, n_subjects) {
  if (is.null(assignments)) {
    stop("assignments required for type='custom'")
  }

  stopifnot(is.list(assignments), length(assignments) >= 2)

  # Check all indices are valid
  all_idx <- unlist(assignments)
  stopifnot(all(all_idx >= 1), all(all_idx <= n_subjects))

  # Warn about overlap
  if (length(all_idx) != length(unique(all_idx))) {
    warning("Some subjects appear in multiple folds")
  }

  list(
    type = "custom",
    k = length(assignments),
    assignments = assignments,
    metadata = list(
      n_subjects = n_subjects,
      coverage = length(unique(all_idx))
    )
  )
}

#' K-fold cross-fitted DKGE contrasts
#'
#' Internal implementation called by dkge_contrast() for method="kfold".
#' Rebuilds the basis excluding each fold and projects held-out data.
#'
#' @param fit dkge object
#' @param contrast_list List of normalized contrasts
#' @param folds Either integer k or dkge_folds object
#' @param ridge Ridge parameter for held-out basis
#' @param parallel Logical; enables future.apply-based parallelism for
#'   per-subject projections (requires future.apply)
#' @param verbose Print progress
#' @param ... Additional arguments
#' @return List with values, metadata, etc.
#' @keywords internal
#' @noRd
.dkge_contrast_kfold <- function(fit, contrast_list, folds, ridge,
                                parallel, verbose, align = FALSE, ...) {
  # Prepare folds
  if (is.numeric(folds) && length(folds) == 1) {
    folds <- dkge_define_folds(fit, type = "subject", k = folds)
  } else if (!inherits(folds, "dkge_folds")) {
    folds <- as_dkge_folds(folds, fit_or_data = fit)
    if (!inherits(folds, "dkge_folds")) {
      stop("folds must be an integer k or convertible via as_dkge_folds().")
    }
  }

  S <- length(fit$Btil)
  q <- nrow(fit$U)
  r <- ncol(fit$U)
  n_contrasts <- length(contrast_list)
  k <- folds$k
  verbose_flag <- .dkge_verbose(verbose)

  if (verbose_flag) {
    message(sprintf("Computing %d contrast(s) via %d-fold CV", n_contrasts, k))
  }

  subject_labels <- fit$subject_ids %||% paste0("subject", seq_len(S))

  fold_info <- .dkge_build_fold_bases(
    fit,
    assignments = folds$assignments,
    ridge = ridge,
    align = align,
    loader_scope = "heldout",
    verbose = verbose
  )

  folds_internal <- fold_info$folds
  c_tilde_list <- lapply(contrast_list, function(ct) backsolve(fit$R, ct, transpose = FALSE))

  values <- vector("list", n_contrasts)
  names(values) <- names(contrast_list)
  fold_alphas <- vector("list", n_contrasts)

  fold_row_names <- vapply(folds_internal, function(fold) paste(subject_labels[fold$subjects], collapse = ","), character(1))

  for (i in seq_along(contrast_list)) {
    values[[i]] <- vector("list", S)
    names(values[[i]]) <- subject_labels
    alpha_mat <- matrix(NA_real_, nrow = length(folds_internal), ncol = r)
    rownames(alpha_mat) <- fold_row_names

    for (fold in folds_internal) {
      U_fold <- fold$basis
      alpha_vec <- as.numeric(t(U_fold) %*% fit$K %*% c_tilde_list[[i]])
      alpha_mat[fold$index, seq_along(alpha_vec)] <- alpha_vec

      holdout <- fold$subjects
      loaders <- fold$loaders
      subject_values <- .dkge_apply(
        holdout,
        function(s) {
          loader <- loaders[[as.character(s)]]
          as.numeric(loader$A %*% alpha_vec)
        },
        parallel = parallel
      )

      for (idx in seq_along(holdout)) {
        s <- holdout[[idx]]
        values[[i]][[subject_labels[[s]]]] <- subject_values[[idx]]
      }
    }

    fold_alphas[[i]] <- alpha_mat
  }

  for (i in seq_along(values)) {
    missing <- which(vapply(values[[i]], is.null, logical(1)))
    if (length(missing) > 0) {
      warning(sprintf("Subjects %s not in any fold for contrast %s",
                      paste(subject_labels[missing], collapse = ","),
                      names(values)[i] %||% paste0("contrast", i)))
    }
  }

  metadata <- list(
    folds = folds,
    fold_bases = lapply(folds_internal, `[[`, "basis"),
    aligned_bases = lapply(folds_internal, `[[`, "basis_aligned"),
    rotations = lapply(folds_internal, `[[`, "rotation"),
    fold_alphas = fold_alphas,
    ridge = ridge,
    procrustes = if (align) list(alignment = fold_info$alignment, consensus = fold_info$consensus) else NULL
  )

  list(
    values = values,
    method = "kfold",
    contrasts = contrast_list,
    metadata = metadata
  )
}

#' Print method for dkge_folds
#'
#' @param x A dkge_folds object
#' @param ... Additional arguments (unused)
#' @export
print.dkge_folds <- function(x, ...) {
  cat("DKGE Fold Definition\n")
  cat("--------------------\n")
  cat(sprintf("Type: %s\n", x$type))
  cat(sprintf("Folds: %d\n", x$k))

  sizes <- vapply(x$assignments, length, integer(1))
  cat(sprintf("Fold sizes: %s\n", paste(sizes, collapse = ", ")))

  if (x$type == "subject" && !is.null(x$metadata$n_subjects)) {
    coverage <- length(unique(unlist(x$assignments)))
    cat(sprintf("Subject coverage: %d/%d\n", coverage, x$metadata$n_subjects))
  }

  invisible(x)
}
