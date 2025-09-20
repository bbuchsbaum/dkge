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
#' @param parallel Use parallel processing
#' @param verbose Print progress
#' @param ... Additional arguments
#' @return List with values, metadata, etc.
#' @keywords internal
#' @export
.dkge_contrast_kfold <- function(fit, contrast_list, folds, ridge,
                                parallel, verbose, ...) {
  # Prepare folds
  if (is.numeric(folds) && length(folds) == 1) {
    folds <- dkge_define_folds(fit, type = "subject", k = folds)
  } else if (!inherits(folds, "dkge_folds")) {
    stop("folds must be an integer k or a dkge_folds object")
  }

  S <- length(fit$Btil)
  q <- nrow(fit$U)
  r <- ncol(fit$U)
  n_contrasts <- length(contrast_list)
  k <- folds$k

  if (verbose) {
    message(sprintf("Computing %d contrast(s) via %d-fold CV", n_contrasts, k))
  }

  # Storage
  values <- vector("list", n_contrasts)
  names(values) <- names(contrast_list)
  fold_bases <- vector("list", k)
  fold_alphas <- vector("list", n_contrasts)

  # Process each fold
  for (fold_idx in seq_len(k)) {
    fold_subjects <- folds$assignments[[fold_idx]]
    if (verbose) {
      message(sprintf("Processing fold %d/%d (%d subjects)", fold_idx, k, length(fold_subjects)))
    }

    # Rebuild Chat excluding this fold
    Chat_minus_fold <- fit$Chat
    for (s in fold_subjects) {
      Chat_minus_fold <- Chat_minus_fold - fit$weights[s] * fit$contribs[[s]]
    }

    if (ridge > 0) {
      Chat_minus_fold <- Chat_minus_fold + ridge * diag(q)
    }
    Chat_minus_fold <- (Chat_minus_fold + t(Chat_minus_fold)) / 2

    # Eigen decomposition for fold-specific basis
    eig_fold <- eigen(Chat_minus_fold, symmetric = TRUE)
    U_fold <- fit$Kihalf %*% eig_fold$vectors[, seq_len(r), drop = FALSE]
    fold_bases[[fold_idx]] <- U_fold

    # Project each contrast
    for (i in seq_along(contrast_list)) {
      c_tilde <- backsolve(fit$R, contrast_list[[i]], transpose = FALSE)
      alpha_fold <- t(U_fold) %*% fit$K %*% c_tilde

      # Initialize storage on first contrast
      if (fold_idx == 1) {
        values[[i]] <- vector("list", S)
        fold_alphas[[i]] <- matrix(NA_real_, k, r)
      }

      fold_alphas[[i]][fold_idx, ] <- alpha_fold

      # Compute values for held-out subjects
      for (s in fold_subjects) {
        Bts <- fit$Btil[[s]]
        A_s <- t(Bts) %*% fit$K %*% U_fold
        v_s <- as.numeric(A_s %*% alpha_fold)
        values[[i]][[s]] <- v_s
      }
    }
  }

  # Verify all subjects were processed
  for (i in seq_along(values)) {
    missing <- which(vapply(values[[i]], is.null, logical(1)))
    if (length(missing) > 0) {
      warning(sprintf("Subjects %s not in any fold for contrast %s",
                     paste(missing, collapse = ","), names(values)[i]))
    }
  }

  list(
    values = values,
    method = "kfold",
    contrasts = contrast_list,
    metadata = list(
      folds = folds,
      fold_bases = fold_bases,
      fold_alphas = fold_alphas,
      ridge = ridge
    )
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
