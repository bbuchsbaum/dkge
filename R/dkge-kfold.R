# dkge-kfold.R
# K-fold cross-fitting for DKGE contrasts

#' Define folds for K-fold cross-validation
#'
#' Create subject-level random folds or validated custom held-out sets.
#'
#' @param fit A `dkge` object or `dkge_data` bundle
#' @param type Type of fold definition:
#'   - `"subject"`: Random assignment of subjects to folds (default)
#'   - `"custom"`: User provides fold assignments directly
#' @param k Number of folds (ignored for type="custom")
#' @param assignments For type="custom", a list of fold assignments
#' @param seed Random seed for reproducible fold assignment
#' @param align Logical; if TRUE (default) compute Procrustes alignment/consensus when folds are evaluated.
#' @param partition Contract for custom assessment sets. `"exact"` (default)
#'   requires a nonoverlapping partition covering every subject. `"repeated"`
#'   permits subjects in more than one assessment set but still requires full
#'   coverage. `"partial"` permits incomplete coverage but remains
#'   nonoverlapping. Defining repeated folds does not imply that every consumer
#'   can aggregate repeated assessments: subject-collapsing contrast and
#'   classification routines reject them rather than silently overwriting or
#'   unequally weighting subjects.
#' @param ... Additional arguments for specific fold types
#'
#' @return A `dkge_folds` object containing:
#'   - `type`: The fold type used
#'   - `k`: Number of folds
#'   - `assignments`: List specifying which data belongs to each fold
#'   - `metadata`: Additional information about fold creation
#'
#' @details
#' **Subject-level folds** (`type = "subject"`): Assigns entire subjects to folds.
#' This maintains subject independence and is appropriate when subjects are
#' exchangeable. Each fold will have approximately S/k subjects.
#'
#' **Custom folds** (`type = "custom"`): Full control over fold assignments.
#' Supply a list where each element specifies held-out subject indices. The
#' default exact-partition contract prevents accidental overlap or omission.
#'
#' @examples
#' \donttest{
#' toy <- dkge_sim_toy(
#'   factors = list(A = list(L = 2), B = list(L = 3)),
#'   active_terms = c("A", "B"), S = 6, P = 20, snr = 5
#' )
#' fit <- dkge(toy$B_list, toy$X_list, kernel = toy$K, rank = 2)
#'
#' # Random subject-level 5-fold
#' folds <- dkge_define_folds(fit, type = "subject", k = 5)
#'
#' # Custom fold specification (explicit holdouts)
#' folds_custom <- dkge_define_folds(
#'   fit,
#'   type = "custom",
#'   assignments = list(fold1 = 1:2, fold2 = 3:4, fold3 = 5:6)
#' )
#' folds_custom$k
#' }
#'
#' @export
dkge_define_folds <- function(fit, type = c("subject", "custom"),
                             k = 5, assignments = NULL,
                             seed = NULL, align = FALSE,
                             partition = c("exact", "repeated", "partial"),
                             ...) {
  type <- match.arg(type)
  partition <- match.arg(partition)

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

  if (!is.null(seed)) {
    seed <- .dkge_validate_integer_scalar(seed, "seed")
  }
  seed_state <- .dkge_seed_enter(seed)
  on.exit(.dkge_seed_exit(seed_state), add = TRUE)

  folds <- switch(type,
    subject = .define_subject_folds(n_subjects, k, subject_ids, seed = seed),
    custom = .validate_custom_folds(assignments, n_subjects, partition)
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
  k <- .dkge_validate_positive_integer(k, "k")
  if (k < 2L || k > n_subjects) {
    .dkge_abort(
      sprintf("`k` supplied value %d; expected an integer from 2 through %d.",
              k, n_subjects),
      "dkge_fold_partition_error"
    )
  }

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

#' Validate custom fold assignments
#' @keywords internal
#' @noRd
.validate_custom_folds <- function(assignments, n_subjects,
                                   partition = c("exact", "repeated", "partial")) {
  partition <- match.arg(partition)
  abort <- function(message) {
    .dkge_abort(message, "dkge_fold_partition_error")
  }
  if (is.null(assignments)) {
    abort("`assignments` is required for `type = \"custom\"`.")
  }
  if (!is.list(assignments) || length(assignments) < 2L) {
    abort("Custom folds require a list containing at least two assessment sets.")
  }

  assignment_names <- names(assignments)
  assignments <- lapply(seq_along(assignments), function(i) {
    idx <- assignments[[i]]
    if (!is.numeric(idx) || !length(idx) || any(!is.finite(idx)) ||
        any(idx != trunc(idx))) {
      abort(sprintf("Custom fold %d must contain finite integer subject indices.", i))
    }
    idx <- as.integer(idx)
    if (anyDuplicated(idx)) {
      abort(sprintf("Custom fold %d repeats a subject within one assessment set.", i))
    }
    if (any(idx < 1L) || any(idx > n_subjects)) {
      abort(sprintf("Custom fold %d contains indices outside 1 through %d.",
                    i, n_subjects))
    }
    sort(idx)
  })
  if (!is.null(assignment_names)) names(assignments) <- assignment_names

  all_idx <- unlist(assignments, use.names = FALSE)
  has_overlap <- anyDuplicated(all_idx) > 0L
  covered <- sort(unique(all_idx))
  full_coverage <- identical(covered, seq_len(n_subjects))

  if (partition != "repeated" && has_overlap) {
    abort(
      "Custom folds must form a nonoverlapping partition; use `partition = \"repeated\"` to allow repeated assessment."
    )
  }
  if (partition != "partial" && !full_coverage) {
    abort(
      "Custom folds must cover every subject; use `partition = \"partial\"` for incomplete assessment."
    )
  }

  list(
    type = "custom",
    k = length(assignments),
    assignments = assignments,
    metadata = list(
      n_subjects = n_subjects,
      coverage = length(covered),
      partition = partition,
      overlap = has_overlap
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
                                parallel, verbose, align = FALSE,
                                missingness = c("none", "rescale", "mask", "shrink"),
                                miss_args = list(), ...) {
  # Prepare folds
  missingness <- match.arg(missingness)
  fold_info_raw <- .dkge_normalize_folds(
    folds, fit, consumer = "K-fold contrasts"
  )
  folds <- fold_info_raw$folds

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
    verbose = verbose,
    missingness = missingness,
    miss_args = miss_args
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
    missingness = missingness,
    miss_args = miss_args,
    pair_counts = lapply(folds_internal, `[[`, "pair_counts"),
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
#' @return `x`, invisibly.
#' @method print dkge_folds
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
