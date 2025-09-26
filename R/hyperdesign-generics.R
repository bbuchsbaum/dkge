# hyperdesign-generics.R
# Minimal S3 generics to accept hyperdesign inputs without touching core solvers.

#' Convert to a DKGE design kernel
#'
#' Coerce arbitrary objects into a DKGE kernel bundle. The default method
#' preserves existing behaviour by accepting matrices or list objects with a
#' `$K` component and optional metadata.
#'
#' @param x Object containing kernel information
#' @param ... Additional arguments passed to methods
#' @return List with entries `K` (q x q matrix) and optional `info`
#' @seealso `vignette("dkge-classification", package = "dkge")` for a worked
#'   example that uses these generics to integrate hyperdesign inputs.
#' @export
as_dkge_kernel <- function(x, ...) {
  UseMethod("as_dkge_kernel")
}

#' @export
as_dkge_kernel.default <- function(x, ...) {
  if (is.matrix(x)) {
    return(list(K = x, info = NULL))
  }
  if (is.list(x) && !is.null(x$K)) {
    return(x)
  }
  stop(
    "Don't know how to convert object of class '",
    paste(class(x), collapse = "/"),
    "' into a DKGE kernel.",
    call. = FALSE
  )
}

#' Convert to DKGE fold assignments
#'
#' Coerce various fold specifications into the `dkge_folds` structure used by
#' dkge cross-fitting helpers. Methods may use `fit_or_data` to resolve subject
#' identifiers when necessary.
#'
#' @param x Object describing fold assignments
#' @param fit_or_data Optional `dkge` or `dkge_data` object used to resolve
#'   subject identifiers
#' @param ... Additional arguments passed to methods
#' @return Object with class `dkge_folds`
#' @seealso `vignette("dkge-classification", package = "dkge")` for an example
#'   of fold conversion in practice.
#' @export
as_dkge_folds <- function(x, fit_or_data = NULL, ...) {
  UseMethod("as_dkge_folds")
}

#' @export
as_dkge_folds.dkge_folds <- function(x, fit_or_data = NULL, ...) {
  x
}

#' @export
as_dkge_folds.default <- function(x, fit_or_data = NULL, ...) {
  subject_ids <- .dkge_resolve_subject_ids(fit_or_data)

  assignments <- .dkge_fold_assignments_from_input(x, subject_ids)

  structure(
    list(
      type = "custom",
      k = length(assignments),
      assignments = assignments,
      subject_ids = subject_ids,
      metadata = list(
        n_subjects = length(subject_ids) %||% NA_integer_,
        coverage = length(unique(unlist(assignments))),
        source = class(x)
      ),
      align = FALSE
    ),
    class = "dkge_folds"
  )
}

# Helpers -----------------------------------------------------------------

.dkge_resolve_subject_ids <- function(fit_or_data) {
  if (inherits(fit_or_data, "dkge")) {
    return(fit_or_data$subject_ids)
  }
  if (inherits(fit_or_data, "dkge_data")) {
    return(fit_or_data$subject_ids)
  }
  NULL
}

.dkge_validate_assignments <- function(assignments, n_subjects = NULL) {
  stopifnot(is.list(assignments), length(assignments) >= 2)
  all_idx <- unlist(assignments, use.names = FALSE)
  stopifnot(length(all_idx) > 0)
  if (!is.null(n_subjects)) {
    stopifnot(all(all_idx >= 1L), all(all_idx <= n_subjects))
  }
  assignments
}

.dkge_fold_assignments_from_input <- function(x, subject_ids) {
  n_subjects <- length(subject_ids)

  if (is.list(x) && all(vapply(x, is.numeric, logical(1)))) {
    assignments <- lapply(x, as.integer)
    return(.dkge_validate_assignments(assignments, n_subjects = n_subjects))
  }

  if (is.data.frame(x) && all(c("subject", "fold") %in% names(x))) {
    if (is.null(subject_ids)) {
      stop("fit_or_data with subject_ids required for data.frame input.", call. = FALSE)
    }
    subj <- as.character(x$subject)
    fold <- x$fold
    idx <- match(subj, subject_ids)
    if (anyNA(idx)) {
      stop("Some subjects in 'x' were not found among dkge subject_ids.", call. = FALSE)
    }
    split(idx, fold)
  } else if (is.vector(x) && !is.list(x)) {
    labs <- x
    if (!is.null(names(labs))) {
      if (is.null(subject_ids)) {
        stop("fit_or_data with subject_ids required for named vector input.", call. = FALSE)
      }
      idx <- match(subject_ids, names(labs))
      if (anyNA(idx)) {
        stop("Some dkge subject_ids were not present in names(x).", call. = FALSE)
      }
      split(seq_along(subject_ids), labs[idx])
    } else {
      if (is.null(subject_ids)) {
        stop("fit_or_data with subject_ids required for unnamed vector input.", call. = FALSE)
      }
      if (length(labs) != n_subjects) {
        stop("Length of x must equal number of subjects.", call. = FALSE)
      }
      split(seq_len(n_subjects), labs)
    }
  } else if (inherits(x, "dkge_fold_over")) {
    stop("Objects of class 'dkge_fold_over' are not supported by default coercion.")
  } else {
    stop(
      "Unsupported folds specification of class '",
      paste(class(x), collapse = "/"),
      "'.",
      call. = FALSE
    )
  }
}
