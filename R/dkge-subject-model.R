# dkge-subject-model.R
# Formula design objects for between-subject DKGE models.

#' Build a between-subject model matrix
#'
#' Wraps formula handling for subject-level DKGE analyses. The result owns the
#' model matrix, term metadata, and subject identifiers used by
#' [dkge_between_rrr()].
#'
#' @param formula A model formula such as `~ group * trait + age + sex`.
#' @param data Data frame containing subject-level variables.
#' @param subject_ids Optional subject identifiers aligned with `data` rows.
#'   Defaults to the `subject_id` column when present, then `rownames(data)`,
#'   then sequential IDs. Identifiers are stored as `rownames` of `data` before
#'   [stats::model.frame()] so `na.action = stats::na.omit` keeps IDs in sync
#'   with the retained rows.
#' @param nuisance Optional character vector of term labels to treat as nuisance
#'   in downstream inference. When `dkge_between_permute(terms = NULL)`, these
#'   terms are excluded from the default test set.
#' @param contrasts.arg Optional contrasts passed to [stats::model.matrix()].
#' @param na.action NA handler for the model frame. Defaults to [stats::na.fail()].
#'
#' @return Object of class `dkge_subject_model`.
#' @examples
#' dat <- data.frame(
#'   subject_id = paste0("s", 1:8),
#'   group = factor(rep(c("A", "B"), each = 4)),
#'   trait = rnorm(8)
#' )
#' dkge_subject_model(~ group * trait, dat)
#' @export
dkge_subject_model <- function(formula,
                               data,
                               subject_ids = NULL,
                               nuisance = NULL,
                               contrasts.arg = NULL,
                               na.action = stats::na.fail) {
  if (!inherits(formula, "formula")) {
    stop("`formula` must be a formula.", call. = FALSE)
  }
  data <- as.data.frame(data)
  n_data <- nrow(data)
  if (is.null(subject_ids)) {
    if ("subject_id" %in% names(data)) {
      ids <- as.character(data$subject_id)
    } else if (!is.null(rownames(data)) && all(nzchar(rownames(data)))) {
      ids <- rownames(data)
    } else {
      ids <- paste0("subj", seq_len(n_data))
    }
  } else {
    ids <- as.character(subject_ids)
    if (length(ids) != n_data) {
      stop("`subject_ids` must be unique, non-empty, and match the data rows.",
           call. = FALSE)
    }
  }
  if (length(ids) != n_data || any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop("`subject_ids` must be unique, non-empty, and match the data rows.",
         call. = FALSE)
  }
  rownames(data) <- ids
  mf <- stats::model.frame(formula, data = data, na.action = na.action)
  tt <- stats::terms(mf)
  X <- stats::model.matrix(tt, data = mf, contrasts.arg = contrasts.arg)
  subject_ids <- rownames(mf)
  if (length(subject_ids) != nrow(X) || any(!nzchar(subject_ids)) ||
      anyDuplicated(subject_ids)) {
    stop("`subject_ids` must be unique, non-empty, and match the model rows.",
         call. = FALSE)
  }
  rownames(X) <- subject_ids

  qrX <- qr(X)
  if (qrX$rank < ncol(X)) {
    stop("Between-subject model matrix is rank deficient; revise the formula or contrasts.",
         call. = FALSE)
  }

  assign <- attr(X, "assign")
  term_labels <- attr(tt, "term.labels")
  term_names <- c("(Intercept)", term_labels)
  term_columns <- stats::setNames(vector("list", length(term_names)), term_names)
  for (i in seq_along(term_names)) {
    term_columns[[i]] <- which(assign == (i - 1L))
  }

  if (!is.null(nuisance)) {
    nuisance <- as.character(nuisance)
    missing <- setdiff(nuisance, term_names)
    if (length(missing)) {
      stop("Nuisance terms not present in the model: ",
           paste(missing, collapse = ", "), call. = FALSE)
    }
  }

  out <- list(
    X = X,
    formula = formula,
    data = data,
    model_frame = mf,
    terms = tt,
    term_labels = term_labels,
    term_names = term_names,
    term_columns = term_columns,
    assign = assign,
    contrasts = attr(X, "contrasts"),
    nuisance = nuisance,
    subject_ids = subject_ids
  )
  class(out) <- c("dkge_subject_model", "list")
  out
}

#' Print a subject-level model specification
#'
#' @param x A `dkge_subject_model` object.
#' @param ... Unused; present for S3 method compatibility.
#' @return `x`, invisibly.
#' @method print dkge_subject_model
#' @export
print.dkge_subject_model <- function(x, ...) {
  cat("<dkge_subject_model>", "\n", sep = "")
  cat("  subjects :", nrow(x$X), "\n")
  cat("  columns  :", ncol(x$X), "\n")
  cat("  terms    :", paste(x$term_names, collapse = ", "), "\n")
  invisible(x)
}
