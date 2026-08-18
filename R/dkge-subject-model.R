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
#' @param subject_ids Optional subject identifiers. Defaults to a `subject_id`
#'   column when present, then row names, then sequential IDs.
#' @param subject_id_col Optional column name in `data` to use as subject
#'   identifiers when `subject_ids` is not supplied. Defaults to `"subject_id"`.
#' @param nuisance Optional character vector of term labels to treat as nuisance
#'   in downstream inference.
#' @param contrasts.arg Optional contrasts passed to [stats::model.matrix()].
#' @param na.action NA handler for the model frame. Defaults to [stats::na.fail()].
#'
#' @return Object of class `dkge_subject_model`.
#' @export
dkge_subject_model <- function(formula,
                               data,
                               subject_ids = NULL,
                               nuisance = NULL,
                               contrasts.arg = NULL,
                               na.action = stats::na.fail,
                               subject_id_col = "subject_id") {
  if (!inherits(formula, "formula")) {
    stop("`formula` must be a formula.", call. = FALSE)
  }
  data <- as.data.frame(data)
  mf <- stats::model.frame(formula, data = data, na.action = na.action)
  tt <- stats::terms(mf)
  X <- stats::model.matrix(tt, data = mf, contrasts.arg = contrasts.arg)

  if (is.null(subject_ids)) {
    if (!is.null(subject_id_col)) {
      if (!is.character(subject_id_col) || length(subject_id_col) != 1L || !nzchar(subject_id_col)) {
        stop("`subject_id_col` must be NULL or a non-empty column name.", call. = FALSE)
      }
      if (subject_id_col %in% names(data)) {
        subject_ids <- as.character(data[[subject_id_col]])
      } else if (!identical(subject_id_col, "subject_id")) {
        stop("`subject_id_col` was not found in `data`: ", subject_id_col, call. = FALSE)
      }
    }
    if (is.null(subject_ids) && !is.null(rownames(data)) && all(nzchar(rownames(data)))) {
      subject_ids <- rownames(data)
    } else if (is.null(subject_ids)) {
      subject_ids <- paste0("subj", seq_len(nrow(X)))
    }
  }
  subject_ids <- as.character(subject_ids)
  if (length(subject_ids) != nrow(X) || any(!nzchar(subject_ids)) || any(duplicated(subject_ids))) {
    stop("`subject_ids` must be unique, non-empty, and match the model rows.", call. = FALSE)
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

#' @export
print.dkge_subject_model <- function(x, ...) {
  cat("<dkge_subject_model>", "\n", sep = "")
  cat("  subjects :", nrow(x$X), "\n")
  cat("  columns  :", ncol(x$X), "\n")
  cat("  terms    :", paste(x$term_names, collapse = ", "), "\n")
  invisible(x)
}
