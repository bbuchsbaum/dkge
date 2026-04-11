# dkge-between-rrr.R
# Reduced-rank regression for between-subject DKGE targets.

#' Reduced-rank regression on a DKGE subject target
#'
#' Fits a directional multivariate second-level model
#' \deqn{Y = X B + E,\quad rank(B) <= r}
#' using a compact SVD of the fitted design-space response. This is intended for
#' subject-level factors, traits, and covariates after DKGE has produced a
#' subject-by-feature target.
#'
#' @param target A `dkge_target` from [dkge_make_target()] or a numeric matrix.
#' @param design A `dkge_subject_model` from [dkge_subject_model()] or a numeric
#'   model matrix.
#' @param rank Reduced rank. Defaults to the maximum estimable rank.
#' @param weights `"none"`, `"target"`, or a list with optional `subject` and
#'   `feature` numeric weights.
#' @param feature_mask Optional logical vector selecting target features.
#' @param tol Numerical tolerance for rank checks.
#'
#' @return Object of class `dkge_between_rrr`.
#' @export
dkge_between_rrr <- function(target,
                             design,
                             rank = NULL,
                             weights = c("none", "target"),
                             feature_mask = NULL,
                             tol = 1e-8) {
  target <- .dkge_as_target(target)
  design <- .dkge_as_subject_model(design, target$subject_ids)
  if (!is.list(weights)) {
    weights <- match.arg(weights)
  }

  ord <- match(target$subject_ids, design$subject_ids)
  if (anyNA(ord)) {
    stop("All target subjects must be present in the between-subject design.",
         call. = FALSE)
  }

  X <- design$X[ord, , drop = FALSE]
  Y <- target$Y

  mask <- target$mask %||% rep(TRUE, ncol(Y))
  if (!is.null(feature_mask)) {
    feature_mask <- as.logical(feature_mask)
    if (length(feature_mask) != ncol(Y)) {
      stop("`feature_mask` must have length ncol(target$Y).", call. = FALSE)
    }
    mask <- mask & feature_mask
  }
  if (!any(mask)) {
    stop("No target features remain after applying the feature mask.", call. = FALSE)
  }
  Y <- Y[, mask, drop = FALSE]
  feature_ids <- target$feature_ids[mask]

  weight_payload <- .dkge_between_weights(weights, target, mask)
  subject_w <- weight_payload$subject
  feature_w <- weight_payload$feature

  core <- .dkge_rrr_fit_core(X = X,
                             Y = Y,
                             rank = rank,
                             subject_weights = subject_w,
                             feature_weights = feature_w,
                             coef_ids = colnames(X),
                             feature_ids = feature_ids,
                             subject_ids = target$subject_ids,
                             tol = tol)

  out <- list(
    coef = core$coef,
    fitted = core$fitted,
    residuals = core$residuals,
    scores_subject = core$scores_subject,
    loadings_brain = core$loadings_brain,
    loadings_design = core$loadings_design,
    design_saliences = core$design_saliences,
    singular_values = core$singular_values,
    rank = core$rank,
    max_rank = core$max_rank,
    target = target,
    design = design,
    target_mask = mask,
    subject_weights = subject_w,
    feature_weights = feature_w,
    call = match.call()
  )
  class(out) <- c("dkge_between_rrr", "list")
  out
}

.dkge_rrr_fit_core <- function(X,
                               Y,
                               rank = NULL,
                               subject_weights = NULL,
                               feature_weights = NULL,
                               coef_ids = colnames(X),
                               feature_ids = colnames(Y),
                               subject_ids = rownames(Y),
                               tol = 1e-8,
  warn_rank = TRUE) {
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  if (!is.numeric(X) || !is.numeric(Y) || any(!is.finite(X)) || any(!is.finite(Y))) {
    stop("RRR inputs `X` and `Y` must be numeric and finite.", call. = FALSE)
  }
  if (nrow(X) != nrow(Y)) {
    stop("RRR inputs `X` and `Y` must have the same number of rows.", call. = FALSE)
  }
  if (ncol(X) < 1L || ncol(Y) < 1L) {
    stop("RRR inputs `X` and `Y` must have at least one column.", call. = FALSE)
  }
  Xw <- X
  Yw <- Y
  if (!is.null(subject_weights)) {
    subject_weights <- as.numeric(subject_weights)
    if (length(subject_weights) != nrow(X) || any(!is.finite(subject_weights)) || any(subject_weights < 0)) {
      stop("Subject weights must be non-negative, finite, and match input rows.",
           call. = FALSE)
    }
    sw <- sqrt(pmax(subject_weights, 0))
    Xw <- Xw * sw
    Yw <- Yw * sw
  }
  if (!is.null(feature_weights)) {
    feature_weights <- as.numeric(feature_weights)
    if (length(feature_weights) != ncol(Y) || any(!is.finite(feature_weights)) || any(feature_weights <= 0)) {
      stop("Feature weights must be positive, finite, and match target columns.",
           call. = FALSE)
    }
    fw <- sqrt(feature_weights)
    Yw <- sweep(Yw, 2L, fw, "*")
  } else {
    fw <- rep(1, ncol(Yw))
  }

  qrX <- qr(Xw, tol = tol)
  if (qrX$rank < ncol(Xw)) {
    stop("Weighted between-subject model matrix is rank deficient.", call. = FALSE)
  }
  Q <- qr.Q(qrX)
  R <- qr.R(qrX)

  C <- crossprod(Q, Yw)
  max_rank <- min(nrow(C), ncol(C))
  if (is.null(rank)) {
    rank <- max_rank
  }
  rank <- as.integer(rank)
  if (length(rank) != 1L || is.na(rank) || rank < 1L) {
    stop("`rank` must be a positive integer.", call. = FALSE)
  }
  if (rank > max_rank) {
    if (warn_rank) {
      warning("Requested rank exceeds the estimable rank; using ", max_rank, ".",
              call. = FALSE)
    }
    rank <- max_rank
  }

  sv <- svd(C, nu = rank, nv = rank)
  U <- sv$u[, seq_len(rank), drop = FALSE]
  V <- sv$v[, seq_len(rank), drop = FALSE]
  d <- sv$d[seq_len(rank)]
  D <- diag(d, nrow = rank)

  A <- backsolve(R, U)
  coef_weighted <- A %*% D %*% t(V)
  coef <- sweep(coef_weighted, 2L, fw, "/")
  rownames(coef) <- coef_ids %||% paste0("x", seq_len(nrow(coef)))
  colnames(coef) <- feature_ids %||% paste0("feature", seq_len(ncol(coef)))

  fitted <- X %*% coef
  residuals <- Y - fitted
  scores_subject <- X %*% A
  colnames(scores_subject) <- paste0("LV", seq_len(rank))
  rownames(scores_subject) <- subject_ids
  colnames(V) <- paste0("LV", seq_len(rank))
  rownames(V) <- colnames(coef)
  colnames(A) <- paste0("LV", seq_len(rank))
  rownames(A) <- rownames(coef)

  list(
    coef = coef,
    fitted = fitted,
    residuals = residuals,
    scores_subject = scores_subject,
    loadings_brain = V,
    loadings_design = A,
    design_saliences = U,
    singular_values = d,
    rank = rank,
    max_rank = max_rank
  )
}

.dkge_as_target <- function(target) {
  if (inherits(target, "dkge_target")) return(target)
  if (is.matrix(target) || is.data.frame(target)) {
    return(dkge_make_target(Y = as.matrix(target)))
  }
  stop("`target` must be a dkge_target or numeric matrix.", call. = FALSE)
}

.dkge_as_subject_model <- function(design, subject_ids) {
  if (inherits(design, "dkge_subject_model")) return(design)
  if (is.matrix(design) || is.data.frame(design)) {
    X <- as.matrix(design)
    if (!is.numeric(X) || any(!is.finite(X))) {
      stop("Design matrix must be numeric and finite.", call. = FALSE)
    }
    if (is.null(rownames(X))) {
      rownames(X) <- subject_ids %||% paste0("subj", seq_len(nrow(X)))
    }
    if (is.null(colnames(X))) {
      colnames(X) <- paste0("x", seq_len(ncol(X)))
    } else {
      blank <- !nzchar(colnames(X))
      colnames(X)[blank] <- paste0("x", which(blank))
    }
    if (length(rownames(X)) != nrow(X) || any(!nzchar(rownames(X))) || any(duplicated(rownames(X)))) {
      stop("Design matrix row names must be unique and non-empty.", call. = FALSE)
    }
    if (length(colnames(X)) != ncol(X) || any(!nzchar(colnames(X))) || any(duplicated(colnames(X)))) {
      stop("Design matrix column names must be unique and non-empty.", call. = FALSE)
    }
    term_columns <- stats::setNames(as.list(seq_len(ncol(X))), colnames(X))
    out <- list(
      X = X,
      formula = NULL,
      data = NULL,
      model_frame = NULL,
      terms = NULL,
      term_labels = colnames(X),
      term_names = colnames(X),
      term_columns = term_columns,
      assign = seq_len(ncol(X)),
      contrasts = NULL,
      nuisance = NULL,
      subject_ids = rownames(X)
    )
    class(out) <- c("dkge_subject_model", "list")
    return(out)
  }
  stop("`design` must be a dkge_subject_model or numeric model matrix.", call. = FALSE)
}

.dkge_between_weights <- function(weights, target, mask) {
  if (is.list(weights)) {
    subject <- weights$subject %||% weights$subject_weights
    feature <- weights$feature %||% weights$feature_weights
  } else {
    weights <- as.character(weights)
    if (length(weights) != 1L || !weights %in% c("none", "target")) {
      stop("`weights` must be 'none', 'target', or a list of weights.", call. = FALSE)
    }
    if (identical(weights, "target")) {
      subject <- target$subject_weights
      feature <- target$feature_weights
      if (is.null(feature) && !is.null(target$precision)) {
        feature <- if (is.matrix(target$precision)) {
          colMeans(target$precision, na.rm = TRUE)
        } else {
          as.numeric(target$precision)
        }
      }
      if (is.null(feature) && !is.null(target$coverage)) {
        feature <- if (is.matrix(target$coverage)) {
          colMeans(target$coverage, na.rm = TRUE)
        } else {
          as.numeric(target$coverage)
        }
      }
    } else {
      subject <- NULL
      feature <- NULL
    }
  }
  if (!is.null(subject)) {
    subject <- as.numeric(subject)
    if (length(subject) != nrow(target$Y) || any(!is.finite(subject)) || any(subject < 0)) {
      stop("Subject weights must be non-negative, finite, and match target rows.",
           call. = FALSE)
    }
  }
  if (!is.null(feature)) {
    feature <- as.numeric(feature)
    if (length(feature) != ncol(target$Y) || any(!is.finite(feature)) || any(feature < 0)) {
      stop("Feature weights must be non-negative, finite, and match target columns.",
           call. = FALSE)
    }
    feature <- feature[mask]
    if (any(feature <= 0)) {
      stop("Selected feature weights must be strictly positive.", call. = FALSE)
    }
  }
  list(subject = subject, feature = feature)
}

.dkge_between_term_rows <- function(design, coef_matrix, term) {
  term <- as.character(term)
  if (term %in% rownames(coef_matrix)) {
    return(match(term, rownames(coef_matrix)))
  }
  cols <- design$term_columns[[term]]
  if (is.null(cols) || !length(cols)) {
    stop("Unknown model term or coefficient: ", term, call. = FALSE)
  }
  cols
}

.dkge_between_term_map_from_coef <- function(coef_matrix,
                                            design,
                                            term,
                                            contrast = NULL,
                                            drop = TRUE) {
  B <- as.matrix(coef_matrix)
  if (!is.null(contrast)) {
    contrast <- as.numeric(contrast)
    if (length(contrast) != nrow(B)) {
      stop("`contrast` must have length equal to the number of design columns.",
           call. = FALSE)
    }
    label <- if (missing(term) || is.null(term)) "contrast" else as.character(term)
    out <- matrix(drop(crossprod(contrast, B)), nrow = 1L,
                  dimnames = list(label, colnames(B)))
    return(if (drop) drop(out) else out)
  }
  rows <- .dkge_between_term_rows(design, B, term)
  out <- B[rows, , drop = FALSE]
  rownames(out) <- rownames(B)[rows]
  if (drop && nrow(out) == 1L) drop(out) else out
}

.dkge_between_new_design_matrix <- function(object, newdata) {
  if (inherits(newdata, "dkge_subject_model")) {
    X <- newdata$X
  } else if (is.data.frame(newdata) && !all(vapply(newdata, is.numeric, logical(1)))) {
    if (is.null(object$design$formula)) {
      stop("Cannot rebuild a formula design matrix because the fitted object does not retain a formula.",
           call. = FALSE)
    }
    mf <- stats::model.frame(object$design$formula,
                             data = newdata,
                             na.action = stats::na.fail)
    X <- stats::model.matrix(object$design$terms,
                             data = mf,
                             contrasts.arg = object$design$contrasts)
  } else if (is.matrix(newdata) || is.data.frame(newdata)) {
    X <- as.matrix(newdata)
  } else {
    stop("`newdata` must be a dkge_subject_model, data.frame, or numeric matrix.",
         call. = FALSE)
  }

  X <- as.matrix(X)
  expected <- rownames(object$coef)
  if (!all(expected %in% colnames(X))) {
    stop("`newdata` design columns must include: ",
         paste(expected, collapse = ", "),
         call. = FALSE)
  }
  X[, expected, drop = FALSE]
}

#' Extract a term-specific coefficient map
#'
#' @param object A `dkge_between_rrr` object.
#' @param term Formula term or model-matrix column name.
#' @param contrast Optional numeric contrast over model-matrix columns. When
#'   supplied, `term` is used only as a label.
#' @param drop Logical; drop single-row results to a vector.
#'
#' @return Numeric vector or matrix of coefficients in target feature space.
#' @export
dkge_term_map <- function(object, term, contrast = NULL, drop = TRUE) {
  stopifnot(inherits(object, "dkge_between_rrr"))
  if (missing(term) || is.null(term)) {
    stop("`term` is required when `contrast` is not supplied.", call. = FALSE)
  }
  .dkge_between_term_map_from_coef(object$coef,
                                   design = object$design,
                                   term = term,
                                   contrast = contrast,
                                   drop = drop)
}

#' @export
coef.dkge_between_rrr <- function(object, ...) {
  object$coef
}

#' @export
fitted.dkge_between_rrr <- function(object, ...) {
  object$fitted
}

#' @export
residuals.dkge_between_rrr <- function(object, ...) {
  object$residuals
}

#' @export
predict.dkge_between_rrr <- function(object,
                                     newdata = NULL,
                                     type = c("response", "scores"),
                                     ...) {
  type <- match.arg(type)
  if (is.null(newdata)) {
    return(switch(type,
                  response = object$fitted,
                  scores = object$scores_subject))
  }

  Xnew <- .dkge_between_new_design_matrix(object, newdata)
  switch(type,
         response = Xnew %*% object$coef,
         scores = Xnew %*% object$loadings_design)
}

#' @export
print.dkge_between_rrr <- function(x, ...) {
  cat("<dkge_between_rrr>", "\n", sep = "")
  cat("  subjects :", nrow(x$target$Y), "\n")
  cat("  features :", ncol(x$coef), "\n")
  cat("  rank     :", x$rank, "\n")
  cat("  terms    :", paste(x$design$term_names, collapse = ", "), "\n")
  invisible(x)
}
