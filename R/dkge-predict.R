
# dkge-predict.R (v0.5 add-on)
# Predict out-of-sample loadings and contrasts given a frozen basis U, K, R.

#' Freeze a DKGE fit into a compact model for prediction
#' @param fit a dkge or dkge_stream object
#' @return list with U, K, R and class 'dkge_model'
#' @export
dkge_freeze <- function(fit) {
  stopifnot(is.list(fit), !is.null(fit$U), !is.null(fit$K), !is.null(fit$R))
  model <- list(U = fit$U, K = fit$K, R = fit$R, effects = fit$effects)
  class(model) <- "dkge_model"
  model
}

.dkge_coerce_beta <- function(x) {
  if (inherits(x, "dkge_subject")) {
    return(as.matrix(x$beta))
  }
  as.matrix(x)
}

.dkge_components <- function(object) {
  if (inherits(object, "dkge_model")) {
    list(U = object$U, K = object$K, R = object$R, effects = object$effects)
  } else {
    list(U = object$U, K = object$K, R = object$R, effects = object$effects)
  }
}

.dkge_align_effects <- function(B, effects) {
  if (is.null(effects)) {
    return(B)
  }
  stopifnot(nrow(B) == length(effects))
  if (!is.null(rownames(B))) {
    idx <- match(effects, rownames(B))
    if (anyNA(idx)) {
      stop("Row names of new betas must include all training effects.")
    }
    B <- B[idx, , drop = FALSE]
  }
  B
}

#' Predict DKGE loadings for new subjects (out-of-sample)
#'
#' @param object dkge | dkge_stream | dkge_model
#' @param B_list list of qxP_s beta matrices for new subjects
#' @return list of P_sxr loadings (A_s) for each subject
#' @export
dkge_predict_loadings <- function(object, B_list) {
  comps <- .dkge_components(object)
  mats <- lapply(B_list, .dkge_coerce_beta)
  lapply(mats, function(Bs) {
    Bs <- .dkge_align_effects(Bs, comps$effects)
    Btil <- t(comps$R) %*% Bs
    project_cpp <- get0("dkge_project_loadings_cpp", mode = "function")
    if (is.function(project_cpp)) {
      project_cpp(Btil, comps$K, comps$U)
    } else {
      t(Btil) %*% comps$K %*% comps$U
    }
  })
}

#' Predict DKGE contrasts for new subjects (out-of-sample)
#'
#' @param object dkge | dkge_stream | dkge_model
#' @param B_list list of qxP_s betas
#' @param contrasts list of named q-vectors or a qxk matrix (columns are contrasts)
#' @param return_loadings logical; if TRUE also return A_list
#' @return list(A_list=..., values = list of per-contrast subject vectors)
#' @export
dkge_predict <- function(object, B_list, contrasts, return_loadings = TRUE) {
  comps <- .dkge_components(object)
  mats <- lapply(B_list, .dkge_coerce_beta)
  mats <- lapply(mats, .dkge_align_effects, effects = comps$effects)
  if (is.matrix(contrasts)) {
    Cmat <- contrasts; cn <- colnames(Cmat); if (is.null(cn)) cn <- paste0("c", seq_len(ncol(Cmat)))
    contrasts <- setNames(lapply(seq_len(ncol(Cmat)), function(j) Cmat[, j]), cn)
  }
  stopifnot(is.list(contrasts) && length(contrasts) > 0)
  if (is.null(names(contrasts)) || any(!nzchar(names(contrasts)))) {
    names(contrasts) <- paste0("c", seq_along(contrasts))
  }
  subj_names <- names(B_list)
  if (is.null(subj_names)) subj_names <- rep("", length(B_list))
  missing <- which(!nzchar(subj_names))
  if (length(missing)) {
    candidate <- vapply(B_list[missing], function(x) {
      if (inherits(x, "dkge_subject")) {
        id <- x$id
        if (is.null(id) || !nzchar(id)) "" else id
      } else ""
    }, character(1))
    subj_names[missing] <- candidate
  }
  if (any(!nzchar(subj_names))) {
    fallback <- paste0("subj", seq_along(B_list))
    subj_names[!nzchar(subj_names)] <- fallback[!nzchar(subj_names)]
  }
  A_list <- dkge_predict_loadings(comps, mats)
  names(A_list) <- subj_names
  # precompute alpha per contrast
  alpha_list <- lapply(contrasts, function(c) {
    ctil <- backsolve(comps$R, c, transpose = FALSE)
    alpha_cpp <- get0("dkge_alpha_cpp", mode = "function")
    if (is.function(alpha_cpp)) {
      alpha_cpp(comps$U, comps$K, comps$R, c)
    } else {
      t(comps$U) %*% comps$K %*% ctil
    }
  })
  vals <- lapply(seq_along(mats), function(i) {
    A <- A_list[[i]]
    res <- sapply(alpha_list, function(a) as.numeric(A %*% a))
    if (is.matrix(res) && length(alpha_list) > 1) {
      colnames(res) <- names(alpha_list)
    } else {
      names(res) <- names(alpha_list)
    }
    res
  })
  names(vals) <- subj_names
  out <- list(values = vals)
  if (return_loadings) out$A_list <- A_list
  out
}


#' Convenience prediction for subject collections
#'
#' Harmonises a variety of beta inputs (matrices, `dkge_subject` objects, or
#' `dkge_data` bundles) before forwarding to [dkge_predict()]. This allows
#' callers to work with tidy inputs without manually assembling `B_list`
#' structures.
#'
#' @param object dkge | dkge_stream | dkge_model.
#' @param betas Subject data. Accepts a matrix, list of matrices,
#'   `dkge_subject` objects, or a `dkge_data` bundle.
#' @param contrasts List or matrix accepted by [dkge_predict()].
#' @param ids Optional subject identifiers overriding those inferred from
#'   `betas`.
#' @param return_loadings Logical; when TRUE, include projected loadings in the
#'   result bundle.
#' @return Output from [dkge_predict()] with harmonised subject names.
#' @export
dkge_predict_subjects <- function(object,
                                  betas,
                                  contrasts,
                                  ids = NULL,
                                  return_loadings = TRUE) {
  prep <- .dkge_prepare_predict_inputs(betas, ids)
  dkge_predict(object,
               B_list = prep$B_list,
               contrasts = contrasts,
               return_loadings = return_loadings)
}

#' @keywords internal
#' @noRd
.dkge_prepare_predict_inputs <- function(betas, ids = NULL) {
  B_list <- NULL
  derived_ids <- character(0)

  add_subject <- function(beta_obj, label = "") {
    if (inherits(beta_obj, "dkge_subject")) {
      mat <- as.matrix(beta_obj$beta)
      lid <- beta_obj$id %||% label %||% ""
    } else {
      mat <- as.matrix(beta_obj)
      lid <- label %||% ""
    }
    list(matrix = mat, id = lid)
  }

  if (inherits(betas, "dkge_data")) {
    B_list <- lapply(betas$betas, as.matrix)
    derived_ids <- betas$subject_ids %||% rep("", length(B_list))
  } else if (inherits(betas, "dkge_subject")) {
    subj <- add_subject(betas)
    B_list <- list(subj$matrix)
    derived_ids <- subj$id
  } else if (is.matrix(betas)) {
    B_list <- list(as.matrix(betas))
    derived_ids <- ""
  } else if (is.list(betas)) {
    B_list <- vector("list", length(betas))
    derived_ids <- rep("", length(betas))
    betas_names <- names(betas)
    for (i in seq_along(betas)) {
      label <- if (!is.null(betas_names) && nzchar(betas_names[i])) betas_names[i] else ""
      subj <- add_subject(betas[[i]], label = label)
      B_list[[i]] <- subj$matrix
      derived_ids[i] <- subj$id %||% label
    }
  } else {
    stop("Unsupported `betas` input; supply matrices, dkge_subject objects, or dkge_data", call. = FALSE)
  }

  if (is.null(B_list) || !length(B_list)) {
    stop("No subjects detected in `betas`", call. = FALSE)
  }

  if (!is.null(ids)) {
    if (length(ids) != length(B_list)) {
      stop("Length of `ids` must match the number of subjects inferred from `betas`.")
    }
    final_ids <- as.character(ids)
  } else {
    final_ids <- derived_ids
  }
  if (length(final_ids) != length(B_list)) {
    final_ids <- rep("", length(B_list))
  }
  missing <- which(!nzchar(final_ids))
  if (length(missing)) {
    fallback <- paste0("subj", seq_along(B_list))
    final_ids[missing] <- fallback[missing]
  }

  names(B_list) <- final_ids
  list(B_list = B_list, ids = final_ids)
}

#' Streaming prediction for new subjects via a loader
#'
#' @param object dkge | dkge_stream | dkge_model
#' @param loader object with n(), B(s) methods (and optional X(s))
#' @param contrasts list or matrix as in dkge_predict()
#' @return list(values=list per subject, A_list=list of loadings)
#' @export
dkge_predict_stream <- function(object, loader, contrasts) {
  comps <- .dkge_components(object)
  # alphas
  if (is.matrix(contrasts)) {
    Cmat <- contrasts; cn <- colnames(Cmat); if (is.null(cn)) cn <- paste0("c", seq_len(ncol(Cmat)))
    contrasts <- setNames(lapply(seq_len(ncol(Cmat)), function(j) Cmat[, j]), cn)
  }
  alpha_list <- lapply(contrasts, function(c) {
    ctil <- backsolve(comps$R, c, transpose = FALSE)
    t(comps$U) %*% comps$K %*% ctil
  })
  S <- loader$n()
  vals <- vector("list", S)
  A_list <- vector("list", S)
  for (s in seq_len(S)) {
    Bs <- .dkge_align_effects(.dkge_coerce_beta(loader$B(s)), comps$effects)
    Btil <- t(comps$R) %*% Bs
    A <- t(Btil) %*% comps$K %*% comps$U
    A_list[[s]] <- A
    res <- sapply(alpha_list, function(a) as.numeric(A %*% a))
    if (is.matrix(res) && length(alpha_list) > 1) {
      colnames(res) <- names(alpha_list)
    } else {
      names(res) <- names(alpha_list)
    }
    vals[[s]] <- res
  }
  names(vals) <- paste0("subj", seq_len(S))
  list(values = vals, A_list = A_list)
}

.dkge_resolve_predict_args <- function(newdata, args) {
  if (!missing(newdata) && !is.null(newdata)) {
    if (is.null(args$B_list) && !is.null(newdata$B_list)) {
      args$B_list <- newdata$B_list
    }
    if (is.null(args$B_list) && !is.null(newdata$betas)) {
      args$B_list <- newdata$betas
    }
    if (is.null(args$contrasts) && !is.null(newdata$contrasts)) {
      args$contrasts <- newdata$contrasts
    }
    if (is.null(args$return_loadings) && !is.null(newdata$return_loadings)) {
      args$return_loadings <- newdata$return_loadings
    }
  }
  args
}

.dkge_predict_dispatch <- function(object, newdata, args) {
  args <- .dkge_resolve_predict_args(newdata, args)
  if (is.null(args$B_list)) {
    stop("Provide `B_list` via `newdata` or the `B_list` argument.", call. = FALSE)
  }
  if (is.null(args$contrasts)) {
    stop("Provide `contrasts` via `newdata` or the `contrasts` argument.", call. = FALSE)
  }
  do.call(dkge_predict, c(list(object = object), args))
}

#' Predict contrasts for new subjects using a DKGE fit
#'
#' S3 front-end that forwards to [dkge_predict()] while accepting `newdata`
#' lists with `betas`/`B_list` and `contrasts` entries.
#'
#' @param object A `dkge` fit.
#' @param newdata Optional list with elements `betas`/`B_list` and `contrasts`.
#' @param ... Additional arguments passed to [dkge_predict()].
#' @return Output of [dkge_predict()].
#' @export
predict.dkge <- function(object, newdata = NULL, ...) {
  .dkge_predict_dispatch(object, newdata, list(...))
}

#' @rdname predict.dkge
#' @export
predict.dkge_model <- function(object, newdata = NULL, ...) {
  .dkge_predict_dispatch(object, newdata, list(...))
}
