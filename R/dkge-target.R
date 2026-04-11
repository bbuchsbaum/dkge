# dkge-target.R
# Subject-level target objects for between-subject DKGE models.

#' Build a subject-by-feature target for between-subject DKGE models
#'
#' Creates the explicit boundary object between the DKGE representation stage
#' and between-subject multivariate models. The returned target stores a
#' subject-by-feature matrix plus optional coverage, precision, and weighting
#' metadata.
#'
#' @param fit Optional `dkge` object used to derive DKGE targets.
#' @param type Target type. `"matrix"` wraps a supplied matrix, `"component_scores"`
#'   averages subject component expressions, and `"transported_maps"` builds
#'   a common-space matrix from contrast values.
#' @param Y Optional subject-by-feature matrix.
#' @param contrast Optional DKGE contrast used when `type = "transported_maps"`.
#' @param contrast_obj Optional precomputed `dkge_contrasts` object.
#' @param transport Optional `dkge_transport_spec` used for transported maps.
#' @param values Optional subject values. May be a matrix or a list of vectors.
#' @param loadings Optional subject loading matrices used by transport.
#' @param centroids,sizes,medoid,mapper Transport inputs overriding `transport`.
#' @param crossfit Cross-fitting method used when computing `contrast`.
#' @param feature_ids,subject_ids Optional identifiers.
#' @param coverage Optional coverage matrix/vector aligned with `Y`.
#' @param precision Optional precision matrix/vector aligned with `Y`.
#' @param subject_weights Optional subject reliability weights.
#' @param feature_weights Optional feature reliability weights.
#' @param geometry Optional feature geometry metadata.
#' @param provenance Optional provenance metadata.
#' @param ... Additional arguments forwarded to contrast/transport helpers.
#'
#' @return Object of class `dkge_target`.
#' @export
dkge_make_target <- function(fit = NULL,
                             type = c("matrix", "component_scores", "transported_maps"),
                             Y = NULL,
                             contrast = NULL,
                             contrast_obj = NULL,
                             transport = NULL,
                             values = NULL,
                             loadings = NULL,
                             centroids = NULL,
                             sizes = NULL,
                             medoid = NULL,
                             mapper = NULL,
                             crossfit = c("none", "analytic", "loso", "kfold"),
                             feature_ids = NULL,
                             subject_ids = NULL,
                             coverage = NULL,
                             precision = NULL,
                             subject_weights = NULL,
                             feature_weights = NULL,
                             geometry = NULL,
                             provenance = NULL,
                             ...) {
  type <- match.arg(type)
  crossfit <- match.arg(crossfit)

  if (!is.null(Y)) {
    Y <- as.matrix(Y)
    type <- "matrix"
  } else if (type == "matrix") {
    if (is.null(values)) {
      stop("`Y` or `values` must be supplied for type = 'matrix'.", call. = FALSE)
    }
    Y <- .dkge_values_to_matrix(values, subject_ids = subject_ids)
  } else if (type == "component_scores") {
    if (is.null(fit) && is.null(loadings)) {
      stop("`fit` or `loadings` is required for component score targets.", call. = FALSE)
    }
    if (is.null(loadings)) {
      loadings <- .dkge_subject_loadings(fit)
    }
    Y <- .dkge_component_score_matrix(loadings, sizes = sizes)
    if (is.null(feature_ids)) {
      feature_ids <- colnames(Y) %||% paste0("component", seq_len(ncol(Y)))
    }
  } else if (type == "transported_maps") {
    payload <- .dkge_target_transport_payload(fit = fit,
                                              contrast = contrast,
                                              contrast_obj = contrast_obj,
                                              transport = transport,
                                              values = values,
                                              loadings = loadings,
                                              centroids = centroids,
                                              sizes = sizes,
                                              medoid = medoid,
                                              mapper = mapper,
                                              crossfit = crossfit,
                                              ...)
    Y <- payload$Y
    geometry <- geometry %||% payload$geometry
    provenance <- provenance %||% payload$provenance
    if (is.null(feature_ids)) feature_ids <- payload$feature_ids
    if (is.null(subject_ids)) subject_ids <- payload$subject_ids
  }

  if (!is.matrix(Y) || !is.numeric(Y)) {
    stop("Target `Y` must be a numeric matrix.", call. = FALSE)
  }
  if (any(!is.finite(Y))) {
    stop("Target `Y` contains non-finite values.", call. = FALSE)
  }

  subject_ids <- subject_ids %||% rownames(Y)
  if (is.null(subject_ids)) {
    if (!is.null(fit) && !is.null(fit$subject_ids) && length(fit$subject_ids) == nrow(Y)) {
      subject_ids <- fit$subject_ids
    } else {
      subject_ids <- paste0("subj", seq_len(nrow(Y)))
    }
  }
  subject_ids <- as.character(subject_ids)
  if (length(subject_ids) != nrow(Y) || any(!nzchar(subject_ids)) || any(duplicated(subject_ids))) {
    stop("`subject_ids` must be unique, non-empty, and match nrow(Y).", call. = FALSE)
  }
  rownames(Y) <- subject_ids

  feature_ids <- feature_ids %||% colnames(Y)
  if (is.null(feature_ids)) {
    feature_ids <- paste0("feature", seq_len(ncol(Y)))
  }
  feature_ids <- as.character(feature_ids)
  if (length(feature_ids) != ncol(Y) || any(!nzchar(feature_ids)) || any(duplicated(feature_ids))) {
    stop("`feature_ids` must be unique, non-empty, and match ncol(Y).", call. = FALSE)
  }
  colnames(Y) <- feature_ids

  coverage <- .dkge_validate_target_sidecar(coverage, nrow(Y), ncol(Y), "coverage")
  precision <- .dkge_validate_target_sidecar(precision, nrow(Y), ncol(Y), "precision")
  subject_weights <- .dkge_validate_target_weights(subject_weights, nrow(Y), "subject_weights")
  feature_weights <- .dkge_validate_target_weights(feature_weights, ncol(Y), "feature_weights")
  mask <- rep(TRUE, ncol(Y))
  if (!is.null(feature_weights)) {
    mask <- mask & feature_weights > 0
  }
  if (!is.null(coverage) && is.vector(coverage)) {
    mask <- mask & as.numeric(coverage) > 0
  }

  out <- list(
    Y = Y,
    subject_ids = subject_ids,
    feature_ids = feature_ids,
    type = type,
    feature_space = switch(type,
                           matrix = "features",
                           component_scores = "components",
                           transported_maps = "transported",
                           "features"),
    geometry = geometry,
    contrast = contrast,
    transport = transport,
    crossfit = crossfit,
    coverage = coverage,
    precision = precision,
    mask = mask,
    subject_weights = subject_weights,
    feature_weights = feature_weights,
    provenance = provenance
  )
  class(out) <- c("dkge_target", "list")
  out
}

.dkge_values_to_matrix <- function(values, subject_ids = NULL) {
  if (is.matrix(values) || is.data.frame(values)) {
    return(as.matrix(values))
  }
  if (!is.list(values) || !length(values)) {
    stop("`values` must be a matrix/data.frame or non-empty list of numeric vectors.", call. = FALSE)
  }
  lens <- vapply(values, length, integer(1))
  if (length(unique(lens)) != 1L) {
    stop("List `values` must have equal vector lengths to stack without transport.", call. = FALSE)
  }
  M <- do.call(rbind, lapply(values, as.numeric))
  rownames(M) <- names(values) %||% subject_ids
  M
}

.dkge_subject_loadings <- function(fit) {
  stopifnot(inherits(fit, "dkge"))
  if (!is.null(fit$Btil)) {
    out <- lapply(fit$Btil, function(Bts) t(Bts) %*% fit$K %*% fit$U)
    names(out) <- fit$subject_ids
    return(out)
  }
  stop("`fit` does not contain Btil; supply `loadings` explicitly.", call. = FALSE)
}

.dkge_component_score_matrix <- function(loadings, sizes = NULL) {
  if (!is.list(loadings) || !length(loadings)) {
    stop("`loadings` must be a non-empty list of matrices.", call. = FALSE)
  }
  if (is.null(sizes)) {
    sizes <- vector("list", length(loadings))
  }
  stopifnot(length(sizes) == length(loadings))
  Y <- Map(function(A, w) {
    A <- as.matrix(A)
    if (is.null(w)) {
      colMeans(A)
    } else {
      w <- as.numeric(w)
      if (length(w) != nrow(A)) stop("Each `sizes` vector must match loading rows.", call. = FALSE)
      w <- pmax(w, 0)
      if (sum(w) <= 0) w <- rep(1, length(w))
      colSums(A * w) / sum(w)
    }
  }, loadings, sizes)
  M <- do.call(rbind, Y)
  rownames(M) <- names(loadings)
  colnames(M) <- colnames(loadings[[1]]) %||% paste0("component", seq_len(ncol(M)))
  M
}

.dkge_target_transport_payload <- function(fit,
                                           contrast,
                                           contrast_obj,
                                           transport,
                                           values,
                                           loadings,
                                           centroids,
                                           sizes,
                                           medoid,
                                           mapper,
                                           crossfit,
                                           ...) {
  if (!is.null(values) && is.matrix(values)) {
    Y <- as.matrix(values)
    return(list(Y = Y, feature_ids = colnames(Y), subject_ids = rownames(Y),
                geometry = NULL, provenance = list(source = "values")))
  }

  if (is.null(fit)) {
    stop("`fit` is required to build transported map targets unless `values` is a matrix.", call. = FALSE)
  }
  stopifnot(inherits(fit, "dkge"))

  if (is.null(loadings)) {
    loadings <- .dkge_subject_loadings(fit)
  }

  if (!is.null(transport) && inherits(transport, "dkge_transport_spec")) {
    centroids <- centroids %||% transport$centroids
    sizes <- sizes %||% transport$sizes
    medoid <- medoid %||% transport$medoid
    mapper <- mapper %||% transport$mapper
    dots <- c(list(...),
              list(epsilon = transport$epsilon,
                   max_iter = transport$max_iter,
                   tol = transport$tol,
                   lambda_emb = transport$lambda_emb,
                   lambda_spa = transport$lambda_spa,
                   sigma_mm = transport$sigma_mm,
                   lambda_size = transport$lambda_size))
  } else {
    dots <- list(...)
  }
  if (is.null(centroids)) {
    centroids <- fit$centroids %||% fit$input$centroids
  }
  if (is.null(medoid)) medoid <- 1L

  if (!is.null(values)) {
    if (is.null(centroids)) {
      Y <- .dkge_values_to_matrix(values)
      return(list(Y = Y, feature_ids = colnames(Y), subject_ids = rownames(Y),
                  geometry = NULL, provenance = list(source = "values")))
    }
    mapper_spec <- .dkge_resolve_mapper_spec(mapper, method = NULL, dots = dots)
    if (is.null(sizes)) sizes <- lapply(loadings, function(A) rep(1, nrow(A)))
    tr <- .dkge_transport_to_medoid(mapper_spec, values, loadings, centroids,
                                    sizes = sizes, medoid = medoid)
    Y <- tr$subj_values
    return(list(Y = Y,
                feature_ids = colnames(Y) %||% paste0("medoid", seq_len(ncol(Y))),
                subject_ids = rownames(Y) %||% names(values),
                geometry = list(centroids = centroids[[medoid]], medoid = medoid),
                provenance = list(source = "values", transport = tr)))
  }

  if (is.null(contrast_obj)) {
    if (is.null(contrast)) {
      stop("`contrast`, `contrast_obj`, or `values` is required for transported map targets.", call. = FALSE)
    }
    method <- if (identical(crossfit, "none")) "analytic" else crossfit
    contrast_obj <- dkge_contrast(fit, contrast, method = method, align = FALSE, ...)
  }

  if (is.null(centroids)) {
    first <- contrast_obj$values[[1]]
    Y <- .dkge_values_to_matrix(first)
    return(list(Y = Y,
                feature_ids = colnames(Y),
                subject_ids = rownames(Y),
                geometry = NULL,
                provenance = list(source = "contrast", contrast_obj = contrast_obj)))
  }

  tr <- do.call(dkge_transport_contrasts_to_medoid,
                c(list(fit = fit,
                       contrast_obj = contrast_obj,
                       medoid = medoid,
                       centroids = centroids,
                       loadings = loadings,
                       sizes = sizes,
                       mapper = mapper),
                  dots))
  mats <- lapply(names(tr), function(nm) {
    M <- tr[[nm]]$subj_values
    colnames(M) <- paste(nm, colnames(M) %||% seq_len(ncol(M)), sep = ":")
    M
  })
  Y <- do.call(cbind, mats)
  list(Y = Y,
       feature_ids = colnames(Y),
       subject_ids = rownames(Y),
       geometry = list(centroids = centroids[[medoid]], medoid = medoid),
       provenance = list(source = "contrast", contrast_obj = contrast_obj,
                         transport = tr))
}

.dkge_validate_target_sidecar <- function(x, n, p, nm) {
  if (is.null(x)) return(NULL)
  if (is.matrix(x) || is.data.frame(x)) {
    x <- as.matrix(x)
    if (!is.numeric(x) || !identical(dim(x), c(n, p))) {
      stop("`", nm, "` matrix must be numeric with the same dimensions as Y.", call. = FALSE)
    }
    return(x)
  }
  x <- as.numeric(x)
  if (length(x) != p) {
    stop("`", nm, "` vector must have length ncol(Y).", call. = FALSE)
  }
  x
}

.dkge_validate_target_weights <- function(x, len, nm) {
  if (is.null(x)) return(NULL)
  x <- as.numeric(x)
  if (length(x) != len || any(!is.finite(x)) || any(x < 0)) {
    stop("`", nm, "` must be a non-negative finite vector with expected length.", call. = FALSE)
  }
  x
}

#' @export
print.dkge_target <- function(x, ...) {
  cat("<dkge_target>", "\n", sep = "")
  cat("  type     :", x$type, "\n")
  cat("  subjects :", nrow(x$Y), "\n")
  cat("  features :", ncol(x$Y), "\n")
  if (!is.null(x$coverage)) cat("  coverage : yes\n")
  if (!is.null(x$precision)) cat("  precision: yes\n")
  invisible(x)
}
