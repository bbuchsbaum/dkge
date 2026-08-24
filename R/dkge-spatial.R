# dkge-spatial.R
# Sparse, model-level spatial regularization for DKGE.

#' Construct a model-level DKGE spatial regularizer
#'
#' Builds a sparse graph-Laplacian regularizer for the spatial columns of each
#' subject beta block. Coordinate-based construction delegates to
#' [adjoin::spatial_laplacian()]; DKGE validates the resulting operator, binds it
#' to the beta-column domain, and applies the resolvent
#' \deqn{H_\lambda = (I + \lambda L)^{-1}}
#' inside the fitted moment and every downstream component or contrast map.
#'
#' Supply either `coords` or `laplacian`. A single matrix defines one shared
#' spatial domain and therefore requires every subject to have the same number
#' and ordering of spatial units. A list allows subject-specific domains and may
#' be named with the fitted subject IDs. When calling [dkge_fit()] on raw lists,
#' use `dkge_data(..., subject_ids = ...)` to declare those IDs explicitly.
#' When both beta columns and `domain` are named, DKGE matches them by name and
#' fails if they do not describe the same units.
#'
#' @param coords Numeric coordinate matrix (`P x d`) or a list of such matrices,
#'   one per subject. Mutually exclusive with `laplacian`.
#' @param laplacian Precomputed graph Laplacian (`P x P`) or a list of
#'   Laplacians. It must be a symmetric combinatorial graph Laplacian: finite,
#'   non-negative diagonal, non-positive off-diagonal, and zero row sums.
#' @param lambda Non-negative Tikhonov smoothing strength. `0` is an exact
#'   identity operation and reproduces an unsmoothed DKGE fit.
#' @param dthresh,nnk,weight_mode,sigma,handle_isolates Arguments forwarded to
#'   [adjoin::spatial_laplacian()] when `coords` is supplied. DKGE does not allow
#'   `handle_isolates = "drop"`, because dropping rows would change the beta
#'   column domain.
#' @param normalized Must currently be `FALSE`. DKGE deliberately uses a
#'   combinatorial Laplacian whose null space contains the constant field;
#'   degree-normalized Laplacians generally do not preserve constants under the
#'   resolvent used here.
#' @param domain Optional spatial-unit labels. Supply one character vector for a
#'   shared matrix or a list aligned with subject-specific matrices. Coordinate
#'   row names or Laplacian dimnames are used when `domain` is omitted.
#'
#' @return An object of class `dkge_spatial_regularizer`, suitable for the
#'   `spatial` argument of [dkge()] or [dkge_fit()].
#' @export
#' @examples
#' coords <- cbind(x = 0:11, y = 0, z = 0)
#' spatial <- dkge_spatial_regularizer(
#'   coords = coords,
#'   lambda = 0.5,
#'   dthresh = 1.01,
#'   nnk = 3,
#'   weight_mode = "binary"
#' )
#' spatial
dkge_spatial_regularizer <- function(
    coords = NULL,
    laplacian = NULL,
    lambda = 1,
    dthresh = 1.42,
    nnk = 27L,
    weight_mode = c("binary", "heat"),
    sigma = dthresh / 2,
    normalized = FALSE,
    handle_isolates = c("keep_zero", "self_loop"),
    domain = NULL) {
  if (is.null(coords) == is.null(laplacian)) {
    .dkge_abort(
      "Supply exactly one of `coords` or `laplacian`.",
      "dkge_spatial_spec_error"
    )
  }
  lambda <- .dkge_spatial_lambda(lambda)
  weight_mode <- match.arg(weight_mode)
  handle_isolates <- match.arg(handle_isolates)
  if (!is.logical(normalized) || length(normalized) != 1L || is.na(normalized)) {
    .dkge_abort("`normalized` must be TRUE or FALSE.",
                "dkge_spatial_spec_error")
  }
  if (isTRUE(normalized)) {
    .dkge_abort(
      paste0(
        "`normalized` must be FALSE: DKGE spatial regularization requires a ",
        "combinatorial Laplacian that preserves constant fields."
      ),
      "dkge_spatial_spec_error"
    )
  }

  source <- if (is.null(coords)) "precomputed" else "adjoin"
  spatial_input <- coords %||% laplacian
  shared <- !(is.list(spatial_input) && !is.data.frame(spatial_input))
  input <- if (is.null(coords)) laplacian else coords
  entries <- if (shared) list(input) else input
  if (!is.list(entries) || !length(entries)) {
    .dkge_abort("Spatial inputs must contain at least one matrix.",
                "dkge_spatial_spec_error")
  }
  entry_names <- if (shared) NULL else names(entries)

  domains <- .dkge_spatial_domains(domain, entries, source = source,
                                   shared = shared)

  if (identical(source, "adjoin")) {
    if (!is.numeric(dthresh) || length(dthresh) != 1L ||
        !is.finite(dthresh) || dthresh <= 0) {
      .dkge_abort("`dthresh` must be one positive finite number.",
                  "dkge_spatial_spec_error")
    }
    if (!is.numeric(sigma) || length(sigma) != 1L ||
        !is.finite(sigma) || sigma <= 0) {
      .dkge_abort("`sigma` must be one positive finite number.",
                  "dkge_spatial_spec_error")
    }
    if (!is.numeric(nnk) || length(nnk) != 1L || !is.finite(nnk) ||
        nnk < 1L || nnk != as.integer(nnk)) {
      .dkge_abort("`nnk` must be one positive integer.",
                  "dkge_spatial_spec_error")
    }
    nnk <- as.integer(nnk)
    entries <- lapply(entries, .dkge_validate_spatial_coords)
    laplacians <- lapply(entries, function(x) {
      adjoin::spatial_laplacian(
        coord_mat = x,
        dthresh = dthresh,
        nnk = nnk,
        weight_mode = weight_mode,
        sigma = sigma,
        normalized = normalized,
        stochastic = FALSE,
        handle_isolates = handle_isolates
      )
    })
  } else {
    laplacians <- entries
  }

  laplacians <- Map(function(L, labels) {
    .dkge_validate_spatial_laplacian(L, domain = labels)
  }, laplacians, domains)
  if (!is.null(entry_names)) {
    names(laplacians) <- entry_names
    names(domains) <- entry_names
  }

  construction <- if (identical(source, "adjoin")) {
    list(
      function_name = "adjoin::spatial_laplacian",
      dthresh = as.numeric(dthresh),
      nnk = nnk,
      weight_mode = weight_mode,
      sigma = as.numeric(sigma),
      normalized = isTRUE(normalized),
      stochastic = FALSE,
      handle_isolates = handle_isolates,
      adjoin_version = as.character(utils::packageVersion("adjoin"))
    )
  } else {
    list(function_name = "precomputed")
  }

  structure(
    list(
      laplacians = laplacians,
      domains = domains,
      lambda = lambda,
      shared = shared,
      source = source,
      construction = construction
    ),
    class = c("dkge_spatial_regularizer", "list")
  )
}

#' @export
print.dkge_spatial_regularizer <- function(x, ...) {
  cat("<dkge_spatial_regularizer>\n")
  cat("  source   :", x$source, "\n")
  cat("  domains  :", if (isTRUE(x$shared)) "shared" else length(x$laplacians), "\n")
  cat("  lambda   :", format(x$lambda), "\n")
  sizes <- vapply(x$laplacians, nrow, integer(1))
  cat("  units    :", paste(sizes, collapse = ", "), "\n")
  invisible(x)
}

.dkge_spatial_lambda <- function(lambda) {
  if (!is.numeric(lambda) || length(lambda) != 1L ||
      !is.finite(lambda) || lambda < 0) {
    .dkge_abort("`lambda` must be one non-negative finite number.",
                "dkge_spatial_spec_error")
  }
  as.numeric(lambda)
}

.dkge_validate_spatial_coords <- function(x) {
  x <- as.matrix(x)
  if (!is.numeric(x) || nrow(x) < 1L || ncol(x) < 1L ||
      any(!is.finite(x))) {
    .dkge_abort(
      "Each coordinate input must be a non-empty finite numeric matrix.",
      "dkge_spatial_spec_error"
    )
  }
  storage.mode(x) <- "double"
  x
}

.dkge_spatial_domain_from_entry <- function(x, source) {
  labels <- if (identical(source, "adjoin")) {
    rownames(x)
  } else {
    rn <- rownames(x)
    cn <- colnames(x)
    if (!is.null(rn) && !is.null(cn) && !identical(rn, cn)) {
      .dkge_abort(
        "Laplacian row and column names must describe the same ordered domain.",
        "dkge_spatial_domain_error"
      )
    }
    rn %||% cn
  }
  if (is.null(labels)) NULL else as.character(labels)
}

.dkge_validate_spatial_domain <- function(domain, n, label = "domain") {
  if (is.null(domain)) return(NULL)
  domain <- as.character(domain)
  if (length(domain) != n || anyNA(domain) || any(!nzchar(domain)) ||
      anyDuplicated(domain)) {
    .dkge_abort(
      sprintf("`%s` must contain %d unique, non-empty labels.", label, n),
      "dkge_spatial_domain_error"
    )
  }
  domain
}

.dkge_spatial_domains <- function(domain, entries, source, shared) {
  n_entries <- length(entries)
  if (is.null(domain)) {
    return(lapply(entries, function(x) {
      .dkge_validate_spatial_domain(
        .dkge_spatial_domain_from_entry(x, source), nrow(x)
      )
    }))
  }
  domains <- if (is.list(domain)) domain else list(domain)
  if (shared && length(domains) != 1L) {
    .dkge_abort("A shared spatial matrix requires one `domain` vector.",
                "dkge_spatial_domain_error")
  }
  if (!shared && length(domains) != n_entries) {
    .dkge_abort(
      "Subject-specific spatial matrices require one `domain` vector per subject.",
      "dkge_spatial_domain_error"
    )
  }
  Map(function(labels, x, i) {
    .dkge_validate_spatial_domain(labels, nrow(x),
                                  sprintf("domain[[%d]]", i))
  }, domains, entries, seq_along(entries))
}

.dkge_validate_spatial_laplacian <- function(L, domain = NULL,
                                              tol = 1e-8) {
  if (!(is.matrix(L) || inherits(L, "Matrix")) ||
      length(dim(L)) != 2L || nrow(L) != ncol(L) || nrow(L) < 1L) {
    .dkge_abort("Each Laplacian must be a non-empty square matrix.",
                "dkge_spatial_laplacian_error")
  }
  L <- Matrix::Matrix(L, sparse = TRUE)
  if (any(!is.finite(L))) {
    .dkge_abort("Spatial Laplacians must contain only finite values.",
                "dkge_spatial_laplacian_error")
  }
  scale <- max(1, max(abs(L)))
  asym <- max(abs(L - Matrix::t(L)))
  if (!is.finite(asym) || asym > tol * scale) {
    .dkge_abort(
      sprintf("Spatial Laplacian must be symmetric (max asymmetry %.3e).", asym),
      "dkge_spatial_laplacian_error"
    )
  }
  L <- Matrix::forceSymmetric(0.5 * (L + Matrix::t(L)), uplo = "U")
  diagonal <- as.numeric(Matrix::diag(L))
  if (any(diagonal < -tol * scale)) {
    .dkge_abort("Spatial Laplacian diagonal entries must be non-negative.",
                "dkge_spatial_laplacian_error")
  }
  nz <- Matrix::summary(L)
  offdiag <- nz$i != nz$j
  if (any(nz$x[offdiag] > tol * scale)) {
    .dkge_abort("Spatial Laplacian off-diagonal entries must be non-positive.",
                "dkge_spatial_laplacian_error")
  }
  constant_error <- max(abs(Matrix::rowSums(L)))
  if (!is.finite(constant_error) || constant_error > tol * scale) {
    .dkge_abort(
      sprintf(
        "Spatial Laplacian must preserve constants (max absolute row sum %.3e).",
        constant_error
      ),
      "dkge_spatial_laplacian_error"
    )
  }
  labels <- .dkge_validate_spatial_domain(domain, nrow(L))
  if (!is.null(labels)) dimnames(L) <- list(labels, labels)
  L
}

.dkge_spatial_named_order <- function(items, subject_ids, what) {
  item_names <- names(items)
  if (is.null(item_names) || any(!nzchar(item_names))) return(items)
  if (anyDuplicated(item_names)) {
    .dkge_abort(sprintf("Named %s entries must have unique subject IDs.", what),
                "dkge_spatial_domain_error")
  }
  idx <- match(subject_ids, item_names)
  if (anyNA(idx)) {
    .dkge_abort(
      sprintf(
        "%s do not cover fit subject(s): %s.",
        what, paste(subject_ids[is.na(idx)], collapse = ", ")
      ),
      "dkge_spatial_domain_error"
    )
  }
  items[idx]
}

.dkge_spatial_operator <- function(L, lambda, domain, B, subject_id, source) {
  P <- ncol(B)
  if (nrow(L) != P) {
    .dkge_abort(
      sprintf(
        "Spatial domain for subject '%s' has %d units but its beta block has %d columns.",
        subject_id, nrow(L), P
      ),
      "dkge_spatial_domain_error"
    )
  }

  beta_domain <- colnames(B)
  domain_mode <- "positional"
  if (!is.null(domain)) {
    if (is.null(beta_domain)) {
      domain_mode <- "declared_positional"
    } else {
      beta_domain <- .dkge_validate_spatial_domain(
        beta_domain, P, sprintf("beta columns for subject '%s'", subject_id)
      )
      idx <- match(beta_domain, domain)
      if (anyNA(idx) || !setequal(beta_domain, domain)) {
        .dkge_abort(
          sprintf(
            "Spatial domain labels do not match beta columns for subject '%s'.",
            subject_id
          ),
          "dkge_spatial_domain_error"
        )
      }
      L <- L[idx, idx, drop = FALSE]
      domain <- beta_domain
      dimnames(L) <- list(domain, domain)
      domain_mode <- "label_matched"
    }
  }

  A <- NULL
  factor <- NULL
  if (lambda > 0) {
    A <- Matrix::Diagonal(P) + lambda * L
    A <- Matrix::forceSymmetric(A, uplo = "U")
    factor <- tryCatch(
      Matrix::Cholesky(A, LDL = FALSE, perm = TRUE),
      error = function(e) {
        .dkge_abort(
          sprintf(
            "Could not factor spatial regularizer for subject '%s': %s",
            subject_id, conditionMessage(e)
          ),
          "dkge_spatial_factor_error"
        )
      }
    )
  }

  nz <- Matrix::summary(L)
  # `L` is stored as a symmetric sparse matrix, whose summary contains one
  # triangular entry per undirected edge.
  edge_nnz <- sum(nz$i != nz$j)
  list(
    L = L,
    A = A,
    factor = factor,
    lambda = lambda,
    n_units = P,
    n_edges = as.integer(edge_nnz),
    domain = domain,
    domain_mode = domain_mode,
    subject_id = subject_id,
    source = source,
    fingerprint = digest::digest(
      list(dim = dim(L), domain = domain, entries = nz),
      algo = "xxhash64", serialize = TRUE
    )
  )
}

.dkge_resolve_spatial <- function(spatial, B_list, subject_ids = NULL) {
  if (is.null(spatial)) return(NULL)
  if (inherits(spatial, "dkge_spatial_fit")) spatial <- spatial$spec
  if (!inherits(spatial, "dkge_spatial_regularizer")) {
    .dkge_abort(
      "`spatial` must be created by `dkge_spatial_regularizer()`.",
      "dkge_spatial_spec_error"
    )
  }
  if (!is.list(B_list) || !length(B_list)) {
    .dkge_abort("Spatial regularization requires a non-empty beta-block list.",
                "dkge_spatial_domain_error")
  }
  S <- length(B_list)
  subject_ids <- as.character(subject_ids %||% names(B_list) %||%
                                paste0("subject", seq_len(S)))
  if (length(subject_ids) != S || anyNA(subject_ids) || any(!nzchar(subject_ids))) {
    subject_ids <- paste0("subject", seq_len(S))
  }

  if (isTRUE(spatial$shared)) {
    laplacians <- rep(spatial$laplacians, S)
    domains <- rep(spatial$domains, S)
  } else {
    if (length(spatial$laplacians) != S) {
      .dkge_abort(
        sprintf(
          "Subject-specific `spatial` has %d domains but the fit has %d subjects.",
          length(spatial$laplacians), S
        ),
        "dkge_spatial_domain_error"
      )
    }
    laplacians <- .dkge_spatial_named_order(
      spatial$laplacians, subject_ids, "spatial Laplacians"
    )
    domains <- .dkge_spatial_named_order(
      spatial$domains, subject_ids, "spatial domains"
    )
  }

  operators <- Map(function(L, domain, B, id) {
    .dkge_spatial_operator(L, spatial$lambda, domain, B, id,
                           source = spatial$source)
  }, laplacians, domains, B_list, subject_ids)
  names(operators) <- subject_ids
  if (isTRUE(spatial$shared) && length(operators) > 1L) {
    fingerprints <- vapply(operators, `[[`, character(1), "fingerprint")
    if (length(unique(fingerprints)) > 1L) {
      .dkge_abort(
        paste0(
          "A shared spatial regularizer requires one common beta-column ",
          "ordering across subjects. Reorder the beta columns or supply a ",
          "subject-specific list of spatial domains."
        ),
        "dkge_spatial_domain_error"
      )
    }
  }

  diagnostics <- do.call(rbind, lapply(operators, function(op) {
    data.frame(
      subject = op$subject_id,
      n_units = op$n_units,
      n_edges = op$n_edges,
      lambda = op$lambda,
      source = op$source,
      domain_mode = op$domain_mode,
      fingerprint = op$fingerprint,
      stringsAsFactors = FALSE
    )
  }))
  rownames(diagnostics) <- NULL

  provenance <- list(
    method = "laplacian_resolvent",
    formula = "H = solve(I + lambda * L)",
    lambda = spatial$lambda,
    shared = isTRUE(spatial$shared),
    source = spatial$source,
    construction = spatial$construction,
    subjects = subject_ids,
    domains = lapply(operators, function(op) {
      list(
        n_units = op$n_units,
        domain_mode = op$domain_mode,
        fingerprint = op$fingerprint
      )
    })
  )

  structure(
    list(
      active = spatial$lambda > 0,
      lambda = spatial$lambda,
      shared = isTRUE(spatial$shared),
      operators = operators,
      diagnostics = diagnostics,
      provenance = provenance,
      spec = spatial
    ),
    class = c("dkge_spatial_fit", "list")
  )
}

.dkge_spatial_apply_betas <- function(B, operator = NULL) {
  B <- as.matrix(B)
  if (is.null(operator) || operator$lambda == 0) return(B)
  if (ncol(B) != operator$n_units) {
    .dkge_abort(
      sprintf(
        "Spatial operator has %d units but the beta block has %d columns.",
        operator$n_units, ncol(B)
      ),
      "dkge_spatial_domain_error"
    )
  }
  beta_domain <- colnames(B)
  if (!is.null(operator$domain) && !is.null(beta_domain) &&
      !identical(as.character(beta_domain), operator$domain)) {
    .dkge_abort(
      sprintf(
        "Beta-column order does not match the bound spatial domain for subject '%s'.",
        operator$subject_id
      ),
      "dkge_spatial_domain_error"
    )
  }
  beta_dimnames <- dimnames(B)
  smoothed <- Matrix::solve(operator$factor, t(B))
  out <- t(as.matrix(smoothed))
  dimnames(out) <- beta_dimnames
  out
}

.dkge_fit_spatial_operator <- function(fit, subject = NULL, n_cols = NULL) {
  spatial <- fit$spatial
  if (is.null(spatial) || !isTRUE(spatial$active)) return(NULL)
  if (is.null(subject)) {
    if (!isTRUE(spatial$shared)) {
      .dkge_abort(
        "This fit has subject-specific spatial domains; supply `subject`.",
        "dkge_spatial_domain_error"
      )
    }
    subject <- 1L
  } else {
    subject <- .dkge_resolve_subject(fit, subject)
  }
  op <- spatial$operators[[subject]]
  if (!is.null(n_cols) && op$n_units != n_cols) {
    .dkge_abort(
      sprintf("Spatial operator has %d units but the block has %d columns.",
              op$n_units, n_cols),
      "dkge_spatial_domain_error"
    )
  }
  op
}

.dkge_apply_fit_spatial <- function(fit, B, subject = NULL) {
  op <- .dkge_fit_spatial_operator(fit, subject = subject, n_cols = ncol(B))
  .dkge_spatial_apply_betas(B, op)
}

.dkge_fit_model_btil <- function(fit, subject, Btil = NULL,
                                  voxel_weights = NULL) {
  s <- .dkge_resolve_subject(fit, subject)
  Btil <- as.matrix(Btil %||% fit$Btil[[s]])
  Bwork <- .dkge_scale_effect_columns(Btil, voxel_weights)
  .dkge_apply_fit_spatial(fit, Bwork, subject = s)
}

.dkge_fit_subject_loadings <- function(fit, U = fit$U,
                                       voxel_weights = NULL) {
  out <- lapply(seq_along(fit$Btil), function(s) {
    w_s <- if (is.list(voxel_weights)) voxel_weights[[s]] else voxel_weights
    Bmodel <- .dkge_fit_model_btil(fit, s, voxel_weights = w_s)
    t(Bmodel) %*% fit$K %*% U
  })
  names(out) <- names(fit$Btil) %||% fit$subject_ids
  out
}

.dkge_spatial_with_lambda <- function(spatial, lambda) {
  if (inherits(spatial, "dkge_spatial_fit")) spatial <- spatial$spec
  if (!inherits(spatial, "dkge_spatial_regularizer")) {
    .dkge_abort("`spatial` must be a DKGE spatial regularizer.",
                "dkge_spatial_spec_error")
  }
  spatial$lambda <- .dkge_spatial_lambda(lambda)
  spatial
}

.dkge_prediction_spatial <- function(object, B_list, spatial = NULL) {
  fitted_spatial <- object$spatial %||% NULL
  active_fit <- !is.null(fitted_spatial) && isTRUE(fitted_spatial$active)
  if (!active_fit) {
    if (!is.null(spatial)) {
      .dkge_abort(
        "Prediction-time smoothing cannot be added to an unsmoothed fit; refit with `spatial`.",
        "dkge_spatial_prediction_error"
      )
    }
    return(NULL)
  }
  if (!is.null(spatial)) {
    return(.dkge_resolve_spatial(spatial, B_list,
                                 subject_ids = names(B_list)))
  }
  if (isTRUE(fitted_spatial$shared)) {
    return(.dkge_resolve_spatial(fitted_spatial$spec, B_list,
                                 subject_ids = names(B_list)))
  }
  .dkge_abort(
    paste0(
      "This model was fit with subject-specific spatial domains. Supply a ",
      "`spatial` regularizer for the prediction subjects."
    ),
    "dkge_spatial_prediction_error"
  )
}
