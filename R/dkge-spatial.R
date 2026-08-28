# dkge-spatial.R
# Sparse, model-level spatial regularization for DKGE.

#' Construct a model-level DKGE spatial regularizer
#'
#' Builds a sparse graph-Laplacian regularizer for the spatial columns of each
#' subject beta block. Coordinate-based construction delegates to
#' [adjoin::spatial_laplacian()]; DKGE validates the resulting operator, binds it
#' to the beta-column domain, and smooths the beta block with the Tikhonov
#' resolvent
#' \deqn{H_\lambda = (I + \lambda L)^{-1},}
#' replacing \eqn{\widetilde B_s} by \eqn{\widetilde B_s H_\lambda} in the
#' fitted moment and in every downstream component or contrast map.
#'
#' @section Effect on the pooled moment:
#' Because the resolvent is applied to the beta block and the moment is then
#' formed from the smoothed block, \eqn{H_\lambda} enters the \eqn{q \times q}
#' moment **twice**:
#' \deqn{M_s(\lambda) = \widetilde B_s H_\lambda \Omega_s H_\lambda \widetilde B_s^{\mathsf T},}
#' which reduces to \eqn{\widetilde B_s H_\lambda^2 \widetilde B_s^{\mathsf T}}
#' under an identity spatial metric. The effective moment-level smoothing is
#' therefore the *squared* resolvent, not \eqn{H_\lambda}; keep that in mind
#' when comparing `lambda` against a target smoothing kernel.
#'
#' Two per-unit weightings sit on **opposite sides** of the smoother. Writing
#' \eqn{W} for [dkge_weights()] voxel weights and \eqn{\Omega_s} for the spatial
#' metric in `Omega_list`, the moment DKGE actually accumulates is
#' \deqn{M_s(\lambda) = (\widetilde B_s W^{1/2}) H_\lambda \Omega_s H_\lambda (\widetilde B_s W^{1/2})^{\mathsf T}.}
#' Voxel weights are applied to the beta columns *before* smoothing, so a
#' down-weighted unit is shrunk first and then diffuses its reduced value to its
#' neighbours -- a zero-weight unit acts as a hole in the field rather than
#' being removed from the graph. `Omega_list` is applied *after* smoothing, as a
#' metric on the smoothed field. To exclude a unit completely, remove the
#' matching beta column, graph node, domain label, and every aligned per-unit
#' input; changing only `coords`, `laplacian`, or a voxel weight is not a domain
#' exclusion operation.
#'
#' @section Choosing `dthresh`:
#' `dthresh` is expressed in the units of `coords`. The defaults
#' (`dthresh = 1.42`, `nnk = 27`) assume **unit-spaced voxel indices**, where
#' 1.42 admits face and edge neighbours of a 3x3x3 neighbourhood. Passing
#' millimetre coordinates with the default threshold isolates every unit, which
#' makes `L` the zero matrix and the resolvent the identity: the fit is then
#' bit-identical to an unregularized one despite a positive `lambda`.
#' Construction warns when a positive-`lambda` graph ends up with no edges.
#' Spatial status is `inactive` for `lambda = 0`, `inert` when all graphs are
#' edgeless, `partial` when only some subject graphs are edgeless, and `active`
#' when every graph can smooth. Both `print()` and
#' `dkge_diagnostics(fit)$spatial$diagnostics` report the status and edge count.
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
#'   [adjoin::spatial_laplacian()] when `coords` is supplied. `dthresh` is a
#'   distance cutoff in the units of `coords` and `nnk` caps the number of
#'   neighbours; see the "Choosing `dthresh`" section, as the defaults assume
#'   unit-spaced voxel indices rather than millimetres. DKGE does not allow
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
  topology <- .dkge_spatial_topology(laplacians, lambda)
  .dkge_spatial_warn_topology(topology, source, dthresh)
  topology_warning <- if (topology$status %in% c("inert", "partial")) {
    .dkge_spatial_warning_receipt(topology)
  } else {
    NULL
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

  out <- structure(
    list(
      laplacians = laplacians,
      domains = domains,
      lambda = lambda,
      shared = shared,
      source = source,
      construction = construction,
      topology = topology,
      topology_warning = topology_warning
    ),
    class = c("dkge_spatial_regularizer", "list")
  )
  .dkge_spatial_spec_payload(out, validate = TRUE)
  out
}

#' Canonical payload for a public spatial regularizer specification
#'
#' A `dkge_spatial_regularizer` is an ordinary mutable R list. Revalidate its
#' complete structural contract whenever it crosses a fit, CV, retune, freeze,
#' or prediction boundary so post-construction mutations cannot be laundered
#' into a fitted estimator.
#'
#' @keywords internal
#' @noRd
.dkge_spatial_spec_payload <- function(spatial, validate = TRUE) {
  required_names <- c(
    "laplacians", "domains", "lambda", "shared", "source",
    "construction", "topology", "topology_warning"
  )
  if (!inherits(spatial, "dkge_spatial_regularizer") ||
      !identical(names(unclass(spatial)), required_names)) {
    .dkge_abort(
      "The spatial regularizer has a foreign or malformed field schema.",
      "dkge_spatial_provenance_error"
    )
  }
  lambda <- .dkge_spatial_lambda(spatial$lambda)
  if (!is.logical(spatial$shared) || length(spatial$shared) != 1L ||
      is.na(spatial$shared) ||
      !is.character(spatial$source) || length(spatial$source) != 1L ||
      is.na(spatial$source) ||
      !spatial$source %in% c("precomputed", "adjoin") ||
      !is.list(spatial$laplacians) || !length(spatial$laplacians) ||
      !is.list(spatial$domains) ||
      length(spatial$domains) != length(spatial$laplacians) ||
      (isTRUE(spatial$shared) && length(spatial$laplacians) != 1L)) {
    .dkge_abort("The spatial regularizer structure is invalid.",
                "dkge_spatial_provenance_error")
  }
  laplacian_names <- names(spatial$laplacians)
  domain_names <- names(spatial$domains)
  names_valid <- if (isTRUE(spatial$shared)) {
    is.null(laplacian_names) && is.null(domain_names)
  } else if (is.null(laplacian_names) && is.null(domain_names)) {
    TRUE
  } else {
    !is.null(laplacian_names) && !is.null(domain_names) &&
      identical(laplacian_names, domain_names) &&
      !anyNA(laplacian_names) && all(nzchar(laplacian_names)) &&
      !anyDuplicated(laplacian_names)
  }
  if (!isTRUE(names_valid)) {
    .dkge_abort(
      "Spatial Laplacian and domain subject-name schemas do not match.",
      "dkge_spatial_provenance_error"
    )
  }

  domains <- Map(function(domain, L) {
    .dkge_validate_spatial_domain(domain, nrow(L), "spatial domain")
  }, spatial$domains, spatial$laplacians)
  domain_dimnames_valid <- Map(function(domain, L) {
    if (is.null(domain)) {
      is.null(rownames(L)) && is.null(colnames(L))
    } else {
      identical(rownames(L), domain) && identical(colnames(L), domain)
    }
  }, domains, spatial$laplacians)
  if (!all(unlist(domain_dimnames_valid, use.names = FALSE))) {
    .dkge_abort(
      "Each spatial domain must exactly equal its Laplacian dimnames.",
      "dkge_spatial_domain_error"
    )
  }
  laplacians <- Map(function(L, domain) {
    .dkge_validate_spatial_laplacian(L, domain = domain)
  }, spatial$laplacians, domains)
  names(domains) <- domain_names
  names(laplacians) <- laplacian_names
  if (isTRUE(validate) &&
      (!identical(spatial$domains, domains) ||
       !identical(spatial$laplacians, laplacians))) {
    .dkge_abort(
      "The spatial regularizer is not in its canonical validated form.",
      "dkge_spatial_provenance_error"
    )
  }

  construction <- spatial$construction
  construction_valid <- if (identical(spatial$source, "precomputed")) {
    identical(construction, list(function_name = "precomputed"))
  } else {
    expected_names <- c(
      "function_name", "dthresh", "nnk", "weight_mode", "sigma",
      "normalized", "stochastic", "handle_isolates", "adjoin_version"
    )
    is.list(construction) &&
      identical(names(construction), expected_names) &&
      identical(construction$function_name, "adjoin::spatial_laplacian") &&
      is.numeric(construction$dthresh) && length(construction$dthresh) == 1L &&
      is.finite(construction$dthresh) && construction$dthresh > 0 &&
      is.numeric(construction$nnk) && length(construction$nnk) == 1L &&
      is.finite(construction$nnk) && construction$nnk >= 1L &&
      construction$nnk == as.integer(construction$nnk) &&
      construction$weight_mode %in% c("binary", "heat") &&
      is.numeric(construction$sigma) && length(construction$sigma) == 1L &&
      is.finite(construction$sigma) && construction$sigma > 0 &&
      identical(construction$normalized, FALSE) &&
      identical(construction$stochastic, FALSE) &&
      construction$handle_isolates %in% c("keep_zero", "self_loop") &&
      is.character(construction$adjoin_version) &&
      length(construction$adjoin_version) == 1L &&
      !is.na(construction$adjoin_version) &&
      nzchar(construction$adjoin_version)
  }
  if (!isTRUE(construction_valid)) {
    .dkge_abort(
      "The spatial regularizer source or construction receipt is invalid.",
      "dkge_spatial_provenance_error"
    )
  }

  topology <- .dkge_spatial_topology(laplacians, lambda)
  warning_expected <- if (topology$status %in% c("inert", "partial")) {
    .dkge_spatial_warning_receipt(topology)
  } else {
    NULL
  }
  warning_valid <- if (is.null(warning_expected)) {
    is.null(spatial$topology_warning)
  } else {
    is.null(spatial$topology_warning) ||
      identical(spatial$topology_warning, warning_expected)
  }
  if (isTRUE(validate) &&
      (!identical(spatial$topology, topology) || !isTRUE(warning_valid))) {
    .dkge_abort(
      "The spatial topology or warning receipt is stale or inconsistent.",
      "dkge_spatial_provenance_error"
    )
  }
  list(
    laplacians = laplacians,
    domains = domains,
    lambda = lambda,
    shared = spatial$shared,
    source = spatial$source,
    construction = construction,
    topology = topology,
    topology_warning = spatial$topology_warning
  )
}

#' Count undirected edges in a validated Laplacian
#'
#' `Matrix::summary()` on a symmetric sparse Laplacian lists one triangular
#' entry per undirected edge, so off-diagonal entries count edges directly.
#'
#' @keywords internal
#' @noRd
.dkge_spatial_edge_count <- function(L) {
  nz <- Matrix::summary(L)
  sum(nz$i != nz$j)
}

#' Summarize whether a spatial specification can have a numerical effect
#'
#' This calculation is sparse and depends only on graph sizes and off-diagonal
#' nonzeros. It deliberately does not guess whether coordinates are expressed
#' in indices or physical units.
#'
#' @keywords internal
#' @noRd
.dkge_spatial_topology <- function(laplacians, lambda) {
  lambda <- .dkge_spatial_lambda(lambda)
  if (!is.list(laplacians) || !length(laplacians)) {
    .dkge_abort("Spatial topology requires at least one Laplacian.",
                "dkge_spatial_spec_error")
  }
  n_units <- vapply(laplacians, nrow, integer(1))
  n_edges <- vapply(laplacians, .dkge_spatial_edge_count, numeric(1))
  labels <- names(laplacians)
  if (is.null(labels)) labels <- rep("", length(laplacians))
  missing_labels <- !nzchar(labels)
  labels[missing_labels] <- as.character(which(missing_labels))
  empty <- which(n_edges == 0)
  requested <- lambda > 0
  effective_by_graph <- requested & n_edges > 0
  status <- if (!requested) {
    "inactive"
  } else if (!any(effective_by_graph)) {
    "inert"
  } else if (!all(effective_by_graph)) {
    "partial"
  } else {
    "active"
  }
  fingerprint <- digest::digest(
    list(n_units = n_units, n_edges = n_edges, labels = labels),
    algo = "xxhash64", serialize = TRUE
  )
  list(
    status = status,
    lambda = lambda,
    requested = requested,
    effective = any(effective_by_graph),
    fully_effective = requested && all(effective_by_graph),
    effective_by_graph = effective_by_graph,
    n_graphs = length(laplacians),
    n_units = n_units,
    n_edges = n_edges,
    empty = empty,
    empty_labels = labels[empty],
    fingerprint = fingerprint
  )
}

.dkge_spatial_warning_receipt <- function(topology) {
  list(
    fingerprint = topology$fingerprint,
    lambda_requested = isTRUE(topology$requested)
  )
}

.dkge_spatial_warning_is_current <- function(spatial, topology) {
  receipt <- spatial$topology_warning
  is.list(receipt) &&
    identical(receipt$fingerprint, topology$fingerprint) &&
    isTRUE(receipt$lambda_requested) &&
    isTRUE(topology$requested)
}

.dkge_spatial_topology_hint <- function(source, dthresh = NULL) {
  hint <- if (identical(source, "adjoin")) {
    threshold <- if (is.numeric(dthresh) && length(dthresh) == 1L &&
                     is.finite(dthresh)) format(dthresh) else "unknown"
    sprintf(
      paste0(" `dthresh` is measured in the units of `coords`; the default ",
             "assumes unit-spaced indices. With this graph dthresh = %s, so ",
             "millimetre coordinates would need a threshold on the millimetre ",
             "scale (e.g. slightly above the voxel spacing)."),
      threshold
    )
  } else {
    " The supplied Laplacian has no off-diagonal entries."
  }
  hint
}

#' Warn when one or more spatial graphs cannot smooth anything
#'
#' @keywords internal
#' @noRd
.dkge_spatial_warn_topology <- function(topology, source, dthresh = NULL) {
  if (!topology$status %in% c("inert", "partial")) {
    return(invisible(NULL))
  }
  if (identical(topology$status, "inert")) {
    who <- if (topology$n_graphs == 1L) {
      "The spatial graph has"
    } else {
      sprintf("All %d spatial graphs have", topology$n_graphs)
    }
    message <- sprintf(
      paste0("%s no edges, so the resolvent is the identity and positive ",
             "lambda = %s will have no effect.%s"),
      who, format(topology$lambda),
      .dkge_spatial_topology_hint(source, dthresh)
    )
    subclass <- "dkge_spatial_inert_warning"
  } else {
    message <- sprintf(
      paste0("%d of %d spatial graphs (%s) have no edges, so their ",
             "resolvents are the identity and positive lambda = %s will have no ",
             "effect for those domains.%s"),
      length(topology$empty), topology$n_graphs,
      paste(topology$empty_labels, collapse = ", "),
      format(topology$lambda),
      .dkge_spatial_topology_hint(source, dthresh)
    )
    subclass <- "dkge_spatial_partial_warning"
  }
  .dkge_warn(message, subclass)
  invisible(NULL)
}

#' @export
print.dkge_spatial_regularizer <- function(x, ...) {
  payload <- .dkge_spatial_spec_payload(x, validate = TRUE)
  topology <- payload$topology
  cat("<dkge_spatial_regularizer>\n")
  cat("  source   :", x$source, "\n")
  cat("  domains  :", if (isTRUE(x$shared)) "shared" else length(x$laplacians), "\n")
  cat("  lambda   :", format(x$lambda), "\n")
  cat("  status   :", topology$status, "\n")
  cat("  units    :", paste(topology$n_units, collapse = ", "), "\n")
  cat("  edges    :", paste(format(topology$n_edges, trim = TRUE), collapse = ", "), "\n")
  if (identical(topology$status, "inert")) {
    cat("  NOTE     : edgeless graphs give H = I; lambda has no effect\n")
  } else if (identical(topology$status, "partial")) {
    cat("  NOTE     : edgeless domains:",
        paste(topology$empty_labels, collapse = ", "), "\n")
  }
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
  edge_nnz <- .dkge_spatial_edge_count(L)
  requested <- lambda > 0
  effective <- requested && edge_nnz > 0
  list(
    L = L,
    A = A,
    factor = factor,
    lambda = lambda,
    requested = requested,
    effective = effective,
    status = if (!requested) "inactive" else if (effective) "active" else "inert",
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

#' Canonical payload for the live spatial solve operator
#'
#' The construction-time `fingerprint` identifies the Laplacian geometry, but
#' it is not sufficient provenance for inference: the fitted estimator applies
#' the cached CHOLMOD factor directly. Bind all live numerical state and check
#' that the factor still represents `I + lambda * L` before it can enter an
#' inferential alignment receipt.
#'
#' @keywords internal
#' @noRd
.dkge_spatial_operator_payload <- function(operator, validate = TRUE,
                                            tolerance = 1e-8) {
  if (is.null(operator)) return(NULL)
  if (!is.list(operator)) {
    .dkge_abort("A spatial operator must be a list.",
                "dkge_spatial_factor_error")
  }
  lambda <- operator$lambda
  n_units <- operator$n_units
  L <- operator$L
  if (!is.numeric(lambda) || length(lambda) != 1L || !is.finite(lambda) ||
      lambda < 0 || !is.numeric(n_units) || length(n_units) != 1L ||
      !is.finite(n_units) || n_units < 1L || n_units != as.integer(n_units) ||
      is.null(L) || !identical(dim(L), rep(as.integer(n_units), 2L))) {
    .dkge_abort("The live spatial operator has invalid dimensions or lambda.",
                "dkge_spatial_factor_error")
  }

  edge_nnz <- as.integer(.dkge_spatial_edge_count(L))
  domain <- .dkge_validate_spatial_domain(
    operator$domain, as.integer(n_units), "live spatial operator domain"
  )
  domain_mode_expected <- if (is.null(domain)) {
    "positional"
  } else {
    c("declared_positional", "label_matched")
  }
  domain_dimnames_valid <- if (is.null(domain)) {
    is.null(rownames(L)) && is.null(colnames(L))
  } else {
    identical(rownames(L), domain) && identical(colnames(L), domain)
  }
  construction_fingerprint <- digest::digest(
    list(dim = dim(L), domain = domain, entries = Matrix::summary(L)),
    algo = "xxhash64", serialize = TRUE
  )
  requested_expected <- lambda > 0
  effective_expected <- requested_expected && edge_nnz > 0L
  status_expected <- if (!requested_expected) {
    "inactive"
  } else if (effective_expected) {
    "active"
  } else {
    "inert"
  }
  if (isTRUE(validate) &&
      (!is.numeric(operator$n_edges) || length(operator$n_edges) != 1L ||
       !is.finite(operator$n_edges) ||
       operator$n_edges != as.integer(operator$n_edges) ||
       !identical(as.integer(operator$n_edges), edge_nnz) ||
       !identical(operator$requested, requested_expected) ||
       !identical(operator$effective, effective_expected) ||
       !identical(operator$status, status_expected) ||
       !is.character(operator$domain_mode) ||
       length(operator$domain_mode) != 1L ||
       !operator$domain_mode %in% domain_mode_expected ||
       !isTRUE(domain_dimnames_valid) ||
       !is.character(operator$subject_id) ||
       length(operator$subject_id) != 1L || is.na(operator$subject_id) ||
       !nzchar(operator$subject_id) ||
       !is.character(operator$source) || length(operator$source) != 1L ||
       is.na(operator$source) || !nzchar(operator$source) ||
       !identical(operator$fingerprint, construction_fingerprint))) {
    .dkge_abort(
      paste0(
        "The live spatial operator effectiveness metadata is inconsistent ",
        "with its Laplacian and lambda."
      ),
      "dkge_spatial_provenance_error"
    )
  }

  factor_expansion <- NULL
  if (lambda == 0) {
    if (isTRUE(validate) &&
        (!is.null(operator$A) || !is.null(operator$factor))) {
      .dkge_abort(
        "A zero-lambda spatial operator must not carry an active solve factor.",
        "dkge_spatial_factor_error"
      )
    }
  } else {
    A_expected <- Matrix::forceSymmetric(
      Matrix::Diagonal(as.integer(n_units)) + lambda * L,
      uplo = "U"
    )
    if (is.null(operator$A) || is.null(operator$factor)) {
      .dkge_abort("An active spatial operator is missing its matrix or factor.",
                  "dkge_spatial_factor_error")
    }
    if (isTRUE(validate)) {
      a_scale <- max(as.numeric(Matrix::norm(A_expected, "F")),
                     .Machine$double.xmin)
      a_error <- as.numeric(Matrix::norm(operator$A - A_expected, "F")) /
        a_scale
      factor_expansion <- tryCatch(
        Matrix::expand(operator$factor),
        error = function(e) NULL
      )
      if (!is.list(factor_expansion) ||
          is.null(factor_expansion$P) || is.null(factor_expansion$L)) {
        .dkge_abort("The live spatial factor cannot be expanded.",
                    "dkge_spatial_factor_error")
      }
      reconstructed <- Matrix::tcrossprod(factor_expansion$L)
      permuted_A <- Matrix::tcrossprod(
        factor_expansion$P %*% operator$A,
        factor_expansion$P
      )
      factor_error <- as.numeric(Matrix::norm(
        permuted_A - reconstructed, "F"
      )) / max(as.numeric(Matrix::norm(permuted_A, "F")),
               .Machine$double.xmin)
      if (!is.finite(a_error) || a_error > tolerance ||
          !is.finite(factor_error) || factor_error > tolerance) {
        .dkge_abort(
          paste0(
            "The live spatial factor is inconsistent with `I + lambda * L` ",
            "(matrix error ", format(a_error, digits = 4L),
            ", factor error ", format(factor_error, digits = 4L), ")."
          ),
          "dkge_spatial_factor_error"
        )
      }
    }
    if (is.null(factor_expansion)) {
      factor_expansion <- Matrix::expand(operator$factor)
    }
  }

  list(
    L = operator$L,
    A = operator$A,
    # Include the exact serialized object that Matrix::solve() consumes, plus
    # its canonical permutation/triangular expansion for auditability.
    factor = operator$factor,
    factor_expansion = factor_expansion,
    lambda = as.numeric(operator$lambda),
    requested = requested_expected,
    effective = effective_expected,
    status = status_expected,
    n_units = as.integer(operator$n_units),
    n_edges = edge_nnz,
    domain = domain,
    domain_mode = operator$domain_mode,
    subject_id = operator$subject_id,
    source = operator$source,
    construction_fingerprint = construction_fingerprint
  )
}

#' Canonical payload for a fitted spatial regularizer
#'
#' Effectiveness flags are control flow, not decorative diagnostics: they
#' decide whether the fitted estimator applies a solve. Recompute them from
#' the live Laplacians and lambda before any early return, and bind the entire
#' fitted spatial object into downstream inferential receipts.
#'
#' @keywords internal
#' @noRd
.dkge_spatial_fit_payload <- function(spatial, validate = TRUE) {
  if (is.null(spatial)) return(NULL)
  if (!inherits(spatial, "dkge_spatial_fit") ||
      !is.list(spatial$operators) || !length(spatial$operators)) {
    .dkge_abort("The fitted spatial state is malformed.",
                "dkge_spatial_provenance_error")
  }
  .dkge_spatial_spec_payload(spatial$spec, validate = validate)
  operators <- lapply(
    spatial$operators, .dkge_spatial_operator_payload, validate = validate
  )
  names(operators) <- names(spatial$operators)
  operator_names <- names(operators)
  if (isTRUE(validate) &&
      (is.null(operator_names) || anyNA(operator_names) ||
       any(!nzchar(operator_names)) || anyDuplicated(operator_names) ||
       !identical(
         unname(vapply(operators, `[[`, character(1), "subject_id")),
         operator_names
       ))) {
    .dkge_abort(
      "Fitted spatial operator names do not match their subject identities.",
      "dkge_spatial_provenance_error"
    )
  }
  lambda <- .dkge_spatial_lambda(spatial$lambda)
  if (isTRUE(validate) && any(vapply(
    operators,
    function(op) !identical(op$lambda, as.numeric(lambda)),
    logical(1)
  ))) {
    .dkge_abort(
      "Fitted spatial operators do not share the recorded top-level lambda.",
      "dkge_spatial_provenance_error"
    )
  }
  topology <- .dkge_spatial_topology(
    lapply(operators, `[[`, "L"), lambda
  )
  if (isTRUE(validate)) {
    if (!inherits(spatial$spec, "dkge_spatial_regularizer")) {
      .dkge_abort("The fitted spatial specification is malformed.",
                  "dkge_spatial_provenance_error")
    }
    spec_laplacians <- if (isTRUE(spatial$shared)) {
      rep(spatial$spec$laplacians, length(operators))
    } else {
      .dkge_spatial_named_order(
        spatial$spec$laplacians, operator_names, "spatial Laplacians"
      )
    }
    spec_domains <- if (isTRUE(spatial$shared)) {
      rep(spatial$spec$domains, length(operators))
    } else {
      .dkge_spatial_named_order(
        spatial$spec$domains, operator_names, "spatial domains"
      )
    }
    spec_matches <- Map(function(op, L_spec, domain_spec) {
      domain_spec <- .dkge_validate_spatial_domain(
        domain_spec, nrow(L_spec), "fitted spatial specification domain"
      )
      if (is.null(domain_spec)) {
        L_expected <- L_spec
        domain_matches <- is.null(op$domain)
      } else {
        domain_matches <- !is.null(op$domain) &&
          setequal(op$domain, domain_spec)
        if (!domain_matches) return(FALSE)
        idx <- match(op$domain, domain_spec)
        L_expected <- L_spec[idx, idx, drop = FALSE]
        dimnames(L_expected) <- list(op$domain, op$domain)
      }
      isTRUE(domain_matches) && identical(op$L, L_expected) &&
        identical(op$source, spatial$spec$source)
    }, operators, spec_laplacians, spec_domains)
    if (!all(unlist(spec_matches, use.names = FALSE))) {
      .dkge_abort(
        paste0(
          "A fitted spatial operator does not match the bound ",
          "specification domain, Laplacian, or source."
        ),
        "dkge_spatial_provenance_error"
      )
    }
  }
  expected_top <- list(
    active = topology$effective,
    requested = topology$requested,
    effective = topology$effective,
    fully_effective = topology$fully_effective,
    status = topology$status,
    lambda = topology$lambda
  )
  observed_top <- unclass(spatial)[names(expected_top)]
  if (isTRUE(validate) && !identical(observed_top, expected_top)) {
    .dkge_abort(
      paste0(
        "The fitted spatial effectiveness summary is inconsistent with ",
        "its live operators."
      ),
      "dkge_spatial_provenance_error"
    )
  }

  diagnostics <- do.call(rbind, lapply(operators, function(op) {
    data.frame(
      subject = op$subject_id,
      n_units = op$n_units,
      n_edges = op$n_edges,
      lambda = op$lambda,
      requested = op$requested,
      effective = op$effective,
      status = op$status,
      source = op$source,
      domain_mode = op$domain_mode,
      fingerprint = op$construction_fingerprint,
      stringsAsFactors = FALSE
    )
  }))
  rownames(diagnostics) <- NULL
  if (isTRUE(validate) && !identical(spatial$diagnostics, diagnostics)) {
    .dkge_abort(
      "The fitted spatial diagnostics are inconsistent with the live operators.",
      "dkge_spatial_provenance_error"
    )
  }

  subject_ids <- unname(vapply(
    operators, `[[`, character(1), "subject_id"
  ))
  expected_provenance <- list(
    method = "laplacian_resolvent",
    formula = "H = solve(I + lambda * L)",
    lambda = lambda,
    status = topology$status,
    requested = topology$requested,
    effective = topology$effective,
    fully_effective = topology$fully_effective,
    shared = isTRUE(spatial$shared),
    source = spatial$spec$source,
    construction = spatial$spec$construction,
    subjects = subject_ids,
    domains = lapply(operators, function(op) {
      list(
        n_units = op$n_units,
        n_edges = op$n_edges,
        effective = op$effective,
        status = op$status,
        domain_mode = op$domain_mode,
        fingerprint = op$construction_fingerprint
      )
    })
  )
  names(expected_provenance$domains) <- names(operators)
  if (isTRUE(validate) &&
      (!inherits(spatial$spec, "dkge_spatial_regularizer") ||
       !identical(spatial$spec$lambda, lambda) ||
       !identical(isTRUE(spatial$spec$shared), isTRUE(spatial$shared)) ||
       !identical(spatial$provenance, expected_provenance))) {
    .dkge_abort(
      "The fitted spatial specification or provenance is inconsistent.",
      "dkge_spatial_provenance_error"
    )
  }

  list(
    active = spatial$active,
    requested = spatial$requested,
    effective = spatial$effective,
    fully_effective = spatial$fully_effective,
    status = spatial$status,
    lambda = spatial$lambda,
    shared = spatial$shared,
    operators = operators,
    diagnostics = spatial$diagnostics,
    provenance = spatial$provenance,
    spec = spatial$spec
  )
}

.dkge_spatial_operator_binding <- function(operator, validate = TRUE) {
  .dkge_object_hash(.dkge_spatial_operator_payload(
    operator, validate = validate
  ))
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
  .dkge_spatial_spec_payload(spatial, validate = TRUE)
  topology <- .dkge_spatial_topology(spatial$laplacians, spatial$lambda)
  if (topology$status %in% c("inert", "partial") &&
      !.dkge_spatial_warning_is_current(spatial, topology)) {
    .dkge_spatial_warn_topology(
      topology,
      spatial$source,
      spatial$construction$dthresh %||% NULL
    )
  }
  spatial$topology <- topology
  spatial["topology_warning"] <- list(
    if (topology$status %in% c("inert", "partial")) {
      .dkge_spatial_warning_receipt(topology)
    } else {
      NULL
    }
  )
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
      requested = op$requested,
      effective = op$effective,
      status = op$status,
      source = op$source,
      domain_mode = op$domain_mode,
      fingerprint = op$fingerprint,
      stringsAsFactors = FALSE
    )
  }))
  rownames(diagnostics) <- NULL
  fitted_topology <- .dkge_spatial_topology(
    lapply(operators, `[[`, "L"), spatial$lambda
  )

  provenance <- list(
    method = "laplacian_resolvent",
    formula = "H = solve(I + lambda * L)",
    lambda = spatial$lambda,
    status = fitted_topology$status,
    requested = fitted_topology$requested,
    effective = fitted_topology$effective,
    fully_effective = fitted_topology$fully_effective,
    shared = isTRUE(spatial$shared),
    source = spatial$source,
    construction = spatial$construction,
    subjects = subject_ids,
    domains = lapply(operators, function(op) {
      list(
        n_units = op$n_units,
        n_edges = op$n_edges,
        effective = op$effective,
        status = op$status,
        domain_mode = op$domain_mode,
        fingerprint = op$fingerprint
      )
    })
  )

  structure(
    list(
      active = fitted_topology$effective,
      requested = fitted_topology$requested,
      effective = fitted_topology$effective,
      fully_effective = fitted_topology$fully_effective,
      status = fitted_topology$status,
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
  if (is.null(operator)) {
    return(B)
  }
  payload <- .dkge_spatial_operator_payload(operator, validate = TRUE)
  if (payload$lambda == 0 || !isTRUE(payload$effective) ||
      identical(payload$n_edges, 0L)) {
    return(B)
  }
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
  if (is.null(spatial)) return(NULL)
  .dkge_spatial_fit_payload(spatial, validate = TRUE)
  if (!isTRUE(spatial$active)) return(NULL)
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

.dkge_spatial_with_lambda <- function(spatial, lambda,
                                      topology_checked = FALSE) {
  if (inherits(spatial, "dkge_spatial_fit")) spatial <- spatial$spec
  if (!inherits(spatial, "dkge_spatial_regularizer")) {
    .dkge_abort("`spatial` must be a DKGE spatial regularizer.",
                "dkge_spatial_spec_error")
  }
  .dkge_spatial_spec_payload(spatial, validate = TRUE)
  spatial$lambda <- .dkge_spatial_lambda(lambda)
  topology <- .dkge_spatial_topology(spatial$laplacians, spatial$lambda)
  warning_current <- .dkge_spatial_warning_is_current(spatial, topology)
  spatial$topology <- topology
  spatial["topology_warning"] <- list(
    if (topology$status %in% c("inert", "partial") &&
        (isTRUE(topology_checked) || warning_current)) {
      .dkge_spatial_warning_receipt(topology)
    } else {
      NULL
    }
  )
  .dkge_spatial_spec_payload(spatial, validate = TRUE)
  spatial
}

.dkge_prediction_spatial <- function(object, B_list, spatial = NULL) {
  fitted_spatial <- object$spatial %||% NULL
  fitted_topology <- NULL
  if (!is.null(fitted_spatial) && length(fitted_spatial$operators)) {
    .dkge_spatial_fit_payload(fitted_spatial, validate = TRUE)
    fitted_topology <- list(effective = fitted_spatial$effective)
  } else if (!is.null(fitted_spatial)) {
    .dkge_spatial_spec_payload(fitted_spatial$spec, validate = TRUE)
    fitted_topology <- .dkge_spatial_topology(
      fitted_spatial$spec$laplacians, fitted_spatial$spec$lambda
    )
    expected_top <- list(
      active = fitted_topology$effective,
      requested = fitted_topology$requested,
      effective = fitted_topology$effective,
      fully_effective = fitted_topology$fully_effective,
      status = fitted_topology$status,
      lambda = fitted_topology$lambda,
      shared = isTRUE(fitted_spatial$spec$shared)
    )
    if (!identical(
      unclass(fitted_spatial)[names(expected_top)], expected_top
    )) {
      .dkge_abort(
        "The frozen spatial prediction summary is inconsistent with its specification.",
        "dkge_spatial_provenance_error"
      )
    }
  }
  active_fit <- !is.null(fitted_topology) &&
    isTRUE(fitted_topology$effective)
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
