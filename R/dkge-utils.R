# dkge-utils.R
# Shared helper utilities for DKGE

#' Null-coalescing helper
#'
#' Returns `b` when `a` is `NULL`, otherwise returns `a`.
#'
#' @name grapes-or-or-grapes
#' @keywords internal
NULL

#' @rdname grapes-or-or-grapes
#' @param a Primary value tested for `NULL`.
#' @param b Fallback value returned when `a` is `NULL`.
#' @usage a \%||\% b
#' @keywords internal
`%||%` <- function(a, b) if (is.null(a)) b else a

#' Signal a stable DKGE error condition
#'
#' Base condition subclasses let callers distinguish input-contract failures
#' without making message text or a new dependency part of the API.
#'
#' @param message User-facing error message.
#' @param subclass Specific condition subclass.
#' @keywords internal
#' @noRd
.dkge_abort <- function(message, subclass = "dkge_input_error") {
  condition <- structure(
    list(message = as.character(message), call = NULL),
    class = c(subclass, "dkge_input_error", "dkge_error", "error", "condition")
  )
  stop(condition)
}

#' Validate one scalar against a named numeric domain
#'
#' Public scalar contracts must be checked before integer coercion can truncate
#' fractional values or vector coercion can silently discard shape information.
#'
#' @keywords internal
#' @noRd
.dkge_validate_scalar_domain <- function(x, arg, domain, predicate,
                                         integer_result = FALSE) {
  supplied <- if (!length(x)) {
    "<empty>"
  } else {
    values <- paste(utils::head(as.character(x), 5L), collapse = ", ")
    if (length(x) > 5L) paste0(values, ", ...") else values
  }
  valid <- is.numeric(x) && is.null(dim(x)) && length(x) == 1L &&
    is.finite(x) && isTRUE(predicate(x))
  if (!valid) {
    .dkge_abort(
      sprintf(
        "`%s` supplied length %d value(s): %s; expected %s.",
        arg, length(x), supplied, domain
      ),
      "dkge_validation_error"
    )
  }
  if (integer_result) as.integer(x) else as.numeric(x)
}

#' @keywords internal
#' @noRd
.dkge_validate_positive_scalar <- function(x, arg) {
  .dkge_validate_scalar_domain(
    x, arg, "a finite positive scalar", function(value) value > 0
  )
}

#' @keywords internal
#' @noRd
.dkge_validate_nonnegative_scalar <- function(x, arg) {
  .dkge_validate_scalar_domain(
    x, arg, "a finite non-negative scalar", function(value) value >= 0
  )
}

#' @keywords internal
#' @noRd
.dkge_validate_integer_scalar <- function(x, arg) {
  .dkge_validate_scalar_domain(
    x, arg, "a finite integer-valued scalar",
    function(value) value == trunc(value), integer_result = TRUE
  )
}

#' @keywords internal
#' @noRd
.dkge_validate_positive_integer <- function(x, arg) {
  .dkge_validate_scalar_domain(
    x, arg, "a strictly positive integer",
    function(value) value > 0 && value == trunc(value),
    integer_result = TRUE
  )
}

#' @keywords internal
#' @noRd
.dkge_validate_probability <- function(x, arg) {
  .dkge_validate_scalar_domain(
    x, arg, "a probability in the closed interval [0, 1]",
    function(value) value >= 0 && value <= 1
  )
}

#' Signal a stable DKGE warning condition
#'
#' @param message User-facing warning message.
#' @param subclass Specific condition subclass.
#' @keywords internal
#' @noRd
.dkge_warn <- function(message, subclass = "dkge_warning") {
  condition <- structure(
    list(message = as.character(message), call = NULL),
    class = unique(c(subclass, "dkge_warning", "warning", "condition"))
  )
  warning(condition)
}

#' Validate optional design-kernel metadata
#'
#' @param info Optional metadata associated with a design kernel.
#' @return `info`, invisibly validated as a list or `NULL`.
#' @keywords internal
#' @noRd
.dkge_validate_kernel_info <- function(info) {
  if (is.null(info)) {
    return(NULL)
  }
  if (!is.list(info)) {
    .dkge_abort(
      paste0(
        "Design-kernel metadata `info` must be a list or NULL; got class ",
        paste(class(info), collapse = "/"), "."
      ),
      "dkge_kernel_info_error"
    )
  }
  if (!is.null(info$info) && !is.list(info$info)) {
    .dkge_abort(
      "Nested design-kernel metadata `info$info` must be a list or NULL.",
      "dkge_kernel_info_error"
    )
  }
  info
}

#' Apply helper with optional parallelism
#'
#' Wraps `lapply()` with an optional future.apply backend so callers can enable
#' `parallel = TRUE` without repeating boilerplate dependency checks.
#'
#' @param X Vector or list to iterate over.
#' @param FUN Function to apply.
#' @param parallel Logical; if `TRUE`, uses `future.apply::future_lapply()`.
#' @param ... Additional arguments passed to the apply backend.
#' @return List of results matching `lapply()` semantics.
#' @keywords internal
.dkge_apply <- function(X, FUN, parallel = FALSE, ...) {
  if (parallel) {
    if (!requireNamespace("future.apply", quietly = TRUE)) {
      stop("parallel=TRUE requires the future.apply package; install it or set parallel=FALSE.",
           call. = FALSE)
    }
    future.apply::future_lapply(X, FUN, ...)
  } else {
    lapply(X, FUN, ...)
  }
}

#' Check whether verbose output should be emitted
#'
#' Uses the per-call `verbose` flag combined with the global
#' `options(dkge.verbose = TRUE)` toggle.
#'
#' @keywords internal
.dkge_verbose <- function(verbose) {

  isTRUE(verbose) && isTRUE(getOption("dkge.verbose", TRUE))
}

# -------------------------------------------------------------------------
# Numerical robustness utilities ------------------------------------------
# -------------------------------------------------------------------------

.dkge_kernel_validation_cache <- new.env(parent = emptyenv())

#' Validate design kernel matrix
#'
#' Enforces finite numeric entries, symmetry (up to tolerance), and positive
#' semidefiniteness. Small numerical asymmetry is corrected by symmetrization.
#'
#' @param K Candidate kernel matrix.
#' @param tol Relative tolerance for symmetry/PSD checks.
#' @return Symmetrized kernel matrix.
#' @keywords internal
#' @noRd
.dkge_validate_kernel <- function(K, tol = 1e-8) {
  if (!is.matrix(K) || !is.numeric(K)) {
    stop("Kernel `K` must be a numeric matrix.", call. = FALSE)
  }
  if (nrow(K) != ncol(K)) {
    stop("Kernel `K` must be square.", call. = FALSE)
  }
  if (any(!is.finite(K))) {
    stop("Kernel `K` contains non-finite values.", call. = FALSE)
  }

  Ksym <- (K + t(K)) / 2
  asym <- max(abs(K - t(K)))
  k_scale <- max(1, max(abs(K)))
  if (asym > tol * k_scale) {
    stop(sprintf(
      "Kernel `K` must be symmetric (max asymmetry %.3e exceeds tolerance %.3e).",
      asym, tol * k_scale
    ), call. = FALSE)
  }
  if (asym > 1e-12 * k_scale) {
    warning("Kernel `K` is not exactly symmetric; using (K + t(K)) / 2.", call. = FALSE)
  }

  # Exact symmetric kernels recur throughout fold alignment and consensus.
  # Hash every entry before reuse, so mutation cannot inherit a stale PSD
  # verdict; only the cubic eigendecomposition is skipped. Slightly asymmetric
  # inputs are not cached so their corrective warning remains visible.
  cache_key <- NULL
  if (asym == 0) {
    cache_key <- digest::digest(list(K = Ksym, tol = tol),
                                algo = "xxhash64", serialize = TRUE)
    cached <- .dkge_kernel_validation_cache[[cache_key]]
    if (!is.null(cached)) {
      return(cached)
    }
  }

  eig_vals <- eigen(Ksym, symmetric = TRUE, only.values = TRUE)$values
  eig_scale <- max(1, max(abs(eig_vals)))
  neg_tol <- tol * eig_scale
  min_eig <- min(eig_vals)
  if (min_eig < -neg_tol) {
    stop(sprintf(
      "Kernel `K` must be positive semidefinite; minimum eigenvalue %.3e is below tolerance %.3e.",
      min_eig, -neg_tol
    ), call. = FALSE)
  }
  if (min_eig < -neg_tol / 10) {
    warning("Kernel `K` has small negative eigenvalues; they will be clamped in kernel roots.", call. = FALSE)
  }

  if (!is.null(cache_key)) {
    assign(cache_key, Ksym, envir = .dkge_kernel_validation_cache)
    keys <- ls(.dkge_kernel_validation_cache, all.names = TRUE)
    if (length(keys) > 32L) {
      rm(list = keys[seq_len(length(keys) - 32L)],
         envir = .dkge_kernel_validation_cache)
    }
  }

  Ksym
}

#' Exact positive-semidefinite kernel geometry
#'
#' Computes square and Moore--Penrose inverse square roots without adding
#' energy to the null space. The numerical support is defined by a relative
#' eigentolerance, so rank is invariant to positive rescaling of the kernel.
#'
#' @param K Finite symmetric positive-semidefinite matrix.
#' @param tol Relative eigentolerance used to define the kernel support.
#' @return Kernel roots, support projectors, eigenstructure, and scalar rank
#'   diagnostics.
#' @keywords internal
#' @noRd
.dkge_kernel_geometry <- function(K, tol = 1e-10) {
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) ||
      tol < 0 || tol >= 1) {
    stop("`tol` must be a finite scalar in [0, 1).", call. = FALSE)
  }

  Ksym <- .dkge_validate_kernel(K)
  ee <- eigen(Ksym, symmetric = TRUE)
  vals_raw <- ee$values
  vals <- pmax(vals_raw, 0)
  spectral_scale <- if (length(vals)) max(vals) else 0
  abs_tol <- tol * spectral_scale
  positive <- if (spectral_scale > 0) vals > abs_tol else rep(FALSE, length(vals))
  vals_support <- ifelse(positive, vals, 0)

  sqrt_vals <- sqrt(vals_support)
  inv_sqrt_vals <- numeric(length(vals_support))
  inv_sqrt_vals[positive] <- 1 / sqrt_vals[positive]
  V <- ee$vectors
  n <- length(vals_support)
  Khalf <- V %*% diag(sqrt_vals, n) %*% t(V)
  Kihalf <- V %*% diag(inv_sqrt_vals, n) %*% t(V)
  V_support <- V[, positive, drop = FALSE]
  support_projector <- if (any(positive)) {
    tcrossprod(V_support)
  } else {
    matrix(0, n, n)
  }
  null_projector <- diag(1, n) - support_projector

  dimnames(Khalf) <- dimnames(Ksym)
  dimnames(Kihalf) <- dimnames(Ksym)
  dimnames(support_projector) <- dimnames(Ksym)
  dimnames(null_projector) <- dimnames(Ksym)

  rank <- sum(positive)
  condition <- if (rank > 0L) {
    max(vals[positive]) / min(vals[positive])
  } else {
    Inf
  }
  spectral_mass <- sum(vals_support)
  effective_rank_pr <- if (spectral_mass > 0) {
    spectral_mass^2 / sum(vals_support^2)
  } else {
    0
  }
  # Clamp only round-off excursions: the participation ratio is bounded by
  # the numerical support rank and is invariant to positive kernel rescaling.
  effective_rank_pr <- pmin(as.numeric(rank), pmax(0, effective_rank_pr))
  effective_rank_fraction <- if (rank > 0L) effective_rank_pr / rank else 0
  leading_eigenvalue_share <- if (spectral_mass > 0) {
    max(vals_support) / spectral_mass
  } else {
    NA_real_
  }
  near_singular <- rank == n && is.finite(condition) && condition >= 1e8

  list(
    K = Ksym,
    Khalf = Khalf,
    Kihalf = Kihalf,
    evals = vals_support,
    evals_raw = vals_raw,
    evecs = V,
    support = positive,
    support_projector = support_projector,
    null_projector = null_projector,
    rank = as.integer(rank),
    nullity = as.integer(n - rank),
    condition = as.numeric(condition),
    effective_rank_pr = as.numeric(effective_rank_pr),
    effective_rank_fraction = as.numeric(effective_rank_fraction),
    leading_eigenvalue_share = as.numeric(leading_eigenvalue_share),
    tolerance = as.numeric(abs_tol),
    relative_tolerance = tol,
    full_rank = rank == n,
    near_singular = near_singular,
    status = if (rank < n) "singular" else if (near_singular) "ill_conditioned" else "well_conditioned"
  )
}

#' Scale-equivariant spectral rank contract
#'
#' The smallest positive normal double protects the exactly-zero case without
#' imposing a fixed data scale. All transformed-moment rank decisions use the
#' same relative threshold, so multiplying betas by a positive constant cannot
#' change the selected rank while the spectrum remains representable.
#'
#' @param values Finite numeric spectrum.
#' @param absolute_tolerance Non-negative absolute tolerance.
#' @param relative_tolerance Non-negative tolerance relative to spectral scale.
#' @return Applied tolerance, scale, positive mask, and numerical rank.
#' @keywords internal
#' @noRd
.dkge_spectral_contract <- function(
    values,
    absolute_tolerance = .Machine$double.xmin,
    relative_tolerance = 1e-8) {
  if (!is.numeric(values) || any(!is.finite(values))) {
    .dkge_abort(
      "A finite numeric spectrum is required for rank selection.",
      "dkge_spectral_error"
    )
  }
  if (!is.numeric(absolute_tolerance) || length(absolute_tolerance) != 1L ||
      !is.finite(absolute_tolerance) || absolute_tolerance < 0 ||
      !is.numeric(relative_tolerance) || length(relative_tolerance) != 1L ||
      !is.finite(relative_tolerance) || relative_tolerance < 0) {
    .dkge_abort(
      "Spectral absolute and relative tolerances must be finite non-negative scalars.",
      "dkge_spectral_error"
    )
  }
  scale <- if (length(values)) max(abs(values)) else 0
  tolerance <- absolute_tolerance + relative_tolerance * scale
  positive <- values > tolerance
  list(
    absolute_tolerance = absolute_tolerance,
    relative_tolerance = relative_tolerance,
    scale = scale,
    tolerance = tolerance,
    positive = positive,
    rank = as.integer(sum(positive))
  )
}

#' Compact public-facing kernel diagnostics
#'
#' @param geometry Result from `.dkge_kernel_geometry()`.
#' @return Scalar diagnostic fields only.
#' @keywords internal
#' @noRd
.dkge_kernel_diagnostics <- function(geometry) {
  geometry[c(
    "rank", "nullity", "condition", "effective_rank_pr",
    "effective_rank_fraction", "leading_eigenvalue_share", "tolerance",
    "relative_tolerance", "full_rank", "near_singular", "status"
  )]
}

#' Check matrix rank for design and/or beta matrices
#'
#' Detects rank deficiency and emits informative warnings identifying the
#' culprit subject and the nature of the problem.
#'
#' @param design Design matrix to check (T x q).
#' @param beta Optional beta matrix to check (q x P).
#' @param subject_id Optional subject identifier for warning messages.
#' @return List with `design_rank` and `beta_rank` (if beta provided).
#' @keywords internal
#' @noRd
.dkge_check_rank <- function(design, beta = NULL, subject_id = NULL) {
  subject_label <- subject_id %||% "(unnamed)"
  result <- list(design_rank = NULL, beta_rank = NULL)


  if (!is.null(design) && is.matrix(design)) {
    qr_design <- qr(design)
    design_rank <- qr_design$rank
    expected_rank <- ncol(design)
    result$design_rank <- design_rank

    if (design_rank < expected_rank) {
      warning(sprintf(
        "Subject '%s': design matrix is rank-deficient (rank %d < %d columns). Effects may be aliased.",
        subject_label, design_rank, expected_rank
      ), call. = FALSE)
    }
  }

  if (!is.null(beta) && is.matrix(beta)) {
    beta_rank <- qr(beta)$rank
    result$beta_rank <- beta_rank

    if (beta_rank < nrow(beta)) {
      warning(sprintf(
        "Subject '%s': beta matrix has reduced rank (%d < %d effects).",
        subject_label, beta_rank, nrow(beta)
      ), call. = FALSE)
    }
  }

  result
}

#' Check matrix condition number against threshold
#'
#' Warns if the condition number exceeds the specified threshold, indicating
#' potential numerical instability.
#'
#' @param M Symmetric matrix to check.
#' @param threshold Condition number threshold (default 1e8).
#' @param name Descriptive name for the matrix (used in warning message).
#' @return The computed condition number.
#' @keywords internal
#' @noRd
.dkge_check_condition <- function(M, threshold = 1e8, name = "matrix") {
  cond <- kappa(M, exact = FALSE)
  if (cond > threshold) {
    warning(sprintf(
      "%s is ill-conditioned (condition number: %.2e > %.2e threshold). Results may be numerically unstable.",
      name, cond, threshold
    ), call. = FALSE)
  }
  cond
}

#' Identify and track voxels with non-finite values
#'
#' Scans each beta matrix for columns containing NA, NaN, or Inf values,
#' emits per-subject warnings, and returns metadata about exclusions.
#'
#' @param B_list List of beta matrices (q x P_s each).
#' @param subject_ids Optional character vector of subject identifiers.
#' @return List with `excluded_voxels` (list of integer vectors per subject),
#'   `excluded_counts` (integer vector), and `total_excluded` (integer).
#' @keywords internal
#' @noRd
.dkge_voxel_exclusion_mask <- function(B_list, subject_ids = NULL) {
  S <- length(B_list)
  excluded_voxels <- vector("list", S)
  excluded_counts <- integer(S)

  for (s in seq_len(S)) {
    B <- B_list[[s]]
    if (is.null(B) || !is.matrix(B) || ncol(B) == 0) {
      excluded_voxels[[s]] <- integer(0)
      excluded_counts[s] <- 0L
      next
    }

    bad_cols <- which(colSums(!is.finite(B)) > 0)
    excluded_voxels[[s]] <- bad_cols
    excluded_counts[s] <- length(bad_cols)

    if (length(bad_cols) > 0) {
      pct <- 100 * length(bad_cols) / ncol(B)
      subject_label <- if (!is.null(subject_ids) && length(subject_ids) >= s) {
        subject_ids[s]
      } else {
        as.character(s)
      }
      warning(sprintf(
        "Subject '%s': %d voxels (%.1f%%) excluded due to NA/NaN/Inf values.",
        subject_label, length(bad_cols), pct
      ), call. = FALSE)
    }
  }

  list(
    excluded_voxels = excluded_voxels,
    excluded_counts = excluded_counts,
    total_excluded = sum(excluded_counts)
  )
}

# -------------------------------------------------------------------------
# Resampling helpers ------------------------------------------------------
# -------------------------------------------------------------------------

#' Validate a resampling replicate count
#'
#' Shared by the between-subject resampling entry points so that `B` is
#' rejected identically everywhere.
#'
#' @param B Candidate number of replicates.
#' @return `B` coerced to a positive integer scalar.
#' @keywords internal
#' @noRd
.dkge_validate_resample_B <- function(B) {
  .dkge_validate_positive_integer(B, "B")
}

#' Enter a seeded RNG scope
#'
#' Records the caller's `.Random.seed` (or its absence) and seeds the stream.
#' Pair with `.dkge_seed_exit()` via `on.exit()` so that a seeded run leaves the
#' caller's RNG state exactly as it found it. A `NULL` seed is a no-op.
#'
#' @param seed Optional seed passed to [set.seed()].
#' @return Opaque state to hand back to `.dkge_seed_exit()`.
#' @keywords internal
#' @noRd
.dkge_seed_enter <- function(seed) {
  if (is.null(seed)) {
    return(list(active = FALSE))
  }
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
  set.seed(seed)
  list(active = TRUE, old_seed = old_seed)
}

#' Leave a seeded RNG scope
#'
#' @param state Value returned by `.dkge_seed_enter()`.
#' @return `NULL`, invisibly.
#' @keywords internal
#' @noRd
.dkge_seed_exit <- function(state) {
  if (!isTRUE(state$active)) {
    return(invisible(NULL))
  }
  if (is.null(state$old_seed)) {
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  } else {
    assign(".Random.seed", state$old_seed, envir = .GlobalEnv)
  }
  invisible(NULL)
}
