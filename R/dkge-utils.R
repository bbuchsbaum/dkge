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
#' All public scalar validators report the argument, the supplied shape/value,
#' and the expected domain before any coercion can truncate or recycle input.
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
  valid <- is.numeric(x) && length(x) == 1L && is.finite(x) &&
    isTRUE(predicate(x))
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
.dkge_validate_finite_scalar <- function(x, arg) {
  .dkge_validate_scalar_domain(
    x, arg, "a finite numeric scalar", function(value) TRUE
  )
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
.dkge_validate_probability <- function(x, arg) {
  .dkge_validate_scalar_domain(
    x, arg, "a probability in the closed interval [0, 1]",
    function(value) value >= 0 && value <= 1
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

#' Central spectral tolerance contract
#'
#' Numerical rank decisions use one scale-equivariant absolute-plus-relative
#' threshold. The absolute term is deliberately the smallest positive normal
#' double: it protects the exactly-zero case without changing rank merely
#' because an otherwise representable matrix is rescaled.
#'
#' @param values Numeric spectrum.
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
    .dkge_abort("A finite numeric spectrum is required for rank selection.",
                "dkge_spectral_error")
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
    rank = sum(positive)
  )
}

#' PSD square root and Moore-Penrose inverse square root
#'
#' Null directions remain exactly null. No ridge or jitter is introduced.
#'
#' @param K Validated positive-semidefinite kernel.
#' @inheritParams .dkge_spectral_contract
#' @return Kernel roots and the applied range-space diagnostics.
#' @keywords internal
#' @noRd
.dkge_psd_roots <- function(
    K,
    absolute_tolerance = .Machine$double.xmin,
    relative_tolerance = 1e-8) {
  K <- .dkge_validate_kernel(K)
  ee <- eigen(K, symmetric = TRUE)
  raw_values <- ee$values
  values <- pmax(raw_values, 0)
  spectral <- .dkge_spectral_contract(
    values,
    absolute_tolerance = absolute_tolerance,
    relative_tolerance = relative_tolerance
  )
  retained <- spectral$positive
  sqrt_values <- numeric(length(values))
  inverse_sqrt_values <- numeric(length(values))
  sqrt_values[retained] <- sqrt(values[retained])
  inverse_sqrt_values[retained] <- 1 / sqrt_values[retained]
  V <- ee$vectors
  n <- length(values)
  Khalf <- V %*% diag(sqrt_values, n) %*% t(V)
  Kihalf <- V %*% diag(inverse_sqrt_values, n) %*% t(V)
  dimnames(Khalf) <- dimnames(K)
  dimnames(Kihalf) <- dimnames(K)
  list(
    Khalf = Khalf,
    Kihalf = Kihalf,
    evals = values,
    raw_evals = raw_values,
    evecs = V,
    rank = spectral$rank,
    nullity = n - spectral$rank,
    retained = retained,
    absolute_tolerance = spectral$absolute_tolerance,
    relative_tolerance = spectral$relative_tolerance,
    tolerance = spectral$tolerance,
    spectral_scale = spectral$scale,
    n_clamped = sum(!retained)
  )
}

#' Select and verify a K-orthonormal basis from a transformed moment
#'
#' @param moment Symmetric transformed q-space moment.
#' @param K Stored design kernel.
#' @param roots Output from [.dkge_psd_roots()].
#' @param requested_rank Requested component count.
#' @param vectors Optional precomputed orthogonal candidate vectors.
#' @param values Optional associated component values.
#' @param label Diagnostic label used in errors.
#' @inheritParams .dkge_spectral_contract
#' @return Basis, spectrum, ranks, and tolerance diagnostics.
#' @keywords internal
#' @noRd
.dkge_select_k_basis <- function(
    moment,
    K,
    roots,
    requested_rank,
    vectors = NULL,
    values = NULL,
    label = "transformed moment",
    absolute_tolerance = roots$absolute_tolerance,
    relative_tolerance = roots$relative_tolerance) {
  moment <- (moment + t(moment)) / 2
  if (is.null(vectors) || is.null(values)) {
    decomposition <- eigen(moment, symmetric = TRUE)
    vectors <- decomposition$vectors
    values <- decomposition$values
  } else {
    decomposition <- list(vectors = vectors, values = values)
  }
  spectral <- .dkge_spectral_contract(
    values,
    absolute_tolerance = absolute_tolerance,
    relative_tolerance = relative_tolerance
  )
  moment_rank <- spectral$rank
  effective_rank <- min(roots$rank, moment_rank, ncol(vectors))
  rank <- min(as.integer(requested_rank), effective_rank)

  if (rank == 0L) {
    U_hat <- matrix(0, nrow = nrow(K), ncol = 0L)
    U <- matrix(0, nrow = nrow(K), ncol = 0L)
    kept_values <- numeric(0)
  } else {
    keep <- which(spectral$positive)[seq_len(rank)]
    U_hat_candidate <- vectors[, keep, drop = FALSE]
    kept_values <- values[keep]
    U <- roots$Kihalf %*% U_hat_candidate
    gram <- crossprod(U, K %*% U)
    gram_error <- max(abs(gram - diag(rank)))
    if (!is.finite(gram_error) || gram_error > 1e-7) {
      .dkge_abort(
        sprintf(
          "%s basis failed the K-orthonormal postcondition (maximum Gram error %.3e).",
          label, gram_error
        ),
        "dkge_kernel_geometry_error"
      )
    }
    U_hat <- roots$Khalf %*% U
  }

  list(
    decomposition = decomposition,
    vectors = vectors,
    values = values,
    kept_values = kept_values,
    U_hat = U_hat,
    U = U,
    rank = rank,
    moment_rank = moment_rank,
    effective_rank = effective_rank,
    rank_reduced = rank < requested_rank,
    moment_tolerance = spectral$tolerance,
    moment_scale = spectral$scale
  )
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
