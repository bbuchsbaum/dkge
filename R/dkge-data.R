
# dkge-data.R
# Subject-level constructors, data bundling, and high-level DKGE entry point.

#' Generate default effect labels
#'
#' @param q Number of design effects.
#' @return Character vector of effect labels (`effect1`, `effect2`, ...).
#' @keywords internal
#' @noRd
.default_effect_names <- function(q) paste0("effect", seq_len(q))

.dkge_check_unique_effects <- function(effects, where = "Effect labels") {
  effects <- as.character(effects)
  if (anyDuplicated(effects)) {
    dup <- unique(effects[duplicated(effects)])
    stop(sprintf("%s must be unique; duplicated: %s.",
                 where, paste(dup, collapse = ", ")),
         call. = FALSE)
  }
  invisible(effects)
}

#' Harmonise beta rows with design-effect names
#'
#' Whichever side carries labels supplies them: design column names win when
#' present, otherwise the beta row names are adopted, and placeholders
#' (`effect1`, ...) are invented only when neither side is labelled. Inventing
#' placeholders while the other side is labelled would manufacture a name
#' clash out of two consistent inputs.
#'
#' @param beta qxP matrix of subject coefficients.
#' @param design Txq design matrix with named columns.
#' @return List containing aligned `beta`, `design`, and the shared effect labels.
#' @keywords internal
#' @noRd
.align_effects <- function(beta, design) {
  stopifnot(is.matrix(beta), is.matrix(design))
  if (nrow(beta) != ncol(design)) {
    .dkge_abort(
      sprintf(
        paste0(
          "Beta matrix effect-row count (%d) must match design column ",
          "count (%d)."
        ),
        nrow(beta), ncol(design)
      ),
      "dkge_dimension_error"
    )
  }
  effects <- colnames(design)
  beta_names <- rownames(beta)
  if (is.null(effects)) {
    usable_beta_names <- !is.null(beta_names) && all(nzchar(beta_names)) &&
      !anyDuplicated(beta_names)
    effects <- if (usable_beta_names) {
      as.character(beta_names)
    } else {
      .default_effect_names(ncol(design))
    }
    colnames(design) <- effects
  } else {
    effects <- as.character(effects)
  }
  .dkge_check_unique_effects(effects, "Effect labels")
  if (is.null(beta_names)) {
    rownames(beta) <- effects
  }
  if (!setequal(rownames(beta), effects)) {
    .dkge_abort(
      "Row names of beta matrix must match design column names (effects).",
      "dkge_effect_alignment_error"
    )
  }
  if (!identical(rownames(beta), effects)) {
    beta <- beta[effects, , drop = FALSE]
  }
  list(beta = beta, design = design, effects = effects)
}

#' Validate optional spatial weights
#'
#' @param omega Weight specification (`NULL`, vector, or matrix).
#' @param clusters Number of spatial units (P).
#' @return Normalised weight object suitable for multiplication.
#' @keywords internal
#' @noRd
.validate_omega <- function(omega, clusters) {
  if (is.null(omega) || length(omega) == 0) {
    return(NULL)
  }
  if (is.vector(omega)) {
    omega <- as.numeric(omega)
    if (length(omega) != clusters) {
      stop("omega vector length must match the number of clusters.", call. = FALSE)
    }
    if (any(!is.finite(omega))) {
      stop("omega vector must contain only finite values.", call. = FALSE)
    }
    if (any(omega < 0)) {
      stop("omega vector entries must be non-negative.", call. = FALSE)
    }
    return(omega)
  }
  if (is.matrix(omega)) {
    if (nrow(omega) != clusters || ncol(omega) != clusters) {
      stop("omega matrix must have dimensions clusters x clusters.", call. = FALSE)
    }
    if (any(!is.finite(omega))) {
      stop("omega matrix must contain only finite values.", call. = FALSE)
    }
    asym <- max(abs(omega - t(omega)))
    scale <- max(1, max(abs(omega)))
    if (asym > 1e-8 * scale) {
      stop("omega matrix must be symmetric.", call. = FALSE)
    }
    omega <- (omega + t(omega)) / 2
    eig <- eigen(omega, symmetric = TRUE, only.values = TRUE)$values
    eig_scale <- max(1, max(abs(eig)))
    tol <- 1e-8 * eig_scale
    if (min(eig) < -tol) {
      stop("omega matrix must be positive semidefinite.", call. = FALSE)
    }
    return(omega)
  }
  stop("Unsupported omega type; supply NULL, a numeric vector, or a matrix.", call. = FALSE)
}

# -------------------------------------------------------------------------
# Subject constructors ----------------------------------------------------
# -------------------------------------------------------------------------

#' Construct a DKGE subject record
#'
#' @param x Source object containing subject-level data:
#'   - matrix: qxP beta coefficients (effects x clusters/voxels)
#'   - NeuroVec: 4D time-series data (TxXxYxZ), betas computed via GLM
#'   - ClusteredNeuroVec: Cluster time-series (TxK), betas computed via GLM
#' @param ... Additional arguments passed to methods. For the matrix method:
#'   `design` (Subject design matrix T_s x q), `id` (Optional subject identifier),
#'   `omega` (Optional cluster weights - numeric vector length P or PxP matrix),
#'   `observed_rows` (Optional observed effect-row indices in a global effect
#'   space; defaults to all local rows), `effect_n` (per-effect trial counts),
#'   `effect_precision` (direct per-effect precision), `effect_noise_cov`
#'   (q-by-q covariance multiplier such as `(X'X)^-1`), `residual_variance`
#'   (per-column residual variances), `noise_trace` (optional precomputed spatial
#'   noise trace), and `split_betas` (two q-by-P half estimates).
#'   For ClusteredNeuroVec method: omega defaults to cluster sizes if not provided
#' @return Object of class `dkge_subject`
#' @export
#' @examples
#' betas <- matrix(rnorm(5 * 200), 5, 200)
#' design <- matrix(rnorm(150 * 5), 150, 5, dimnames = list(NULL, paste0("eff", 1:5)))
#' subj <- dkge_subject(betas, design, id = "sub01")
#' str(subj)
dkge_subject <- function(x, ...) {
  UseMethod("dkge_subject")
}

#' @export
dkge_subject.dkge_subject <- function(x, ...) {
  x
}

.validate_observed_rows <- function(observed_rows, q) {
  if (is.null(observed_rows)) {
    return(seq_len(q))
  }
  if (is.logical(observed_rows)) {
    if (length(observed_rows) != q) {
      stop("logical `observed_rows` must have length equal to the number of effects.",
           call. = FALSE)
    }
    observed_rows <- which(observed_rows)
    if (!length(observed_rows)) {
      stop("`observed_rows` must select at least one effect row.", call. = FALSE)
    }
    return(observed_rows)
  }
  observed_rows <- as.integer(observed_rows)
  if (!length(observed_rows)) {
    stop("`observed_rows` must select at least one effect row.", call. = FALSE)
  }
  if (anyNA(observed_rows) || any(observed_rows < 1L) || any(observed_rows > q)) {
    stop("`observed_rows` must contain valid 1-based effect-row indices.",
         call. = FALSE)
  }
  sort(unique(observed_rows))
}

.validate_effect_vector <- function(x, effects, name,
                                    nonnegative = TRUE,
                                    allow_na = FALSE,
                                    unit = "design effect") {
  if (is.null(x)) return(NULL)
  x_names <- names(x)
  x <- as.numeric(x)
  if (length(x) != length(effects)) {
    stop(sprintf("`%s` must have one entry per %s (got %d, expected %d).",
                 name, unit, length(x), length(effects)),
         call. = FALSE)
  }
  if (!is.null(x_names)) {
    idx <- match(effects, x_names)
    if (anyNA(idx)) {
      stop(sprintf("Named `%s` must contain every %s.", name, unit),
           call. = FALSE)
    }
    x <- x[idx]
  }
  if ((!allow_na && anyNA(x)) || any(!is.finite(x[!is.na(x)]))) {
    stop(sprintf("`%s` must contain only finite values%s.", name,
                 if (allow_na) " or NA" else ""), call. = FALSE)
  }
  if (nonnegative && any(x < 0, na.rm = TRUE)) {
    stop(sprintf("`%s` must be non-negative.", name), call. = FALSE)
  }
  names(x) <- effects
  x
}

.validate_effect_covariance <- function(x, effects, name = "effect_noise_cov") {
  if (is.null(x)) return(NULL)
  x <- as.matrix(x)
  q <- length(effects)
  if (!identical(dim(x), c(q, q))) {
    stop(sprintf("`%s` must be a q x q matrix.", name), call. = FALSE)
  }
  if (any(!is.finite(x))) {
    stop(sprintf("`%s` must contain only finite values.", name), call. = FALSE)
  }
  if (!is.null(rownames(x)) || !is.null(colnames(x))) {
    if (is.null(rownames(x)) || is.null(colnames(x)) ||
        !setequal(rownames(x), effects) || !setequal(colnames(x), effects)) {
      stop(sprintf("Dimnames of `%s` must match the design effects.", name),
           call. = FALSE)
    }
    x <- x[effects, effects, drop = FALSE]
  }
  asym <- max(abs(x - t(x)))
  scale <- max(1, max(abs(x)))
  if (asym > 1e-8 * scale) {
    stop(sprintf("`%s` must be symmetric.", name), call. = FALSE)
  }
  x <- (x + t(x)) / 2
  ev <- eigen(x, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) < -1e-8 * max(1, max(abs(ev)))) {
    stop(sprintf("`%s` must be positive semidefinite.", name), call. = FALSE)
  }
  dimnames(x) <- list(effects, effects)
  x
}

.validate_split_betas <- function(x, effects, P, cluster_ids = NULL) {
  if (is.null(x)) return(NULL)
  if (!is.list(x) || length(x) != 2L) {
    stop("`split_betas` must be a two-element list of q x P matrices.",
         call. = FALSE)
  }
  out <- lapply(x, function(B) {
    B <- as.matrix(B)
    if (!identical(dim(B), c(length(effects), P)) || any(!is.finite(B))) {
      stop("Each `split_betas` entry must be a finite q x P matrix.",
           call. = FALSE)
    }
    if (!is.null(rownames(B))) {
      idx <- match(effects, rownames(B))
      if (anyNA(idx)) stop("Split-beta row names must match design effects.", call. = FALSE)
      B <- B[idx, , drop = FALSE]
    }
    rownames(B) <- effects
    if (!is.null(cluster_ids)) colnames(B) <- cluster_ids
    B
  })
  names(out) <- names(x) %||% c("first", "second")
  out
}

#' @export
dkge_subject.matrix <- function(x, design, id = NULL, omega = NULL,
                                observed_rows = NULL,
                                effect_n = NULL,
                                effect_precision = NULL,
                                effect_noise_cov = NULL,
                                residual_variance = NULL,
                                residual_df = NULL,
                                noise_trace = NULL,
                                split_betas = NULL,
                                ...) {
  stopifnot(is.matrix(x), is.matrix(design))
  aligned <- .align_effects(x, design)
  beta <- aligned$beta
  design <- aligned$design
  effects <- aligned$effects
  P <- ncol(beta)
  omega <- .validate_omega(omega, P)
  observed_rows <- .validate_observed_rows(observed_rows, length(effects))
  effect_n <- .validate_effect_vector(effect_n, effects, "effect_n")
  effect_precision <- .validate_effect_vector(effect_precision, effects,
                                              "effect_precision")
  effect_noise_cov <- .validate_effect_covariance(effect_noise_cov, effects)
  split_betas <- .validate_split_betas(split_betas, effects, P, colnames(beta))
  residual_variance <- .validate_effect_vector(
    residual_variance,
    colnames(beta) %||% paste0("cluster_", seq_len(P)),
    "residual_variance",
    allow_na = TRUE,
    unit = "cluster/voxel column of `x`"
  )
  if (!is.null(residual_df)) {
    residual_df <- as.numeric(residual_df)
    if (length(residual_df) != 1L || !is.finite(residual_df) || residual_df < 0) {
      stop("`residual_df` must be a finite non-negative scalar.", call. = FALSE)
    }
  }
  if (!is.null(noise_trace)) {
    noise_trace <- as.numeric(noise_trace)
    if (length(noise_trace) != 1L || !is.finite(noise_trace) || noise_trace < 0) {
      stop("`noise_trace` must be a finite non-negative scalar.", call. = FALSE)
    }
  }

  # Check for rank deficiency (warnings only, does not block construction)
  .dkge_check_rank(design, beta, subject_id = id)
  if (is.null(colnames(beta))) {
    colnames(beta) <- paste0("cluster_", seq_len(P))
  }
  structure(list(
    id = if (is.null(id)) NA_character_ else as.character(id),
    beta = beta,
    design = design,
    omega = omega,
    effects = effects,
    observed_rows = observed_rows,
    effect_n = effect_n,
    effect_precision = effect_precision,
    effect_noise_cov = effect_noise_cov,
    residual_variance = residual_variance,
    residual_df = residual_df,
    noise_trace = noise_trace,
    split_betas = split_betas,
    n_clusters = P,
    cluster_ids = colnames(beta)
  ), class = "dkge_subject")
}

#' @export
dkge_subject.list <- function(x, ...) {
  stopifnot(!is.null(x$beta), !is.null(x$design))
  dkge_subject(as.matrix(x$beta), design = as.matrix(x$design), id = x$id,
               omega = x$omega, observed_rows = x$observed_rows,
               effect_n = x$effect_n,
               effect_precision = x$effect_precision,
               effect_noise_cov = x$effect_noise_cov,
               residual_variance = x$residual_variance,
               residual_df = x$residual_df,
               noise_trace = x$noise_trace,
               split_betas = x$split_betas,
               ...)
}

#' @export
dkge_subject.default <- function(x, ...) {
  stop("Unsupported object type for dkge_subject().", call. = FALSE)
}

#' @export
dkge_subject.NeuroVec <- function(x, design, id = NULL, omega = NULL, mask = NULL, compute_betas = TRUE, ...) {
  if (!requireNamespace("neuroim2", quietly = TRUE)) {
    stop("Install 'neuroim2' to build dkge subjects from NeuroVec objects.", call. = FALSE)
  }

  # Check if x is a 4D time-series or already contains betas
  dims <- dim(x)

  if (length(dims) == 4 && compute_betas) {
    # This is a 4D time-series, need to compute betas
    if (!requireNamespace("fmrireg", quietly = TRUE)) {
      stop("Install 'fmrireg' to compute GLM betas from NeuroVec time-series.", call. = FALSE)
    }

    # Convert to matrix (voxels x time)
    mat <- neuroim2::as.matrix(x)

    # Apply mask if provided
    if (!is.null(mask)) {
      if (inherits(mask, "LogicalNeuroVol")) {
        mask_vals <- neuroim2::values(mask)
        idx <- which(as.logical(mask_vals))
      } else if (is.numeric(mask)) {
        idx <- as.integer(mask)
      } else {
        stop("mask must be a LogicalNeuroVol or numeric indices.", call. = FALSE)
      }
      if (!length(idx)) stop("Mask selects zero voxels.", call. = FALSE)
      mat <- mat[idx, , drop = FALSE]
    }

    # Transpose to get time x voxels
    mat <- t(mat)

    # Compute betas using GLM
    fit <- fmrireg::fmri_ols_fit(design, mat)
    beta <- fit$betas  # q x voxels

  } else {
    # Assume x already contains betas
    mat <- neuroim2::as.matrix(x)  # voxels x effects/time

    if (!is.null(mask)) {
      if (inherits(mask, "LogicalNeuroVol")) {
        mask_vals <- neuroim2::values(mask)
        idx <- which(as.logical(mask_vals))
      } else if (is.numeric(mask)) {
        idx <- as.integer(mask)
      } else {
        stop("mask must be a LogicalNeuroVol or numeric indices.", call. = FALSE)
      }
      if (!length(idx)) stop("Mask selects zero voxels.", call. = FALSE)
      mat <- mat[idx, , drop = FALSE]
    }

    beta <- t(mat)  # effects x voxels
  }

  if (is.null(colnames(beta))) {
    colnames(beta) <- paste0("voxel_", seq_len(ncol(beta)))
  }

  dkge_subject.matrix(beta, design = design, id = id, omega = omega, ...)
}

#' @export
dkge_subject.ClusteredNeuroVec <- function(x, design, id = NULL, omega = NULL, ...) {
  if (!requireNamespace("neuroim2", quietly = TRUE)) {
    stop("Install 'neuroim2' to build dkge subjects from ClusteredNeuroVec objects.", call. = FALSE)
  }

  # Extract the TxK cluster time-series matrix using neuroim2 method
  cluster_ts <- neuroim2::as.matrix(x)  # This returns TxK matrix with cluster labels as colnames

  # For GLM fitting, we need to compute betas from time-series
  # This requires the design matrix
  if (!requireNamespace("fmrireg", quietly = TRUE)) {
    stop("Install 'fmrireg' to compute GLM betas from ClusteredNeuroVec.", call. = FALSE)
  }

  # Fit GLM to get betas (qxK)
  fit <- fmrireg::fmri_ols_fit(design, cluster_ts)
  beta <- fit$betas  # qxK matrix

  # Use cluster labels from the ClusteredNeuroVec if available
  if (is.null(colnames(beta)) || all(grepl("^Cluster_", colnames(beta)))) {
    # If colnames are generic or missing, try to get meaningful labels
    if (!is.null(x@cvol@label_map) && length(x@cvol@label_map) == ncol(beta)) {
      colnames(beta) <- names(x@cvol@label_map)
    } else if (is.null(colnames(beta))) {
      colnames(beta) <- paste0("cluster_", seq_len(ncol(beta)))
    }
  }

  # If omega not provided, could use cluster sizes as weights
  if (is.null(omega) && !is.null(x@cvol)) {
    # Get cluster sizes as potential weights
    cluster_sizes <- table(x@cl_map[x@cl_map > 0])
    if (length(cluster_sizes) == ncol(beta)) {
      omega <- as.numeric(cluster_sizes)
    }
  }

  dkge_subject.matrix(beta, design = design, id = id, omega = omega, ...)
}

# -------------------------------------------------------------------------
# Subject bundle ---------------------------------------------------------
# -------------------------------------------------------------------------

#' Check whether an object is a dkge_subject
#'
#' @param x Object to test.
#' @return Logical scalar indicating whether `x` inherits `dkge_subject`.
#' @keywords internal
#' @noRd
.is_dkge_subject <- function(x) inherits(x, "dkge_subject")

#' Normalise subject identifiers
#'
#' @param subjects List of `dkge_subject` objects.
#' @param provided Optional identifiers supplied by the caller.
#' @return Character vector of subject IDs applied to the subjects.
#' @keywords internal
#' @noRd
.normalize_subject_ids <- function(subjects, provided = NULL) {
  n <- length(subjects)
  ids <- vapply(subjects, function(s) s$id, character(1))
  if (!is.null(provided)) {
    stopifnot(length(provided) == n)
    ids <- as.character(provided)
  }
  missing_idx <- which(is.na(ids) | ids == "")
  if (length(missing_idx)) {
    fmt <- paste0("%0", max(2L, nchar(n)), "d")
    defaults <- paste0("sub", sprintf(fmt, seq_len(n)))
    ids[missing_idx] <- defaults[missing_idx]
  }
  # Duplicate IDs would let name-keyed lookups (observation masks, between-
  # subject designs) silently reuse one subject's metadata for another.
  if (anyDuplicated(ids)) {
    dupes <- unique(ids[duplicated(ids)])
    stop(sprintf(
      "Subject identifiers must be unique; duplicated: %s.",
      paste(dupes, collapse = ", ")
    ), call. = FALSE)
  }
  ids
}

#' Permute per-effect provenance to a new effect order
#'
#' @param provenance Provenance list keyed by the current effect order.
#' @param perm Integer vector such that `old_effects[perm]` is the new order.
#' @param effects New effect labels.
#' @keywords internal
#' @noRd
.dkge_reorder_provenance <- function(provenance, perm, effects) {
  if (is.null(provenance)) return(NULL)
  provenance$effect_ids <- effects
  if (!is.null(provenance$obs_mask)) {
    provenance$obs_mask <- lapply(provenance$obs_mask, function(mask) {
      mask <- as.logical(mask)[perm]
      names(mask) <- effects
      mask
    })
    provenance$observed_rows <- lapply(provenance$obs_mask,
                                       function(mask) unname(which(mask)))
  }
  if (!is.null(provenance$pair_counts)) {
    pc <- provenance$pair_counts[perm, perm, drop = FALSE]
    dimnames(pc) <- list(effects, effects)
    provenance$pair_counts <- pc
  }
  if (!is.null(provenance$coverage)) {
    cov <- provenance$coverage[perm, , drop = FALSE]
    cov$effect <- effects
    rownames(cov) <- NULL
    provenance$coverage <- cov
  }
  provenance
}

#' Bundle subject-level inputs for DKGE
#'
#' @param betas List of subject records (matrices or `dkge_subject` objects)
#' @param designs Optional list of design matrices (ignored when `betas` already contain subjects)
#' @param omega Optional list of cluster weights
#' @param subject_ids Optional subject identifiers
#' @param effects Optional character vector pinning the global effect order.
#'   It must be a permutation of the effects observed across subjects (the
#'   union when coverage is partial). Supply
#'   `dkge_effect_grid()$cell_labels` here so that the bundle, the design
#'   kernel, and the grid all index effects identically; without it the union
#'   is ordered by first appearance across subjects, which depends on the order
#'   of `betas`.
#' @return An object of class `dkge_data`
#' @export
#' @examples
#' betas <- replicate(3, matrix(rnorm(5 * 80), 5, 80), simplify = FALSE)
#' designs <- replicate(
#'   3,
#'   matrix(
#'     rnorm(150 * 5), 150, 5,
#'     dimnames = list(NULL, paste0("eff", 1:5))
#'   ),
#'   simplify = FALSE
#' )
#' data <- dkge_data(betas, designs)
#' data$effects
#'
#' # Pin a specific effect order
#' dkge_data(betas, designs, effects = paste0("eff", 5:1))$effects
dkge_data <- function(betas, designs = NULL, omega = NULL, subject_ids = NULL,
                      effects = NULL) {
  if (inherits(betas, "dkge_subject")) {
    subjects <- list(dkge_subject(betas))
  } else if (is.list(betas) && length(betas) > 0 && all(vapply(betas, .is_dkge_subject, logical(1)))) {
    subjects <- lapply(betas, dkge_subject)
  } else {
    stopifnot(is.list(betas), is.list(designs), length(betas) == length(designs), length(betas) > 0)
    n <- length(betas)
    if (is.null(omega)) {
      omega <- vector("list", n)
    } else if (!is.list(omega)) {
      omega <- as.list(omega)
    }
    stopifnot(length(omega) == n)
    subjects <- vector("list", n)
    for (i in seq_len(n)) {
      id <- if (is.null(subject_ids)) NULL else subject_ids[[i]]
      subjects[[i]] <- dkge_subject(betas[[i]], design = designs[[i]], id = id, omega = omega[[i]])
    }
  }

  # Require minimum 2 subjects for group analysis
  if (length(subjects) < 2) {
    stop("At least 2 subjects required for group analysis.", call. = FALSE)
  }

  ids <- .normalize_subject_ids(subjects, provided = subject_ids)
  for (i in seq_along(subjects)) subjects[[i]]$id <- ids[[i]]
  effects_list <- lapply(subjects, `[[`, "effects")
  for (i in seq_along(effects_list)) {
    .dkge_check_unique_effects(
      effects_list[[i]],
      sprintf("Subject '%s' effect labels", subjects[[i]]$id %||% "(unnamed)")
    )
  }
  reference_effects <- as.character(effects_list[[1]])
  identical_effects <- all(vapply(effects_list, function(e) identical(as.character(e), reference_effects), logical(1)))

  provenance <- NULL
  if (!identical_effects) {
    aligned <- .dkge_align_subjects_to_union(subjects)
    subjects <- aligned$subjects
    provenance <- aligned$provenance
    reference_effects <- provenance$effect_ids
  } else {
    provenance <- .dkge_full_coverage_provenance(subjects)
  }

  effects_ref <- as.character(reference_effects)

  if (!is.null(effects)) {
    effects <- as.character(effects)
    if (anyDuplicated(effects)) {
      stop("`effects` must not contain duplicates.", call. = FALSE)
    }
    if (!setequal(effects, effects_ref)) {
      missing_eff <- setdiff(effects_ref, effects)
      extra_eff <- setdiff(effects, effects_ref)
      stop(sprintf(
        "`effects` must be a permutation of the effects observed across subjects.%s%s",
        if (length(missing_eff)) sprintf(" Missing: %s.", paste(missing_eff, collapse = ", ")) else "",
        if (length(extra_eff)) sprintf(" Unknown: %s.", paste(extra_eff, collapse = ", ")) else ""
      ), call. = FALSE)
    }
    perm <- match(effects, effects_ref)
    provenance <- .dkge_reorder_provenance(provenance, perm, effects)
    effects_ref <- effects
  }

  for (i in seq_along(subjects)) {
    subj <- subjects[[i]]
    perm <- match(effects_ref, as.character(subj$effects))
    if (anyNA(perm)) {
      stop(sprintf("Subject '%s' effects do not match the reference effects.",
                   subj$id), call. = FALSE)
    }
    design_effects <- colnames(subj$design)
    if (!identical(design_effects, effects_ref)) {
      match_idx <- match(effects_ref, design_effects)
      if (anyNA(match_idx)) {
        stop(sprintf("Subject '%s' design columns do not match reference effects.", subjects[[i]]$id), call. = FALSE)
      }
      subj$design <- subj$design[, match_idx, drop = FALSE]
    }
    beta_effects <- rownames(subj$beta)
    if (!identical(beta_effects, effects_ref)) {
      match_idx <- match(effects_ref, beta_effects)
      if (anyNA(match_idx)) {
        stop(sprintf("Subject '%s' beta rows do not match reference effects.", subjects[[i]]$id), call. = FALSE)
      }
      subj$beta <- subj$beta[match_idx, , drop = FALSE]
    }
    if (!identical(perm, seq_along(effects_ref))) {
      # Per-effect metadata must follow the beta rows, otherwise reliability
      # weighting and analytic debiasing index the wrong effects.
      if (!is.null(subj$effect_n)) subj$effect_n <- subj$effect_n[perm]
      if (!is.null(subj$effect_precision)) {
        subj$effect_precision <- subj$effect_precision[perm]
      }
      if (!is.null(subj$effect_noise_cov)) {
        subj$effect_noise_cov <- subj$effect_noise_cov[perm, perm, drop = FALSE]
      }
      if (!is.null(subj$split_betas)) {
        subj$split_betas <- lapply(subj$split_betas,
                                   function(B) B[perm, , drop = FALSE])
      }
    }
    if (!is.null(subj$effect_n)) names(subj$effect_n) <- effects_ref
    if (!is.null(subj$effect_precision)) names(subj$effect_precision) <- effects_ref
    if (!is.null(subj$effect_noise_cov)) {
      dimnames(subj$effect_noise_cov) <- list(effects_ref, effects_ref)
    }
    if (!is.null(subj$split_betas)) {
      subj$split_betas <- lapply(subj$split_betas, function(B) {
        rownames(B) <- effects_ref
        B
      })
    }
    rownames(subj$beta) <- effects_ref
    colnames(subj$design) <- effects_ref
    subj$effects <- effects_ref
    mask <- as.logical(provenance$obs_mask[[i]] %||%
                         rep(TRUE, length(effects_ref)))
    subj$observed_rows <- unname(which(mask))
    if (any(!mask)) {
      subj$beta[!mask, ] <- 0
      subj$design[, !mask] <- 0
      if (!is.null(subj$split_betas)) {
        subj$split_betas <- lapply(subj$split_betas, function(B) {
          B[!mask, ] <- 0
          B
        })
      }
    }
    subjects[[i]] <- subj
  }

  # An effect that nobody observes contributes an identically zero row/column
  # to every moment; downstream eigensolves then report a spurious null
  # direction instead of a data problem.
  observed_any <- Reduce(`|`, lapply(subjects, function(s) {
    m <- rep(FALSE, length(effects_ref))
    m[s$observed_rows] <- TRUE
    m
  }))
  if (!all(observed_any)) {
    stop(sprintf(
      "Effect(s) observed by no subject: %s. Drop them from the design or supply data for them.",
      paste(effects_ref[!observed_any], collapse = ", ")
    ), call. = FALSE)
  }

  structure(list(
    subjects = subjects,
    betas = lapply(subjects, `[[`, "beta"),
    designs = lapply(subjects, `[[`, "design"),
    omega = lapply(subjects, `[[`, "omega"),
    subject_ids = ids,
    effects = effects_ref,
    q = length(effects_ref),
    n_subjects = length(subjects),
    cluster_ids = lapply(subjects, `[[`, "cluster_ids"),
    observed_rows = lapply(subjects, `[[`, "observed_rows"),
    effect_n = lapply(subjects, `[[`, "effect_n"),
    effect_precision = lapply(subjects, `[[`, "effect_precision"),
    effect_noise_cov = lapply(subjects, `[[`, "effect_noise_cov"),
    residual_variance = lapply(subjects, `[[`, "residual_variance"),
    residual_df = lapply(subjects, `[[`, "residual_df"),
    noise_trace = lapply(subjects, `[[`, "noise_trace"),
    split_betas = lapply(subjects, `[[`, "split_betas"),
    provenance = provenance
  ), class = "dkge_data")
}

# -------------------------------------------------------------------------
# High-level entry point --------------------------------------------------
# -------------------------------------------------------------------------

#' Fit DKGE across multiple subjects
#'
#' This high-level wrapper prepares subject-level inputs and calls [dkge_fit()]
#' to estimate the shared Design-Kernel Group Embedding (DKGE) basis. Provide a
#' design kernel in effect space together with per-subject GLM coefficients and
#' design matrices, and receive a `dkge` object ready for LOSO contrasts, medoid
#' transport, or out-of-sample prediction.
#'
#' @param betas Either (i) a list of qxP_s beta matrices (one per subject), (ii)
#'   a list/tibble of [dkge_subject()] objects, or (iii) a pre-built `dkge_data`
#'   bundle. Rows must align with design effects.
#' @param designs Optional list of T_sxq design matrices from the subject GLMs.
#'   Each column corresponds to an effect/regressor; column names become the
#'   canonical effect labels enforced across subjects. Ignored when `betas`
#'   already carries `dkge_subject`/`dkge_data` entries.
#' @param K Design kernel that expresses similarity or smoothness between
#'   effects in the design space (e.g. identity for nominal factors, RBF for
#'   ordinal factors, Kronecker combinations for interactions). Supply either a
#'   qxq numeric matrix or the list returned by [design_kernel()], whose `K`
#'   element is extracted automatically. Matches the `K` parameter name used by
#'   [dkge_fit()].
#' @param Omega_list Optional list overriding per-subject spatial weights. Each
#'   element may be `NULL`, a length-P_s numeric vector, or a P_sxP_s matrix.
#'   When betas are voxelwise (e.g. coming from `NeuroVec`), these weights
#'   operate on spatial units (voxels) rather than clusters; equal weights are
#'   assumed when omitted. Matches the `Omega_list` parameter name used by
#'   [dkge_fit()].
#' @param kernel Deprecated. Use `K` instead.
#' @param omega Deprecated. Use `Omega_list` instead.
#' @param subject_ids Optional subject identifiers used when raw matrices are
#'   provided. Ignored when `betas` carries IDs already.
#' @param keep_inputs Logical; when `TRUE` (default) the returned object stores the
#'   canonicalised `dkge_data` bundle under `$input` for later inspection or
#'   prediction.
#' @param cpca_blocks Optional integer vector specifying the effect rows that
#'   span a CPCA design subspace. Ignored when `cpca_part = "none"` or when
#'   `cpca_T` is provided.
#' @param cpca_T Optional qxq0 matrix giving the CPCA design basis explicitly.
#'   Overrides `cpca_blocks` when supplied.
#' @param cpca_part Which CPCA-filtered component to fit: `"none"` (default)
#'   performs the standard DKGE fit; `"design"` uses only the CPCA design part;
#'   `"resid"` uses the residual part; `"both"` fits the design part but also
#'   stores the residual basis.
#' @param cpca_ridge Optional ridge applied to the CPCA-filtered matrices before
#'   eigen-decomposition.
#' @param weights Optional [`dkge_weights()`] specification that overrides the
#'   default adaptive weighting behaviour when accumulating covariance. Use this
#'   to inject custom spatial reliabilities or to reuse a precomputed weighting
#'   recipe produced by [dkge_weights()].
#' @param effect_weights Optional [dkge_effect_weights()] specification for
#'   subject-by-effect count or precision weighting.
#' @param debias Effect-moment estimator passed to [dkge_fit()]: observed
#'   (`"none"`), analytic finite-trial subtraction (`"analytic"`), or a
#'   split-half cross-moment (`"split_half"`).
#' @param effect_scaling Effect-space scaling passed to [dkge_fit()]. Use
#'   `"none"` when input rows already share an absolute scale, such as
#'   AUC-minus-chance cell maps.
#' @param effects Optional character vector pinning the global effect order,
#'   forwarded to [dkge_data()]. Ignored (and checked for agreement) when
#'   `betas` is already a `dkge_data` bundle.
#' @param ... Additional arguments forwarded to [dkge_fit()] (e.g. `rank`,
#'   `ridge`, `w_method`, `weights`).
#'
#' @details The design matrices `X_s` supplied here are the same T_sxq regressors
#'   used to fit the subject-level GLMs; their columns must align with the effects
#'   encoded by the kernel. The kernel `K` captures how effects relate to one
#'   another-identity recovers standard OLS scaling, while structured kernels (e.g.
#'   RBF for ordinal factors, circulant for wrapped factors, Kronecker products for
#'   interactions) encourage shared smoothness or coupling between design effects.
#'
#'   DKGE itself operates entirely in this low-dimensional design space: (1) the
#'   pooled Gram matrix across subjects yields a shared Cholesky factor `R`; (2) each beta
#'   matrix is row-standardised; (3) compressed covariance is accumulated in the
#'   K-metric with optional subject and effect weighting; and (4) a tiny
#'   symmetric eigenproblem produces
#'   the K-orthonormal group basis. The input harmonisation performed by this
#'   wrapper ensures consistent effect naming, subject identifiers, and spatial
#'   weights so downstream utilities such as [dkge_loso_contrast()] can work without
#'   additional bookkeeping. Use [dkge_subject()] to build subjects from raw
#'   matrices, `NeuroVec`, or `ClusteredNeuroVec` sources prior to calling `dkge()`.
#'
#' @return A `dkge` object containing the learned basis (`$U`), eigenvalues,
#'   pooled design Cholesky factor (`$R`), compressed covariance matrix, subject weights, and
#'   metadata (subject IDs, effect names, cluster identifiers) derived from the
#'   input bundle.
#' @seealso [dkge_subject()], [dkge_data()], [design_kernel()], [dkge_fit()],
#'   [dkge_loso_contrast()]
#' @export
dkge <- function(betas, designs = NULL, K = NULL, Omega_list = NULL,
                 subject_ids = NULL,
                 keep_inputs = TRUE, cpca_blocks = NULL, cpca_T = NULL,
                 cpca_part = c("none", "design", "resid", "both"),
                 cpca_ridge = 0, weights = NULL,
                 kernel = NULL, omega = NULL,
                 effect_weights = NULL,
                 debias = c("none", "analytic", "split_half"),
                 effect_scaling = c("pooled_design", "none"),
                 effects = NULL, ...) {
  # Deprecated aliases
  if (!is.null(kernel) && is.null(K)) {
    warning("Argument 'kernel' is deprecated; use 'K' instead.", call. = FALSE)
    K <- kernel
  }
  if (!is.null(omega) && is.null(Omega_list)) {
    warning("Argument 'omega' is deprecated; use 'Omega_list' instead.", call. = FALSE)
    Omega_list <- omega
  }

  cpca_part <- match.arg(cpca_part)
  debias <- match.arg(debias)
  effect_scaling <- match.arg(effect_scaling)
  omega_override <- NULL
  if (inherits(betas, "dkge_data")) {
    data <- betas
    omega_override <- Omega_list
    if (!is.null(effects) && !identical(as.character(effects), data$effects)) {
      stop("`effects` cannot reorder an existing `dkge_data` bundle; rebuild it with dkge_data(effects = ...).",
           call. = FALSE)
    }
  } else {
    data <- dkge_data(betas, designs = designs, omega = Omega_list,
                      subject_ids = subject_ids, effects = effects)
  }

  kernel_obj <- as_dkge_kernel(K)
  K_mat <- kernel_obj$K
  stopifnot(is.matrix(K_mat), nrow(K_mat) == data$q, ncol(K_mat) == data$q)

  # Pass the kernel bundle so `.dkge_fit_prepare()` can permute `info`
  # together with `K`. Stamping the original `info` back on afterwards
  # left `kernel_info$cells` in kernel order after a data-order reorder.
  fit <- dkge_fit(data, K = kernel_obj, Omega_list = omega_override,
                     cpca_blocks = cpca_blocks, cpca_T = cpca_T,
                     cpca_part = cpca_part, cpca_ridge = cpca_ridge,
                     weights = weights, effect_weights = effect_weights,
                     debias = debias, effect_scaling = effect_scaling, ...)
  if (keep_inputs) fit$input <- data
  fit
}
