
# dkge-data.R
# Subject-level constructors, data bundling, and high-level DKGE entry point.

#' Return the first non-NULL value
#'
#' @param a Primary value to test.
#' @param b Fallback value when `a` is `NULL`.
#' @return `a` if it is not `NULL`, otherwise `b`.
#' @keywords internal
#' @noRd
`%||%` <- function(a, b) if (!is.null(a)) a else b

#' Generate default effect labels
#'
#' @param q Number of design effects.
#' @return Character vector of effect labels (`effect1`, `effect2`, ...).
#' @keywords internal
#' @noRd
.default_effect_names <- function(q) paste0("effect", seq_len(q))

#' Harmonise beta rows with design-effect names
#'
#' @param beta q×P matrix of subject coefficients.
#' @param design T×q design matrix with named columns.
#' @return List containing aligned `beta`, `design`, and the shared effect labels.
#' @keywords internal
#' @noRd
.align_effects <- function(beta, design) {
  stopifnot(is.matrix(beta), is.matrix(design))
  effects <- colnames(design)
  if (is.null(effects)) {
    effects <- .default_effect_names(ncol(design))
    colnames(design) <- effects
  } else {
    effects <- as.character(effects)
  }
  if (is.null(rownames(beta))) {
    rownames(beta) <- effects
  }
  if (!setequal(rownames(beta), effects)) {
    stop("Row names of beta matrix must match design column names (effects).", call. = FALSE)
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
    stopifnot(length(omega) == clusters)
    return(omega)
  }
  if (is.matrix(omega)) {
    stopifnot(nrow(omega) == clusters, ncol(omega) == clusters)
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
#'   - matrix: q×P beta coefficients (effects × clusters/voxels)
#'   - NeuroVec: 4D time-series data (T×X×Y×Z), betas computed via GLM
#'   - ClusteredNeuroVec: Cluster time-series (T×K), betas computed via GLM
#' @param ... Additional arguments passed to methods. For the matrix method:
#'   `design` (Subject design matrix T_s × q), `id` (Optional subject identifier),
#'   `omega` (Optional cluster weights - numeric vector length P or P×P matrix).
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

#' @export
dkge_subject.matrix <- function(x, design, id = NULL, omega = NULL, ...) {
  stopifnot(is.matrix(x), is.matrix(design))
  aligned <- .align_effects(x, design)
  beta <- aligned$beta
  design <- aligned$design
  effects <- aligned$effects
  P <- ncol(beta)
  omega <- .validate_omega(omega, P)
  if (is.null(colnames(beta))) {
    colnames(beta) <- paste0("cluster_", seq_len(P))
  }
  structure(list(
    id = if (is.null(id)) NA_character_ else as.character(id),
    beta = beta,
    design = design,
    omega = omega,
    effects = effects,
    n_clusters = P,
    cluster_ids = colnames(beta)
  ), class = "dkge_subject")
}

#' @export
dkge_subject.list <- function(x, ...) {
  stopifnot(!is.null(x$beta), !is.null(x$design))
  dkge_subject(as.matrix(x$beta), design = as.matrix(x$design), id = x$id, omega = x$omega, ...)
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

    # Convert to matrix (voxels × time)
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

    # Transpose to get time × voxels
    mat <- t(mat)

    # Compute betas using GLM
    fit <- fmrireg::fmri_ols_fit(design, mat)
    beta <- fit$betas  # q × voxels

  } else {
    # Assume x already contains betas
    mat <- neuroim2::as.matrix(x)  # voxels × effects/time

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

    beta <- t(mat)  # effects × voxels
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

  # Extract the T×K cluster time-series matrix using neuroim2 method
  cluster_ts <- neuroim2::as.matrix(x)  # This returns T×K matrix with cluster labels as colnames

  # For GLM fitting, we need to compute betas from time-series
  # This requires the design matrix
  if (!requireNamespace("fmrireg", quietly = TRUE)) {
    stop("Install 'fmrireg' to compute GLM betas from ClusteredNeuroVec.", call. = FALSE)
  }

  # Fit GLM to get betas (q×K)
  fit <- fmrireg::fmri_ols_fit(design, cluster_ts)
  beta <- fit$betas  # q×K matrix

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
  for (i in seq_len(n)) subjects[[i]]$id <- ids[i]
  ids
}

#' Bundle subject-level inputs for DKGE
#'
#' @param betas List of subject records (matrices or `dkge_subject` objects)
#' @param designs Optional list of design matrices (ignored when `betas` already contain subjects)
#' @param omega Optional list of cluster weights
#' @param subject_ids Optional subject identifiers
#' @return An object of class `dkge_data`
#' @export
#' @examples
#' betas <- replicate(3, matrix(rnorm(5 * 80), 5, 80), simplify = FALSE)
#' designs <- replicate(3, matrix(rnorm(150 * 5), 150, 5, dimnames = list(NULL, paste0("eff", 1:5))), simplify = FALSE)
#' data <- dkge_data(betas, designs)
#' data$effects
dkge_data <- function(betas, designs = NULL, omega = NULL, subject_ids = NULL) {
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

  ids <- .normalize_subject_ids(subjects, provided = subject_ids)
  effects_ref <- subjects[[1]]$effects
  for (i in seq_along(subjects)) {
    subj <- subjects[[i]]
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
    rownames(subj$beta) <- effects_ref
    colnames(subj$design) <- effects_ref
    subj$effects <- effects_ref
    subjects[[i]] <- subj
  }

  structure(list(
    betas = lapply(subjects, `[[`, "beta"),
    designs = lapply(subjects, `[[`, "design"),
    omega = lapply(subjects, `[[`, "omega"),
    subject_ids = ids,
    effects = effects_ref,
    q = length(effects_ref),
    n_subjects = length(subjects),
    cluster_ids = lapply(subjects, `[[`, "cluster_ids")
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
#' @param betas Either (i) a list of q×P_s beta matrices (one per subject), (ii)
#'   a list/tibble of [dkge_subject()] objects, or (iii) a pre-built `dkge_data`
#'   bundle. Rows must align with design effects.
#' @param designs Optional list of T_s×q design matrices from the subject GLMs.
#'   Each column corresponds to an effect/regressor; column names become the
#'   canonical effect labels enforced across subjects. Ignored when `betas`
#'   already carries `dkge_subject`/`dkge_data` entries.
#' @param kernel Design kernel that expresses similarity or smoothness between
#'   effects in the design space (e.g. identity for nominal factors, RBF for
#'   ordinal factors, Kronecker combinations for interactions). Supply either a
#'   q×q numeric matrix or the list returned by [design_kernel()], whose `K`
#'   element is extracted automatically.
#' @param omega Optional list overriding per-subject spatial weights. Each element
#'   may be `NULL`, a length-P_s numeric vector, or a P_s×P_s matrix. When betas
#'   are voxelwise (e.g. coming from `NeuroVec`), these weights operate on spatial
#'   units (voxels) rather than clusters; equal weights are assumed when omitted.
#' @param subject_ids Optional subject identifiers used when raw matrices are
#'   provided. Ignored when `betas` carries IDs already.
#' @param keep_inputs Logical; when `TRUE` (default) the returned object stores the
#'   canonicalised `dkge_data` bundle under `$input` for later inspection or
#'   prediction.
#' @param cpca_blocks Optional integer vector specifying the effect rows that
#'   span a CPCA design subspace. Ignored when `cpca_part = "none"` or when
#'   `cpca_T` is provided.
#' @param cpca_T Optional q×q0 matrix giving the CPCA design basis explicitly.
#'   Overrides `cpca_blocks` when supplied.
#' @param cpca_part Which CPCA-filtered component to fit: `"none"` (default)
#'   performs the standard DKGE fit; `"design"` uses only the CPCA design part;
#'   `"resid"` uses the residual part; `"both"` fits the design part but also
#'   stores the residual basis.
#' @param cpca_ridge Optional ridge applied to the CPCA-filtered matrices before
#'   eigen-decomposition.
#' @param ... Additional arguments forwarded to [dkge_fit()] (e.g. `rank`,
#'   `ridge`, `w_method`).
#'
#' @details The design matrices `X_s` supplied here are the same T_s×q regressors
#'   used to fit the subject-level GLMs; their columns must align with the effects
#'   encoded by the kernel. The kernel `K` captures how effects relate to one
#'   another—identity recovers standard OLS scaling, while structured kernels (e.g.
#'   RBF for ordinal factors, circulant for wrapped factors, Kronecker products for
#'   interactions) encourage shared smoothness or coupling between design effects.
#'
#'   DKGE itself operates entirely in this low-dimensional design space: (1) the
#'   pooled Gram matrix across subjects yields a shared Cholesky factor `R`; (2) each beta
#'   matrix is row-standardised; (3) compressed covariance is accumulated in the
#'   K-metric with optional subject weighting; and (4) a tiny eigenproblem produces
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
dkge <- function(betas, designs = NULL, kernel, omega = NULL, subject_ids = NULL,
                 keep_inputs = TRUE, cpca_blocks = NULL, cpca_T = NULL,
                 cpca_part = c("none", "design", "resid", "both"),
                 cpca_ridge = 0, ...) {
  cpca_part <- match.arg(cpca_part)
  omega_override <- NULL
  if (inherits(betas, "dkge_data")) {
    data <- betas
    omega_override <- omega
  } else {
    data <- dkge_data(betas, designs = designs, omega = omega, subject_ids = subject_ids)
  }

  K <- if (is.list(kernel) && !is.null(kernel$K)) kernel$K else kernel
  stopifnot(is.matrix(K), nrow(K) == data$q, ncol(K) == data$q)
  kernel_info <- if (is.list(kernel) && !is.null(kernel$info)) kernel$info else NULL

  fit <- dkge_fit(data, K = K, Omega_list = omega_override,
                     cpca_blocks = cpca_blocks, cpca_T = cpca_T,
                     cpca_part = cpca_part, cpca_ridge = cpca_ridge, ...)
  fit$kernel_info <- kernel_info
  if (keep_inputs) fit$input <- data
  fit
}
