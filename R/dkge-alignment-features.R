# dkge-alignment-features.R
# Over-ranked functional features and covariance-conditional residualization.

.dkge_alignment_feature_control_defaults <- function() {
  list(
    kernel_tolerance = 1e-10,
    rank_tolerance = 1e-9,
    max_condition = 1e8,
    min_feature_rank = 2L,
    min_residual_rank = 2L,
    min_effective_rank = 1.25,
    min_retained_energy = 0.05,
    max_row_collapse_fraction = 0.25,
    row_norm_tolerance = 1e-7,
    reconstruction_tolerance = 1e-8,
    covariance_orthogonality_tolerance = 1e-8
  )
}

#' Numerical gates for functional-alignment features
#'
#' These controls deliberately fail closed when a low-dimensional contrast
#' family consumes the useful functional feature space or when the conditional
#' residualization problem is not identified. No generalized inverse is used
#' to conceal a singular contrast covariance.
#'
#' @param kernel_tolerance Relative eigentolerance defining `image(K)`.
#' @param rank_tolerance Relative singular/eigenvalue tolerance for numerical
#'   ranks.
#' @param max_condition Largest permitted condition number for the contrast
#'   covariance.
#' @param min_feature_rank Minimum observed rank for independent features.
#' @param min_residual_rank Minimum numerical rank after residualization.
#' @param min_effective_rank Minimum participation-ratio rank after
#'   residualization.
#' @param min_retained_energy Minimum observed residual-to-original squared
#'   Frobenius energy ratio.
#' @param max_row_collapse_fraction Largest permitted fraction of initially
#'   nonzero rows that collapse numerically after residualization.
#' @param row_norm_tolerance Relative row-norm threshold used to diagnose
#'   collapse.
#' @param reconstruction_tolerance Relative tolerance for reconstructing each
#'   LOSO contrast as `G_s gamma_s`.
#' @param covariance_orthogonality_tolerance Relative tolerance for the
#'   analytic conditional-covariance oracle.
#' @return A typed alignment-feature control object.
#' @export
dkge_alignment_feature_control <- function(
    kernel_tolerance = 1e-10,
    rank_tolerance = 1e-9,
    max_condition = 1e8,
    min_feature_rank = 2L,
    min_residual_rank = 2L,
    min_effective_rank = 1.25,
    min_retained_energy = 0.05,
    max_row_collapse_fraction = 0.25,
    row_norm_tolerance = 1e-7,
    reconstruction_tolerance = 1e-8,
    covariance_orthogonality_tolerance = 1e-8) {
  out <- list(
    kernel_tolerance = kernel_tolerance,
    rank_tolerance = rank_tolerance,
    max_condition = max_condition,
    min_feature_rank = min_feature_rank,
    min_residual_rank = min_residual_rank,
    min_effective_rank = min_effective_rank,
    min_retained_energy = min_retained_energy,
    max_row_collapse_fraction = max_row_collapse_fraction,
    row_norm_tolerance = row_norm_tolerance,
    reconstruction_tolerance = reconstruction_tolerance,
    covariance_orthogonality_tolerance = covariance_orthogonality_tolerance
  )
  scalar_nonnegative <- c(
    "kernel_tolerance", "rank_tolerance", "min_effective_rank",
    "min_retained_energy", "max_row_collapse_fraction",
    "row_norm_tolerance", "reconstruction_tolerance",
    "covariance_orthogonality_tolerance"
  )
  for (nm in scalar_nonnegative) {
    value <- out[[nm]]
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value < 0) {
      .dkge_abort(sprintf("`%s` must be one finite non-negative number.", nm),
                  "dkge_alignment_feature_control_error")
    }
  }
  if (out$kernel_tolerance >= 1 || out$rank_tolerance >= 1 ||
      out$row_norm_tolerance >= 1 || out$reconstruction_tolerance >= 1 ||
      out$covariance_orthogonality_tolerance >= 1) {
    .dkge_abort("Relative alignment-feature tolerances must be smaller than one.",
                "dkge_alignment_feature_control_error")
  }
  if (out$min_retained_energy > 1 ||
      out$max_row_collapse_fraction > 1) {
    .dkge_abort("Energy and row-collapse fractions must lie in [0, 1].",
                "dkge_alignment_feature_control_error")
  }
  if (!is.numeric(out$max_condition) || length(out$max_condition) != 1L ||
      !is.finite(out$max_condition) || out$max_condition < 1) {
    .dkge_abort("`max_condition` must be one finite number at least one.",
                "dkge_alignment_feature_control_error")
  }
  for (nm in c("min_feature_rank", "min_residual_rank")) {
    value <- out[[nm]]
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value < 1 || value != as.integer(value)) {
      .dkge_abort(sprintf("`%s` must be one positive integer.", nm),
                  "dkge_alignment_feature_control_error")
    }
    out[[nm]] <- as.integer(value)
  }
  structure(out, class = c("dkge_alignment_feature_control", "list"))
}

.dkge_resolve_alignment_feature_control <- function(control = NULL) {
  if (is.null(control)) return(dkge_alignment_feature_control())
  if (inherits(control, "dkge_alignment_feature_control")) return(control)
  if (!is.list(control)) {
    .dkge_abort("`control` must be created by `dkge_alignment_feature_control()`.",
                "dkge_alignment_feature_control_error")
  }
  defaults <- .dkge_alignment_feature_control_defaults()
  unknown <- setdiff(names(control), names(defaults))
  if (length(unknown)) {
    .dkge_abort(
      sprintf("Unknown alignment-feature control field(s): %s.",
              paste(unknown, collapse = ", ")),
      "dkge_alignment_feature_control_error"
    )
  }
  do.call(dkge_alignment_feature_control,
          utils::modifyList(defaults, control, keep.null = TRUE))
}

.dkge_alignment_spectrum_diagnostics <- function(values, tolerance) {
  values <- as.numeric(values)
  positive_scale <- if (length(values)) max(values, 0) else 0
  threshold <- tolerance * positive_scale
  keep <- if (positive_scale > 0) values > threshold else rep(FALSE, length(values))
  kept <- pmax(values[keep], 0)
  rank <- length(kept)
  effective_rank <- if (rank && sum(kept^2) > 0) {
    sum(kept)^2 / sum(kept^2)
  } else {
    0
  }
  list(
    rank = as.integer(rank),
    effective_rank = as.numeric(min(rank, max(0, effective_rank))),
    condition = if (rank) max(kept) / min(kept) else Inf,
    threshold = threshold,
    keep = keep,
    values = values
  )
}

.dkge_validate_psd_matrix <- function(x, label, tolerance) {
  x <- as.matrix(x)
  if (!is.numeric(x) || nrow(x) < 1L || nrow(x) != ncol(x) ||
      any(!is.finite(x))) {
    .dkge_abort(sprintf("`%s` must be a finite non-empty square matrix.", label),
                "dkge_alignment_covariance_error")
  }
  scale <- max(abs(x))
  asymmetry <- max(abs(x - t(x)))
  if (asymmetry > tolerance * max(scale, .Machine$double.xmin)) {
    .dkge_abort(
      sprintf("`%s` is not symmetric within tolerance (error %.3e).",
              label, asymmetry),
      "dkge_alignment_covariance_error"
    )
  }
  x <- (x + t(x)) / 2
  eigenvalues <- eigen(x, symmetric = TRUE, only.values = TRUE)$values
  spectral_scale <- max(abs(eigenvalues), 0)
  negative_tolerance <- tolerance * spectral_scale
  if (length(eigenvalues) && min(eigenvalues) < -negative_tolerance) {
    .dkge_abort(
      sprintf("`%s` is not positive semidefinite (minimum eigenvalue %.3e).",
              label, min(eigenvalues)),
      "dkge_alignment_covariance_error"
    )
  }
  eigenvalues[eigenvalues < 0] <- 0
  list(
    matrix = x,
    spectrum = .dkge_alignment_spectrum_diagnostics(eigenvalues, tolerance)
  )
}

.dkge_matrix_singular_diagnostics <- function(x, tolerance) {
  singular_values <- svd(as.matrix(x), nu = 0L, nv = 0L)$d
  squared <- singular_values^2
  .dkge_alignment_spectrum_diagnostics(squared, tolerance)
}

#' Remove a contrast family from functional features conditionally
#'
#' For feature row `g`, generating directions `Gamma`, and effect-coordinate
#' covariance `Sigma_G`, this computes
#' `f = g - Sigma_G Gamma (Gamma' Sigma_G Gamma)^-1 Gamma' g`.
#' With row matrices this is `F = G - (G Gamma) C^-1 Gamma' Sigma_G`.
#' Under a separable spatial/effect covariance, the returned analytic oracle
#' verifies `Cov(F_p, v_q) = 0` for every parcel pair, up to the arbitrary
#' spatial scalar `rho[p,q]`.
#'
#' @param G `P` by `k` kernel-image feature matrix.
#' @param gamma `k` by `m` matrix of contrast-generating directions.
#' @param Sigma_G `k` by `k` effect-coordinate covariance.
#' @param control Numerical gates from [dkge_alignment_feature_control()].
#' @return A typed result containing residual features, reconstructed contrast
#'   values, the checked conditional coefficient, covariance oracle, and
#'   diagnostics.
#' @export
dkge_residualize_alignment_features <- function(G, gamma, Sigma_G,
                                                control = NULL) {
  started <- proc.time()[["elapsed"]]
  control <- .dkge_resolve_alignment_feature_control(control)
  G <- as.matrix(G)
  gamma <- as.matrix(gamma)
  if (!is.numeric(G) || nrow(G) < 1L || ncol(G) < 1L ||
      any(!is.finite(G))) {
    .dkge_abort("`G` must be a finite non-empty numeric matrix.",
                "dkge_alignment_feature_error")
  }
  if (!is.numeric(gamma) || nrow(gamma) != ncol(G) || ncol(gamma) < 1L ||
      any(!is.finite(gamma))) {
    .dkge_abort("`gamma` must be a finite k-by-m matrix aligned to `G`.",
                "dkge_alignment_rank_error")
  }
  covariance <- .dkge_validate_psd_matrix(
    Sigma_G, "Sigma_G", control$rank_tolerance
  )
  Sigma_G <- covariance$matrix
  if (nrow(Sigma_G) != ncol(G)) {
    .dkge_abort("`Sigma_G` must have one row and column per feature coordinate.",
                "dkge_alignment_covariance_error")
  }

  gamma_rank <- .dkge_matrix_singular_diagnostics(
    gamma, control$rank_tolerance
  )
  if (gamma_rank$rank < ncol(gamma)) {
    .dkge_abort(
      sprintf(
        paste0("The contrast family has estimable span %d but contains %d ",
               "columns; remove redundant/non-identifiable contrasts."),
        gamma_rank$rank, ncol(gamma)
      ),
      "dkge_alignment_rank_error"
    )
  }
  available_dimension <- ncol(G) - gamma_rank$rank
  if (available_dimension < control$min_residual_rank) {
    .dkge_abort(
      sprintf(
        paste0("Residual feature dimension is %d = rank(K) %d minus ",
               "contrast span %d, below the required floor %d."),
        available_dimension, ncol(G), gamma_rank$rank,
        control$min_residual_rank
      ),
      "dkge_alignment_rank_error"
    )
  }

  Sigma_gamma <- Sigma_G %*% gamma
  contrast_covariance <- crossprod(gamma, Sigma_gamma)
  contrast_covariance <- (contrast_covariance + t(contrast_covariance)) / 2
  contrast_covariance_check <- .dkge_validate_psd_matrix(
    contrast_covariance, "contrast covariance gamma' Sigma_G gamma",
    control$rank_tolerance
  )
  c_spectrum <- contrast_covariance_check$spectrum
  if (c_spectrum$rank < ncol(gamma)) {
    .dkge_abort(
      sprintf(
        "Contrast covariance rank is %d for %d generating directions.",
        c_spectrum$rank, ncol(gamma)
      ),
      "dkge_alignment_covariance_error"
    )
  }
  if (!is.finite(c_spectrum$condition) ||
      c_spectrum$condition > control$max_condition) {
    .dkge_abort(
      sprintf(
        "Contrast covariance condition %.3e exceeds the gate %.3e.",
        c_spectrum$condition, control$max_condition
      ),
      "dkge_alignment_covariance_error"
    )
  }

  chol_covariance <- tryCatch(
    chol(contrast_covariance_check$matrix),
    error = function(e) NULL
  )
  if (is.null(chol_covariance)) {
    .dkge_abort(
      "Contrast covariance failed the checked Cholesky solve.",
      "dkge_alignment_covariance_error"
    )
  }
  rhs <- crossprod(gamma, Sigma_G)
  coefficient <- backsolve(
    chol_covariance,
    forwardsolve(t(chol_covariance), rhs)
  )
  values <- G %*% gamma
  features <- G - values %*% coefficient
  residual_covariance <- Sigma_G - Sigma_gamma %*% coefficient
  residual_covariance <- (residual_covariance + t(residual_covariance)) / 2
  residual_covariance_check <- .dkge_validate_psd_matrix(
    residual_covariance, "residual feature covariance",
    control$rank_tolerance * 10
  )
  residual_spectrum <- residual_covariance_check$spectrum
  observed_spectrum <- .dkge_matrix_singular_diagnostics(
    features, control$rank_tolerance
  )
  if (residual_spectrum$rank < control$min_residual_rank ||
      observed_spectrum$rank < control$min_residual_rank) {
    .dkge_abort(
      sprintf(
        paste0("Residual numerical rank failed: covariance rank %d, observed ",
               "rank %d, required %d."),
        residual_spectrum$rank, observed_spectrum$rank,
        control$min_residual_rank
      ),
      "dkge_alignment_rank_error"
    )
  }
  effective_rank <- min(residual_spectrum$effective_rank,
                        observed_spectrum$effective_rank)
  if (effective_rank < control$min_effective_rank) {
    .dkge_abort(
      sprintf("Residual effective rank %.3f is below the gate %.3f.",
              effective_rank, control$min_effective_rank),
      "dkge_alignment_rank_error"
    )
  }

  original_energy <- sum(G^2)
  residual_energy <- sum(features^2)
  retained_energy <- if (original_energy > 0) {
    residual_energy / original_energy
  } else {
    0
  }
  if (!is.finite(retained_energy) ||
      retained_energy < control$min_retained_energy) {
    .dkge_abort(
      sprintf("Residual feature energy %.3f is below the gate %.3f.",
              retained_energy, control$min_retained_energy),
      "dkge_alignment_feature_collapse_error"
    )
  }
  original_row_norm <- sqrt(rowSums(G^2))
  residual_row_norm <- sqrt(rowSums(features^2))
  row_scale <- max(original_row_norm, 0)
  active_rows <- original_row_norm > control$row_norm_tolerance * row_scale
  if (!any(active_rows)) {
    .dkge_abort("Every functional-feature row is numerically zero.",
                "dkge_alignment_feature_collapse_error")
  }
  row_ratio <- residual_row_norm[active_rows] /
    pmax(original_row_norm[active_rows], .Machine$double.xmin)
  collapsed_fraction <- mean(row_ratio <= control$row_norm_tolerance)
  if (collapsed_fraction > control$max_row_collapse_fraction) {
    .dkge_abort(
      sprintf("Residualization collapsed %.1f%% of active rows; gate is %.1f%%.",
              100 * collapsed_fraction,
              100 * control$max_row_collapse_fraction),
      "dkge_alignment_feature_collapse_error"
    )
  }

  covariance_cross <- residual_covariance %*% gamma
  covariance_scale <- max(norm(Sigma_gamma, "F"), .Machine$double.xmin)
  covariance_error <- norm(covariance_cross, "F") / covariance_scale
  if (!is.finite(covariance_error) ||
      covariance_error > control$covariance_orthogonality_tolerance) {
    .dkge_abort(
      sprintf(
        "Conditional covariance oracle error %.3e exceeds tolerance %.3e.",
        covariance_error, control$covariance_orthogonality_tolerance
      ),
      "dkge_alignment_covariance_error"
    )
  }

  theoretical_energy <- sum(pmax(residual_spectrum$values, 0)) /
    max(sum(pmax(covariance$spectrum$values, 0)), .Machine$double.xmin)
  diagnostics <- list(
    kernel_support_dimension = ncol(G),
    contrast_count = ncol(gamma),
    contrast_span_rank = gamma_rank$rank,
    contrast_covariance_rank = c_spectrum$rank,
    contrast_covariance_condition = c_spectrum$condition,
    available_dimension = available_dimension,
    residual_covariance_rank = residual_spectrum$rank,
    residual_observed_rank = observed_spectrum$rank,
    residual_effective_rank = effective_rank,
    retained_energy = retained_energy,
    theoretical_retained_covariance = theoretical_energy,
    collapsed_row_fraction = collapsed_fraction,
    covariance_orthogonality_error = covariance_error,
    runtime_seconds = proc.time()[["elapsed"]] - started,
    gates_passed = TRUE,
    conditional_independence = paste0(
      "Cov(F[p,], v[q,]) = 0 for every parcel pair under separable ",
      "spatial/effect covariance. Gaussian independence is claimed only ",
      "conditional on a fixed generating direction gamma."
    )
  )
  structure(
    list(
      features = features,
      values = values,
      coefficient = coefficient,
      contrast_covariance = contrast_covariance,
      residual_covariance = residual_covariance,
      covariance_cross_oracle = covariance_cross,
      diagnostics = diagnostics,
      control = control
    ),
    class = c("dkge_residualized_alignment_features", "list")
  )
}

.dkge_compact_kernel_factor <- function(K, tolerance) {
  geometry <- .dkge_kernel_geometry(K, tol = tolerance)
  if (geometry$rank < 1L) {
    .dkge_abort("The design kernel has empty numerical support.",
                "dkge_alignment_rank_error")
  }
  support <- which(geometry$support)
  L <- sweep(
    geometry$evecs[, support, drop = FALSE], 2L,
    sqrt(geometry$evals[support]), "*"
  )
  rownames(L) <- rownames(geometry$K)
  colnames(L) <- paste0("kernel_image_", seq_len(ncol(L)))
  reconstruction_error <- norm(L %*% t(L) - geometry$K, "F") /
    max(norm(geometry$K, "F"), .Machine$double.xmin)
  list(
    factor = L,
    support_vectors = geometry$evecs[, support, drop = FALSE],
    geometry = geometry,
    reconstruction_error = reconstruction_error
  )
}

.dkge_alignment_fit_binding <- function(fit) {
  .dkge_object_hash(list(
    fit_class = class(fit),
    subject_ids = as.character(fit$subject_ids %||% seq_along(fit$Btil)),
    effects = fit$effects,
    kernel = list(
      K = fit$K,
      Khalf = fit$Khalf,
      Kihalf = fit$Kihalf,
      support_projector = fit$kernel_support_projector,
      info = fit$kernel_info,
      diagnostics = fit$kernel_diagnostics,
      rank = fit$kernel_rank,
      nullity = fit$kernel_nullity,
      condition = fit$kernel_condition
    ),
    data = list(
      R = fit$R,
      Omega = fit$Omega,
      effect_scaling = fit$effect_scaling,
      Braw = fit$Braw,
      Btil = fit$Btil,
      subjects = fit$subjects,
      provenance = fit$provenance
    ),
    estimator = list(
      U = fit$U,
      KU = fit$KU,
      Chat = fit$Chat,
      Chat_sym = fit$Chat_sym,
      evals = fit$evals,
      sdev = fit$sdev,
      eig_vectors_full = fit$eig_vectors_full,
      eig_values_full = fit$eig_values_full,
      rank = fit$rank %||% ncol(fit$U),
      rank_requested = fit$rank_requested,
      effective_rank = fit$effective_rank,
      rank_reduced = fit$rank_reduced,
      solver = fit$solver %||% "pooled",
      cpca = fit$cpca,
      jd = fit$jd,
      representation = fit$representation,
      representation_reasons = fit$representation_reasons,
      ridge_input = fit$ridge_input %||% 0
    ),
    pooling = list(
      subject_weights = fit$weights,
      subject_weight_scores_raw = fit$subject_weight_scores_raw,
      subject_weight_usable = fit$subject_weight_usable,
      w_method = fit$w_method,
      w_tau = fit$w_tau,
      weight_spec = fit$weight_spec,
      effect_weight_spec = fit$effect_weight_spec,
      effect_precision = fit$effect_precision,
      effect_precision_diagnostics = fit$effect_precision_diagnostics,
      voxel_weights = fit$voxel_weights,
      voxel_weights_subject = fit$voxel_weights_subject,
      voxel_weights_prior = fit$voxel_weights_prior,
      voxel_weights_adapt = fit$voxel_weights_adapt,
      missingness = fit$missingness,
      miss_args = fit$miss_args,
      debias = fit$debias,
      pool_cache = fit$pool_cache
    ),
    moments = list(
      contribs = fit$contribs,
      effect_moment = fit$effect_moment,
      effect_moments = fit$effect_moments,
      effect_moments_raw = fit$effect_moments_raw,
      noise_moments = fit$noise_moments,
      pair_counts = fit$pair_counts,
      pair_weight = fit$pair_weight,
      pair_ess = fit$pair_ess,
      diagnostics = fit$moment_diagnostics
    ),
    spatial = .dkge_spatial_fit_payload(fit$spatial, validate = TRUE)
  ))
}

.dkge_alignment_contrast_binding <- function(contrast_obj) {
  .dkge_object_hash(list(
    method = contrast_obj$method,
    contrasts = contrast_obj$contrasts,
    values = contrast_obj$values,
    receipts = contrast_obj$metadata$alignment_receipts %||% NULL
  ))
}

.dkge_alignment_contrast_family_binding <- function(contrasts) {
  if (inherits(contrasts, "dkge_contrasts")) {
    contrasts <- contrasts$contrasts
  }
  if (!is.list(contrasts) || !length(contrasts)) {
    .dkge_abort("A non-empty normalized contrast family is required.",
                "dkge_alignment_feature_error")
  }
  .dkge_object_hash(list(
    contrast_ids = names(contrasts),
    contrasts = lapply(contrasts, function(x) {
      list(
        values = as.numeric(x),
        effect_ids = names(x),
        scope = attr(x, "dkge_scope", exact = TRUE),
        term = attr(x, "dkge_term", exact = TRUE)
      )
    })
  ))
}

.dkge_alignment_features_payload <- function(x) {
  list(
    schema_version = x$schema_version,
    feature_source = x$feature_source,
    source_verified = x$source_verified,
    estimator_source = x$estimator_source,
    subject_ids = x$subject_ids,
    contrast_ids = x$contrast_ids,
    cluster_ids = x$cluster_ids,
    features = x$features,
    kernel_factor = x$kernel_factor,
    kernel_rank = x$kernel_rank,
    kernel_hash = x$kernel_hash,
    generating_directions = x$generating_directions,
    feature_covariances = x$feature_covariances,
    diagnostics = x$diagnostics,
    control = x$control,
    fit_binding = x$fit_binding,
    contrast_binding = x$contrast_binding,
    contrast_family_binding = x$contrast_family_binding,
    provenance = x$provenance,
    eligibility = x$eligibility
  )
}

.dkge_validate_alignment_features <- function(x, fit = NULL,
                                              contrast_obj = NULL) {
  if (!inherits(x, "dkge_alignment_features")) {
    .dkge_abort("Expected a typed `dkge_alignment_features` object.",
                "dkge_alignment_feature_error")
  }
  if (!is.null(fit)) {
    .dkge_assert_crossfit_estimator_supported(
      fit, "Alignment-feature validation"
    )
  }
  expected <- .dkge_object_hash(.dkge_alignment_features_payload(x))
  if (!identical(expected, x$structural_hash)) {
    .dkge_abort("Alignment features were mutated after construction.",
                "dkge_alignment_feature_error")
  }
  if (!is.null(fit) &&
      !identical(x$fit_binding, .dkge_alignment_fit_binding(fit))) {
    .dkge_abort("Alignment features and the live DKGE fit have mixed provenance.",
                "dkge_alignment_feature_error")
  }
  if (!is.null(contrast_obj) &&
      !identical(x$contrast_binding,
                 .dkge_alignment_contrast_binding(contrast_obj))) {
    .dkge_abort("Alignment features were built for a different contrast result.",
                "dkge_alignment_feature_error")
  }
  if (!is.null(contrast_obj) &&
      !identical(x$contrast_family_binding,
                 .dkge_alignment_contrast_family_binding(contrast_obj))) {
    .dkge_abort("Alignment features were built for a different contrast family.",
                "dkge_alignment_feature_error")
  }
  invisible(x)
}

.dkge_order_subject_list <- function(x, subject_ids, label) {
  if (!is.list(x) || length(x) != length(subject_ids)) {
    .dkge_abort(
      sprintf("`%s` must contain one entry per fitted subject.", label),
      "dkge_alignment_feature_error"
    )
  }
  nms <- names(x)
  if (!is.null(nms)) {
    if (length(nms) != length(subject_ids) || anyNA(nms) ||
        any(!nzchar(nms)) || anyDuplicated(nms) ||
        !setequal(nms, subject_ids)) {
      .dkge_abort(
        sprintf("Named `%s` entries must match fitted subject IDs exactly.", label),
        "dkge_alignment_feature_error"
      )
    }
    idx <- match(subject_ids, nms)
    x <- x[idx]
  }
  names(x) <- subject_ids
  x
}

.dkge_validate_independent_beta <- function(B, fit, s) {
  B <- as.matrix(B)
  q <- nrow(fit$K)
  if (!is.numeric(B) || nrow(B) != q || any(!is.finite(B))) {
    .dkge_abort(
      sprintf("Independent beta block for subject '%s' must be finite and q-by-P.",
              fit$subject_ids[[s]]),
      "dkge_alignment_feature_error"
    )
  }
  effects <- fit$effects %||% rownames(fit$K)
  if (!is.null(effects) && !is.null(rownames(B))) {
    idx <- match(effects, rownames(B))
    if (anyNA(idx) || !setequal(effects, rownames(B))) {
      .dkge_abort("Independent beta effect labels do not match the fitted design.",
                  "dkge_alignment_feature_error")
    }
    B <- B[idx, , drop = FALSE]
  }
  expected_clusters <- colnames(fit$Btil[[s]]) %||% seq_len(ncol(fit$Btil[[s]]))
  if (ncol(B) != length(expected_clusters)) {
    .dkge_abort(
      sprintf("Independent beta block for subject '%s' has the wrong parcel count.",
              fit$subject_ids[[s]]),
      "dkge_alignment_feature_error"
    )
  }
  if (!is.null(colnames(B)) && !is.null(colnames(fit$Btil[[s]]))) {
    idx <- match(expected_clusters, colnames(B))
    if (anyNA(idx) || !setequal(expected_clusters, colnames(B))) {
      .dkge_abort("Independent beta parcel labels do not match the fitted support.",
                  "dkge_alignment_feature_error")
    }
    B <- B[, idx, drop = FALSE]
  }
  rownames(B) <- effects %||% rownames(B)
  colnames(B) <- as.character(expected_clusters)
  B
}

.dkge_independent_feature_diagnostics <- function(G, control, kernel_rank) {
  spectrum <- .dkge_matrix_singular_diagnostics(G, control$rank_tolerance)
  if (kernel_rank < control$min_feature_rank ||
      spectrum$rank < control$min_feature_rank) {
    .dkge_abort(
      sprintf("Independent functional features have rank %d; required floor is %d.",
              spectrum$rank, control$min_feature_rank),
      "dkge_alignment_rank_error"
    )
  }
  norms <- sqrt(rowSums(G^2))
  scale <- max(norms, 0)
  zero_fraction <- if (scale > 0) {
    mean(norms <= control$row_norm_tolerance * scale)
  } else {
    1
  }
  if (zero_fraction > control$max_row_collapse_fraction) {
    .dkge_abort(
      sprintf("Independent features contain %.1f%% collapsed rows; gate is %.1f%%.",
              100 * zero_fraction,
              100 * control$max_row_collapse_fraction),
      "dkge_alignment_feature_collapse_error"
    )
  }
  list(
    kernel_support_dimension = kernel_rank,
    available_dimension = kernel_rank,
    residual_observed_rank = spectrum$rank,
    residual_effective_rank = spectrum$effective_rank,
    retained_energy = 1,
    collapsed_row_fraction = zero_fraction,
    gates_passed = TRUE
  )
}

.dkge_validate_receipt_preprocessing <- function(fit, receipt, s) {
  p <- receipt$preprocessing
  Btil <- fit$Btil[[s]]
  if (!is.list(p) || !identical(p$schema_version, "2.1.0") ||
      !identical(p$standardized_beta_source, "fit$Btil") ||
      !identical(p$beta_hash, .dkge_object_hash(Btil)) ||
      !identical(p$ruler_hash,
                 .dkge_object_hash(fit$R %||% diag(nrow(Btil)))) ||
      !identical(p$voxel_weights_hash,
                 .dkge_object_hash(p$voxel_weights %||% NULL)) ||
      !isTRUE(p$effect_separability$transform_verified)) {
    .dkge_abort(
      sprintf("Unverified alignment preprocessing for subject '%s'.",
              receipt$subject_id),
      "dkge_alignment_preprocessing_error"
    )
  }
  if (!identical(p$effect_scaling, fit$effect_scaling) ||
      !p$effect_scaling %in% c("pooled_design", "none")) {
    .dkge_abort("The recorded effect-space transformation is unverified.",
                "dkge_alignment_preprocessing_error")
  }
  expected_cluster_order <- colnames(Btil) %||% seq_len(ncol(Btil))
  if (!identical(as.character(p$cluster_order),
                 as.character(expected_cluster_order)) ||
      !identical(p$cluster_order_hash,
                 .dkge_object_hash(p$cluster_order))) {
    .dkge_abort("Alignment receipt parcel order is invalid.",
                "dkge_alignment_preprocessing_error")
  }
  op <- .dkge_fit_spatial_operator(fit, subject = s, n_cols = ncol(Btil))
  active <- !is.null(op) && isTRUE(op$lambda > 0)
  if (!identical(isTRUE(p$spatial$active), active) ||
      (active && (!identical(p$spatial$fingerprint, op$fingerprint) ||
                  !identical(p$spatial$operator_binding,
                             .dkge_spatial_operator_binding(op)) ||
                  !isTRUE(all.equal(p$spatial$lambda, as.numeric(op$lambda),
                                    tolerance = 0))))) {
    .dkge_abort("The recorded spatial preprocessing has mixed provenance.",
                "dkge_alignment_preprocessing_error")
  }
  invisible(TRUE)
}

.dkge_alignment_receipts_for_features <- function(fit, contrast_obj) {
  receipts <- contrast_obj$metadata$alignment_receipts %||% NULL
  if (!inherits(receipts, "dkge_alignment_receipts") || !length(receipts)) {
    .dkge_abort(
      paste0("Same-data residualized features require typed LOSO/K-fold ",
             "alignment receipts."),
      "dkge_alignment_receipt_error"
    )
  }
  subject_ids <- as.character(fit$subject_ids %||% seq_along(fit$Btil))
  if (length(receipts) != length(subject_ids) ||
      !identical(unname(vapply(receipts, `[[`, character(1), "subject_id")),
                 subject_ids)) {
    .dkge_abort("Alignment receipts do not match the fitted subject cohort.",
                "dkge_alignment_receipt_error")
  }
  for (s in seq_along(receipts)) {
    receipt <- receipts[[s]]
    if (!inherits(receipt, "dkge_alignment_receipt") ||
        !identical(receipt$basis_hash, .dkge_object_hash(receipt$basis)) ||
        !identical(receipt$loadings_hash, .dkge_object_hash(receipt$loadings)) ||
        !isTRUE(receipt$inference$eligible)) {
      .dkge_abort(
        sprintf("Alignment receipt for subject '%s' is invalid or ineligible.",
                subject_ids[[s]]),
        "dkge_alignment_receipt_error"
      )
    }
    .dkge_validate_receipt_preprocessing(fit, receipt, s)
  }
  receipts
}

.dkge_resolve_effect_noise_covariances <- function(fit, effect_noise_cov) {
  subject_ids <- as.character(fit$subject_ids %||% seq_along(fit$Btil))
  if (is.null(effect_noise_cov)) {
    effect_noise_cov <- lapply(fit$subjects %||% vector("list", length(subject_ids)),
                               `[[`, "effect_noise_cov")
    source <- "fit_subject_effect_noise_cov"
  } else {
    source <- "explicit_effect_noise_cov"
  }
  effect_noise_cov <- .dkge_order_subject_list(
    effect_noise_cov, subject_ids, "effect_noise_cov"
  )
  missing <- vapply(effect_noise_cov, is.null, logical(1))
  if (any(missing)) {
    .dkge_abort(
      sprintf("Missing effect covariance for subject(s): %s.",
              paste(subject_ids[missing], collapse = ", ")),
      "dkge_alignment_covariance_error"
    )
  }
  q <- nrow(fit$K)
  effects <- fit$effects %||% rownames(fit$K)
  out <- lapply(seq_along(effect_noise_cov), function(s) {
    Lambda <- as.matrix(effect_noise_cov[[s]])
    if (!is.numeric(Lambda) || !identical(dim(Lambda), c(q, q)) ||
        any(!is.finite(Lambda))) {
      .dkge_abort(
        sprintf("Effect covariance for subject '%s' must be finite q-by-q.",
                subject_ids[[s]]),
        "dkge_alignment_covariance_error"
      )
    }
    if (!is.null(effects) && !is.null(rownames(Lambda))) {
      idx <- match(effects, rownames(Lambda))
      if (anyNA(idx) || is.null(colnames(Lambda)) ||
          !setequal(effects, colnames(Lambda))) {
        .dkge_abort("Effect covariance labels do not match fitted effects.",
                    "dkge_alignment_covariance_error")
      }
      Lambda <- Lambda[idx, match(effects, colnames(Lambda)), drop = FALSE]
    }
    .dkge_validate_psd_matrix(
      Lambda, sprintf("effect_noise_cov[[%d]]", s), 1e-9
    )$matrix
  })
  names(out) <- subject_ids
  list(covariances = out, source = source)
}

.dkge_alignment_feature_preprocessing <- function(x) {
  feature_provenance <- x$provenance
  feature_provenance$feature_object_hash <- x$structural_hash
  list(
    source = paste0("kernel_image_", x$feature_source),
    feature_source = x$feature_source,
    estimator_source = x$estimator_source,
    fit_binding = x$fit_binding,
    contrast_binding = x$contrast_binding,
    contrast_family_binding = x$contrast_family_binding,
    recompute_under_null = FALSE,
    feature_provenance = feature_provenance,
    inferentially_eligible = isTRUE(x$eligibility$eligible),
    eligibility_status = x$eligibility$status,
    contract = x$provenance$inferential_contract
  )
}

#' Construct over-ranked functional-alignment features
#'
#' Functional correspondence is represented in compact coordinates for
#' `image(K)`: if `K = L_K L_K'`, subject features are
#' `G_s = Bmodel_s' L_K`. This retains `rank(K)` coordinates independently of
#' the low fitted estimation rank. With no parcel preprocessing,
#' `Btil_s' Khalf = G_s V_K'`, so all Euclidean feature distances are exactly
#' preserved by the compact representation.
#'
#' `feature_source = "independent"` is the default and requires separate beta
#' maps plus an explicit independent-data identifier. The independent maps are
#' transformed by the design-only fitted ruler but never by beta-adaptive voxel
#' weights or fitted spatial smoothing.
#'
#' `feature_source = "same_data_residualized"` is an opt-in approximate mode.
#' It reconstructs each LOSO contrast as `G_s gamma_s`, then conditions the
#' whole feature matrix on the joint contrast family using
#' `Sigma_G,s = L_K' R' Lambda_s R L_K`. Its covariance oracle is exact under
#' the stated separable model conditional on fixed `gamma_s`; because
#' `gamma_s` depends on the estimated rank-truncated held-out basis, frozen-plan
#' sign-flip inference remains approximate. Full re-estimation under every
#' valid sign action is the exact same-data reference. The preprocessing receipt
#' verifies the implemented effect/parcel transformations; it records
#' separability as a model assumption and does not empirically establish it.
#' Both feature modes require contrasts from a pooled, non-CPCA fit. CPCA and
#' JD estimators fail closed until their exact training-fold estimator can be
#' replayed.
#'
#' @param fit A fitted `dkge` object.
#' @param contrast_obj A `dkge_contrasts` result to which the features are
#'   immutably bound.
#' @param feature_source Either independent features (recommended/default) or
#'   covariance-residualized same-data features.
#' @param independent_betas For independent mode, one raw q-by-P beta matrix per
#'   fitted subject, on the same effect and parcel ordering.
#' @param independent_data_hash Stable non-empty identifier for the independent
#'   acquisition/training data. Supplying the primary fitted betas is rejected.
#' @param effect_noise_cov Optional list of q-by-q subject effect covariances for
#'   same-data residualization. By default these are read from `fit$subjects`.
#' @param control Numerical gates from [dkge_alignment_feature_control()].
#' @return An immutable `dkge_alignment_features` object. Its `$features` list
#'   can be supplied directly to the transport fitter through the typed object.
#' @export
dkge_alignment_features <- function(
    fit,
    contrast_obj,
    feature_source = c("independent", "same_data_residualized"),
    independent_betas = NULL,
    independent_data_hash = NULL,
    effect_noise_cov = NULL,
    control = NULL) {
  if (!inherits(fit, "dkge") || !inherits(contrast_obj, "dkge_contrasts")) {
    .dkge_abort("`fit` and `contrast_obj` must be DKGE fitted/contrast objects.",
                "dkge_alignment_feature_error")
  }
  .dkge_assert_crossfit_estimator_supported(
    fit, "`dkge_alignment_features()`"
  )
  feature_source <- match.arg(feature_source)
  control <- .dkge_resolve_alignment_feature_control(control)
  if (is.null(fit$K) || is.null(fit$R) || !length(fit$Btil)) {
    .dkge_abort("The DKGE fit lacks K, R, or standardized beta blocks.",
                "dkge_alignment_feature_error")
  }
  subject_ids <- as.character(fit$subject_ids %||% seq_along(fit$Btil))
  fit$subject_ids <- subject_ids
  contrast_ids <- names(contrast_obj$values) %||%
    names(contrast_obj$contrasts) %||%
    paste0("contrast", seq_along(contrast_obj$values))
  if (!length(contrast_ids) || any(!nzchar(contrast_ids)) ||
      anyDuplicated(contrast_ids)) {
    .dkge_abort("Contrast identifiers must be non-empty and unique.",
                "dkge_alignment_feature_error")
  }
  kernel <- .dkge_compact_kernel_factor(
    fit$K, tolerance = control$kernel_tolerance
  )
  L <- kernel$factor
  if (ncol(L) < control$min_feature_rank) {
    .dkge_abort(
      sprintf("Kernel support dimension %d is below the feature floor %d.",
              ncol(L), control$min_feature_rank),
      "dkge_alignment_rank_error"
    )
  }
  features <- vector("list", length(subject_ids))
  directions <- vector("list", length(subject_ids))
  covariances <- vector("list", length(subject_ids))
  diagnostics <- vector("list", length(subject_ids))
  cluster_ids <- lapply(fit$Btil, function(B) {
    as.character(colnames(B) %||% seq_len(ncol(B)))
  })
  names(features) <- names(directions) <- names(covariances) <-
    names(diagnostics) <- names(cluster_ids) <- subject_ids

  if (identical(feature_source, "independent")) {
    if (!is.character(independent_data_hash) ||
        length(independent_data_hash) != 1L || is.na(independent_data_hash) ||
        !nzchar(independent_data_hash)) {
      .dkge_abort(
        "Independent alignment requires one non-empty `independent_data_hash`.",
        "dkge_alignment_feature_error"
      )
    }
    independent_betas <- .dkge_order_subject_list(
      independent_betas, subject_ids, "independent_betas"
    )
    independent_betas <- lapply(seq_along(independent_betas), function(s) {
      .dkge_validate_independent_beta(independent_betas[[s]], fit, s)
    })
    names(independent_betas) <- subject_ids
    independent_content_hash <- .dkge_object_hash(
      unname(lapply(independent_betas, unname))
    )
    fitted_content_hash <- if (length(fit$Braw)) {
      .dkge_object_hash(unname(lapply(fit$Braw, unname)))
    } else {
      NA_character_
    }
    if (identical(independent_content_hash, fitted_content_hash)) {
      .dkge_abort(
        "The claimed independent beta maps are identical to the fitted data.",
        "dkge_alignment_feature_error"
      )
    }
    if (!fit$effect_scaling %in% c("pooled_design", "none")) {
      .dkge_abort("The fitted effect transformation is unverified.",
                  "dkge_alignment_preprocessing_error")
    }
    independent_model <- if (identical(fit$effect_scaling, "pooled_design")) {
      .dkge_row_standardize(independent_betas, fit$R)
    } else {
      independent_betas
    }
    for (s in seq_along(subject_ids)) {
      G <- t(independent_model[[s]]) %*% L
      colnames(G) <- colnames(L)
      rownames(G) <- cluster_ids[[s]]
      diagnostics[[s]] <- .dkge_independent_feature_diagnostics(
        G, control, ncol(L)
      )
      Z <- t(independent_model[[s]]) %*% fit$Khalf
      diagnostics[[s]]$full_square_root_equivalence_error <-
        norm(Z - G %*% t(kernel$support_vectors), "F") /
        max(norm(Z, "F"), .Machine$double.xmin)
      diagnostics[[s]]$kernel_reconstruction_error <-
        kernel$reconstruction_error
      features[[s]] <- G
      directions[s] <- list(NULL)
      covariances[s] <- list(NULL)
    }
    source_verified <- TRUE
    provenance <- list(
      independent_data_hash = independent_data_hash,
      independent_content_hash = independent_content_hash,
      primary_content_hash = fitted_content_hash,
      standardization = if (identical(fit$effect_scaling, "pooled_design")) {
        "design-only R_transpose"
      } else {
        "identity"
      },
      voxel_weighting = "none",
      spatial_preprocessing = "none",
      inferential_contract = paste0(
        "Independent functional correspondence. Group exactness still depends ",
        "on the contrast estimator; a same-data rank-truncated estimator is ",
        "classified as approximate by the frozen court."
      )
    )
  } else {
    if (!is.null(independent_betas) || !is.null(independent_data_hash)) {
      .dkge_abort(
        "Independent-data arguments are not used in same-data residualized mode.",
        "dkge_alignment_feature_error"
      )
    }
    receipts <- .dkge_alignment_receipts_for_features(fit, contrast_obj)
    covariance_info <- .dkge_resolve_effect_noise_covariances(
      fit, effect_noise_cov
    )
    R <- fit$R
    for (s in seq_along(subject_ids)) {
      receipt <- receipts[[s]]
      alpha_names <- names(receipt$alphas)
      if (is.null(alpha_names)) alpha_names <- contrast_ids
      idx <- match(contrast_ids, alpha_names)
      if (anyNA(idx)) {
        .dkge_abort(
          sprintf("Receipt for subject '%s' does not cover the contrast family.",
                  subject_ids[[s]]),
          "dkge_alignment_receipt_error"
        )
      }
      Bmodel <- .dkge_fit_model_btil(
        fit, s, voxel_weights = receipt$preprocessing$voxel_weights
      )
      G <- t(Bmodel) %*% L
      Gamma <- do.call(cbind, lapply(receipt$alphas[idx], function(alpha) {
        as.numeric(crossprod(L, receipt$basis %*% as.numeric(alpha)))
      }))
      colnames(Gamma) <- contrast_ids
      rownames(Gamma) <- colnames(L)
      Lambda <- covariance_info$covariances[[s]]
      Sigma_G <- crossprod(L, t(R) %*% Lambda %*% R %*% L)
      residual <- dkge_residualize_alignment_features(
        G, Gamma, Sigma_G, control = control
      )
      expected <- do.call(cbind, lapply(contrast_obj$values, function(values) {
        as.numeric(values[[s]])
      }))
      colnames(expected) <- contrast_ids
      reconstruction_error <- max(abs(residual$values - expected)) /
        max(1, max(abs(expected)))
      if (!is.finite(reconstruction_error) ||
          reconstruction_error > control$reconstruction_tolerance) {
        .dkge_abort(
          sprintf(
            "G_s gamma_s reconstruction failed for subject '%s' (error %.3e).",
            subject_ids[[s]], reconstruction_error
          ),
          "dkge_alignment_reconstruction_error"
        )
      }
      colnames(residual$features) <- colnames(L)
      rownames(residual$features) <- cluster_ids[[s]]
      residual$diagnostics$reconstruction_error <- reconstruction_error
      Z <- t(Bmodel) %*% fit$Khalf
      residual$diagnostics$full_square_root_equivalence_error <-
        norm(Z - G %*% t(kernel$support_vectors), "F") /
        max(norm(Z, "F"), .Machine$double.xmin)
      residual$diagnostics$kernel_reconstruction_error <-
        kernel$reconstruction_error
      features[[s]] <- residual$features
      directions[[s]] <- Gamma
      covariances[[s]] <- Sigma_G
      diagnostics[[s]] <- residual$diagnostics
    }
    source_verified <- TRUE
    provenance <- list(
      covariance_source = covariance_info$source,
      covariance_hash = .dkge_object_hash(covariance_info$covariances),
      preprocessing = paste(
        "receipt-validated transforms under a declared, not empirically",
        "verified, effect-separability assumption"
      ),
      residualization = "joint covariance-conditional contrast family",
      inferential_contract = paste0(
        "Same-data residual features are covariance-orthogonal to the entire ",
        "contrast vector under separability, conditional on fixed generating ",
        "directions. Those directions depend on estimated rank-truncated fold ",
        "bases, so frozen-plan inference is approximate; only full null-action ",
        "re-estimation is exact."
      )
    )
  }

  estimator_source <- "same_data_rank_truncated"
  eligibility <- dkge_alignment_eligibility(
    feature_source = feature_source,
    estimator_source = estimator_source,
    recompute_under_null = FALSE,
    source_verified = source_verified
  )
  out <- list(
    schema_version = "1.0.0",
    feature_source = feature_source,
    source_verified = source_verified,
    estimator_source = estimator_source,
    subject_ids = subject_ids,
    contrast_ids = contrast_ids,
    cluster_ids = cluster_ids,
    features = features,
    kernel_factor = L,
    kernel_rank = ncol(L),
    kernel_hash = .dkge_object_hash(fit$K),
    generating_directions = directions,
    feature_covariances = covariances,
    diagnostics = diagnostics,
    control = control,
    fit_binding = .dkge_alignment_fit_binding(fit),
    contrast_binding = .dkge_alignment_contrast_binding(contrast_obj),
    contrast_family_binding = .dkge_alignment_contrast_family_binding(
      contrast_obj
    ),
    provenance = provenance,
    eligibility = eligibility
  )
  out$structural_hash <- .dkge_object_hash(.dkge_alignment_features_payload(out))
  structure(out, class = c("dkge_alignment_features", "list"))
}

#' @export
print.dkge_alignment_features <- function(x, ...) {
  .dkge_validate_alignment_features(x)
  available <- vapply(x$diagnostics, `[[`, numeric(1), "available_dimension")
  cat("<dkge_alignment_features>\n")
  cat("  source          :", x$feature_source, "\n")
  cat("  subjects        :", length(x$subject_ids), "\n")
  cat("  kernel dimension:", x$kernel_rank, "\n")
  cat("  available range :", paste(range(available), collapse = "-"), "\n")
  cat("  eligibility     :", x$eligibility$status, "\n")
  invisible(x)
}
