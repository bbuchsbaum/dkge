# dkge-fit-core.R
# Internal staging functions that modularise the dkge_fit lifecycle.

# --- Fit Object Fields ---
# The assembled dkge object returned by dkge_fit() contains the following fields:
#
# STABLE (public contract - do not rename without version bump):
#   U                  q x r    K-orthonormal group latent basis
#   K                  q x q    design kernel (positive semidefinite)
#   R                  q x q    pooled design Cholesky factor
#   eig_vectors_full   q x q    full eigenvectors of pooled compressed covariance
#   eig_values_full    q        full eigenvalues (descending)
#   evals              q        full eigenvalues (descending, length = q; same as eig_values_full)
#   Btil               list[S]  compressed beta matrices (q x P_s) per subject
#   subject_ids        character  subject identifiers
#   effects            character  effect names
#   rank               integer    number of latent dimensions kept
#   Chat_sym           q x q    symmetric compressed covariance (used by tests)
#   scores_matrix      q x r    q-space spectral factor (alias of s)
#   v                  sum(P_s) x r or NULL  exact block loadings when available
#   s                  q x r    retained positive q-space spectral factor
#   sdev               r        standard deviations per dimension
#   KU                 q x r    K %*% U (precomputed product)
#
# INTERNAL (may change between versions):
#   Khalf              q x q    K^{1/2} factor
#   Kihalf             q x q    K^{-1/2} factor
#   Chat               q x q    raw compressed covariance accumulator
#   contribs           list     per-subject contribution matrices
#   weights            numeric  per-subject weights
#   Omega              list[S]  per-subject AR/noise covariance structures
#   provenance         list     data provenance metadata
#   kernel_info        list     kernel construction metadata
#   block_indices      list     multivarious block index mapping
#   X_concat           matrix   concatenated design matrices (NULL if keep_X=FALSE)
#   cpca               list     CPCA partitioning info (when cpca_part != "none")
#   solver             character  solver used
#   jd                 list     joint-diagonalisation info (when applicable)
#   weight_spec        list     weight specification used
#   voxel_weights      numeric  aggregated voxel weights
#   voxel_weights_subject  list  per-subject voxel weights
#   voxel_weights_prior    numeric  prior voxel weights
#   voxel_weights_adapt    numeric  adaptive voxel weights
#   w_method           character  weighting method
#   w_tau              numeric    weighting tau parameter
#   pair_counts        matrix     effect-pair coverage counts used by missingness
#   missingness        character  partial-effect missingness policy used at fit time
#   miss_args          list       missingness policy arguments
#   ridge_input        numeric    ridge penalty applied to input
#   rank_requested     integer    rank requested by caller
#   effective_rank     integer    effective rank after regularisation
#   rank_reduced       logical    whether rank was reduced from requested
#   representation     character  block_biprojector or qspace_moment
#   reconstruction     list       explicit reconstruction domain and target

#' Kernel roots for the fitting metric
#'
#' Unlike the historical helper used by post-fit alignment, this uses a
#' Moore-Penrose inverse on the numerical support of K. This is what makes the
#' fitted basis K-orthonormal for positive-semidefinite as well as positive-
#' definite kernels.
#'
#' @keywords internal
#' @noRd
.dkge_fit_kernel_roots <- function(K) {
  Ksym <- (K + t(K)) / 2
  eig <- eigen(Ksym, symmetric = TRUE)
  scale <- max(abs(eig$values))
  tol <- if (scale == 0) 0 else max(1e-12, nrow(K) * .Machine$double.eps) * scale
  support <- eig$values > tol
  vals <- ifelse(support, eig$values, 0)
  inv_sqrt <- numeric(length(vals))
  inv_sqrt[support] <- 1 / sqrt(vals[support])
  V <- eig$vectors
  Khalf <- V %*% diag(sqrt(vals), length(vals)) %*% t(V)
  Kihalf <- V %*% diag(inv_sqrt, length(vals)) %*% t(V)
  list(
    Khalf = (Khalf + t(Khalf)) / 2,
    Kihalf = (Kihalf + t(Kihalf)) / 2,
    eigen = eig,
    support = support,
    support_rank = sum(support),
    support_tolerance = tol
  )
}

#' Scale-relative threshold for retaining positive moment eigenpairs
#'
#' @keywords internal
#' @noRd
.dkge_positive_eigen_tolerance <- function(values) {
  scale <- max(abs(values))
  if (!is.finite(scale) || scale == 0) return(0)
  max(1e-12, length(values) * .Machine$double.eps) * scale
}

#' Enforce structural missingness before any effect-coordinate transformation
#'
#' @keywords internal
#' @noRd
.dkge_sanitize_unobserved_effects <- function(dataset) {
  masks <- .dkge_obs_masks_from_provenance(
    dataset$provenance %||% NULL,
    dataset$subject_ids,
    dataset$q
  )
  if (is.null(masks)) {
    masks <- replicate(dataset$n_subjects, rep(TRUE, dataset$q), simplify = FALSE)
  }

  for (s in seq_len(dataset$n_subjects)) {
    mask <- as.logical(masks[[s]])
    if (length(mask) != dataset$q) mask <- rep(TRUE, dataset$q)
    if (any(!mask)) {
      dataset$betas[[s]][!mask, ] <- 0
      dataset$designs[[s]][, !mask] <- 0
      if (!is.null(dataset$split_betas[[s]])) {
        dataset$split_betas[[s]] <- lapply(dataset$split_betas[[s]], function(B) {
          B[!mask, ] <- 0
          B
        })
      }
      if (!is.null(dataset$effect_noise_cov[[s]])) {
        dataset$effect_noise_cov[[s]][!mask, ] <- 0
        dataset$effect_noise_cov[[s]][, !mask] <- 0
      }
      if (!is.null(dataset$effect_information[[s]])) {
        dataset$effect_information[[s]][!mask, ] <- 0
        dataset$effect_information[[s]][, !mask] <- 0
      }
      if (!is.null(dataset$effect_score[[s]])) {
        dataset$effect_score[[s]][!mask, ] <- 0
      }
    }
    if (!is.null(dataset$subjects) && length(dataset$subjects) >= s) {
      dataset$subjects[[s]]$beta <- dataset$betas[[s]]
      dataset$subjects[[s]]$design <- dataset$designs[[s]]
      dataset$subjects[[s]]$split_betas <- dataset$split_betas[[s]]
      dataset$subjects[[s]]$effect_noise_cov <- dataset$effect_noise_cov[[s]]
      dataset$subjects[[s]]$effect_information <- dataset$effect_information[[s]]
      dataset$subjects[[s]]$effect_score <- dataset$effect_score[[s]]
      dataset$subjects[[s]]$observed_rows <- which(mask)
    }
  }
  list(dataset = dataset, obs_masks = masks)
}

#' Prepare DKGE inputs for fitting
#'
#' Harmonises data, resolves kernel metadata, and constructs the pooled
#' row-standardised betas used by downstream fitting stages.
#'
#' @keywords internal
#' @noRd
.dkge_fit_prepare <- function(data,
                              designs = NULL,
                              K = NULL,
                              Omega_list = NULL,
                              weights = NULL,
                              effect_weights = NULL,
                              rank = NULL,
                              effect_scaling = c("pooled_design", "none")) {
  effect_scaling <- match.arg(effect_scaling)
  if (inherits(data, "dkge_data")) {
    dataset <- data
    if (!is.null(Omega_list)) {
      stopifnot(length(Omega_list) == dataset$n_subjects)
      dataset$omega <- Map(function(om, B) .validate_omega(om, ncol(B)),
                           Omega_list, dataset$betas)
    }
  } else {
    dataset <- dkge_data(data, designs = designs, omega = Omega_list)
  }

  sanitized <- .dkge_sanitize_unobserved_effects(dataset)
  dataset <- sanitized$dataset
  obs_masks <- sanitized$obs_masks
  betas <- dataset$betas
  bad_finite <- which(vapply(betas, function(B) any(!is.finite(B)), logical(1)))
  if (length(bad_finite)) {
    labels <- (names(betas) %||% seq_along(betas))[bad_finite]
    stop(sprintf(
      paste0("Beta matrices contain non-finite values (NA/NaN/Inf) for subject(s): %s. ",
             "Clean or drop these voxels before fitting; .dkge_voxel_exclusion_mask() ",
             "identifies the offending columns."),
      paste(labels, collapse = ", ")
    ), call. = FALSE)
  }
  q_vals <- vapply(betas, nrow, integer(1))
  if (length(unique(q_vals)) > 1) {
    stop(sprintf(
      "All subjects must have the same number of design effects (q = nrow(B)). Found: %s",
      paste(names(betas) %||% seq_along(betas), q_vals, sep = "=", collapse = ", ")
    ), call. = FALSE)
  }
  designs <- dataset$designs
  Omega_list <- dataset$omega
  subject_ids <- dataset$subject_ids
  effects <- dataset$effects
  provenance <- dataset$provenance %||% NULL
  q <- dataset$q
  S <- dataset$n_subjects

  kernel_info <- NULL
  if (is.list(K) && !is.null(K$K)) {
    kernel_info <- if (!is.null(K$info)) K$info else NULL
    K <- K$K
  }

  stopifnot(is.matrix(K), nrow(K) == q, ncol(K) == q)
  kernel_rows <- rownames(K)
  kernel_cols <- colnames(K)
  labels_match <- !is.null(kernel_rows) && !is.null(kernel_cols) &&
    !anyDuplicated(kernel_rows) && !anyDuplicated(kernel_cols) &&
    setequal(kernel_rows, effects) && setequal(kernel_cols, effects)
  if (labels_match) {
    K <- K[effects, effects, drop = FALSE]
  }
  K <- .dkge_validate_kernel(K)

  rank_requested <- if (is.null(rank)) q else rank
  rank <- max(1L, min(rank_requested, q))

  if (is.null(Omega_list)) {
    Omega_list <- vector("list", S)
  }
  Omega_list <- Map(function(om, B) .validate_omega(om, ncol(B)),
                    Omega_list, betas)
  dataset$omega <- Omega_list

  if (identical(effect_scaling, "pooled_design")) {
    ruler <- .dkge_compute_shared_ruler(
      designs,
      information_list = dataset$effect_information %||% NULL
    )
    Btil <- .dkge_row_standardize(betas, ruler$R)
  } else {
    R_identity <- diag(1, q)
    if (!is.null(effects)) {
      dimnames(R_identity) <- list(effects, effects)
    }
    ruler <- list(R = R_identity, G_pool = R_identity)
    Btil <- betas
  }
  kernels <- .dkge_fit_kernel_roots(K)
  weight_spec <- if (is.null(weights)) dkge_weights(adapt = "none") else weights
  stopifnot(inherits(weight_spec, "dkge_weights"))
  effect_weight_spec <- effect_weights %||% dkge_effect_weights("none")
  if (!inherits(effect_weight_spec, "dkge_effect_weights")) {
    stop("`effect_weights` must be created by `dkge_effect_weights()`.",
         call. = FALSE)
  }

  kernel_payload <- .dkge_weight_kernel_payload(K, kernel_info)
  weight_eval <- .dkge_resolve_voxel_weights(weight_spec, Btil, kernel_payload)

  list(
    dataset = dataset,
    Btil = Btil,
    ruler = ruler,
    kernels = kernels,
    kernel_info = kernel_info,
    K = K,
    weight_spec = weight_spec,
    effect_weight_spec = effect_weight_spec,
    weight_eval = weight_eval,
    subject_ids = subject_ids,
    effects = effects,
    provenance = provenance,
    obs_masks = obs_masks,
    effect_scaling = effect_scaling,
    rank = rank,
    rank_requested = rank_requested,
    q = q,
    S = S
  )
}

#' Accumulate compressed covariance and weights
#'
#' Consumes the prepared DKGE payload and returns the compressed covariance in
#' the K-metric alongside subject and voxel weights.
#'
#' @keywords internal
#' @noRd
.dkge_fit_accumulate <- function(prepped,
                                 w_method,
                                 w_tau,
                                 missingness = c("none", "rescale", "mask", "shrink"),
                                 miss_args = list(),
                                 debias = c("none", "analytic", "split_half")) {
  missingness <- match.arg(missingness)
  debias <- match.arg(debias)
  Btil <- prepped$Btil
  Omega_list <- prepped$dataset$omega
  kernels <- prepped$kernels
  obs_masks <- prepped$obs_masks %||% .dkge_obs_masks_from_provenance(
    prepped$provenance, prepped$subject_ids, prepped$q
  )

  voxel_weights <- prepped$weight_eval$total
  voxel_weights_subject <- prepped$weight_eval$total_subject
  voxel_payload <- voxel_weights_subject %||% voxel_weights

  if (is.null(obs_masks)) {
    obs_masks <- replicate(prepped$S, rep(TRUE, prepped$q), simplify = FALSE)
  }
  subject_weights <- .dkge_subject_weights(
    prepped$dataset$betas,
    Omega_list,
    kernels$Khalf,
    w_method,
    w_tau,
    obs_masks = obs_masks,
    R = prepped$ruler$R,
    voxel_weights = voxel_payload
  )
  effect_precision <- .dkge_resolve_effect_precision(
    prepped$dataset,
    prepped$effect_weight_spec,
    obs_masks = obs_masks
  )
  effect_precision_diagnostics <- attr(effect_precision, "diagnostics")
  subjects <- prepped$dataset$subjects
  if (is.null(subjects)) {
    subjects <- lapply(seq_len(prepped$S), function(s) list(
      id = prepped$subject_ids[[s]],
      effect_noise_cov = prepped$dataset$effect_noise_cov[[s]],
      residual_variance = prepped$dataset$residual_variance[[s]],
      residual_df = prepped$dataset$residual_df[[s]],
      noise_trace = prepped$dataset$noise_trace[[s]]
    ))
  }

  accum <- .dkge_build_moment_pool(
    subjects = subjects,
    B_list = prepped$dataset$betas,
    Omega_list = Omega_list,
    voxel_weights = voxel_payload,
    obs_masks = obs_masks,
    subject_weights = subject_weights,
    effect_precision = effect_precision,
    effect_method = prepped$effect_weight_spec$method,
    R = prepped$ruler$R,
    Khalf = kernels$Khalf,
    missingness = missingness,
    miss_args = miss_args,
    debias = debias
  )
  Chat <- accum$Chat

  list(
    Chat = Chat,
    Chat_sym = (Chat + t(Chat)) / 2,
    contribs = accum$contribs,
    subject_weights = subject_weights,
    voxel_weights = voxel_weights,
    voxel_weights_subject = voxel_weights_subject,
    pair_counts = accum$pair_counts,
    pair_weight = accum$pair_weight,
    pair_ess = accum$pair_ess,
    effect_precision = effect_precision,
    effect_precision_diagnostics = effect_precision_diagnostics,
    effect_weight_spec = prepped$effect_weight_spec,
    effect_moment = accum$pooled,
    effect_moments = accum$moments,
    effect_moments_raw = accum$moments_raw,
    noise_moments = accum$noise_moments,
    moment_diagnostics = accum$diagnostics,
    debias = debias,
    missingness = missingness,
    miss_args = miss_args,
    obs_masks = obs_masks
  )
}

#' Solve the DKGE eigen problem
#'
#' Handles CPCA branches, ridge adjustments, and eigen decomposition in the
#' compressed K-metric space.
#'
#' @keywords internal
#' @noRd
.dkge_fit_solve <- function(prepped,
                            accum,
                            rank,
                            cpca_part,
                            cpca_blocks,
                            cpca_T,
                            cpca_ridge,
                            ridge,
                            solver = "pooled",
                            jd_control,
                            jd_mask,
                            jd_init) {
  solver <- match.arg(solver, c("pooled", "jd"))
  q <- prepped$q
  K <- prepped$K
  Chat <- accum$Chat
  contribs <- accum$contribs
  weights <- accum$subject_weights

  cpca_info <- NULL
  contribs_design <- NULL
  contribs_resid <- NULL
  if (cpca_part != "none") {
    if (!is.null(cpca_T)) {
      T_mat <- as.matrix(cpca_T)
      stopifnot(nrow(T_mat) == q)
    } else {
      stopifnot(!is.null(cpca_blocks), length(cpca_blocks) >= 1)
      T_mat <- diag(1, q)[, unique(cpca_blocks), drop = FALSE]
    }
    split <- dkge_cpca_split_chat(Chat, T_mat, K)
    cpca_info <- list(
      part = cpca_part,
      blocks = cpca_blocks,
      T = T_mat,
      ridge = cpca_ridge,
      P_hat = split$P_hat,
      Chat_design_raw = split$Chat_design,
      Chat_resid_raw = split$Chat_resid
    )
    if (cpca_part %in% c("design", "both")) {
      Chat_design <- split$Chat_design
      if (cpca_ridge > 0) Chat_design <- Chat_design + cpca_ridge * diag(q)
      Chat_design <- (Chat_design + t(Chat_design)) / 2
      cpca_info$Chat_design <- Chat_design
    }
    if (cpca_part %in% c("resid", "both")) {
      Chat_resid <- split$Chat_resid
      if (cpca_ridge > 0) Chat_resid <- Chat_resid + cpca_ridge * diag(q)
      Chat_resid <- (Chat_resid + t(Chat_resid)) / 2
      cpca_info$Chat_resid <- Chat_resid
    }
    if (cpca_part == "design") {
      Chat <- cpca_info$Chat_design
      contribs_design <- lapply(contribs, function(S) {
        M <- cpca_info$P_hat %*% S %*% cpca_info$P_hat
        (M + t(M)) / 2
      })
    } else if (cpca_part == "resid") {
      Chat <- cpca_info$Chat_resid
      Iq <- diag(1, q)
      P_resid <- Iq - cpca_info$P_hat
      contribs_resid <- lapply(contribs, function(S) {
        M <- P_resid %*% S %*% P_resid
        (M + t(M)) / 2
      })
    } else if (cpca_part == "both") {
      Chat <- cpca_info$Chat_design
      Iq <- diag(1, q)
      P_resid <- Iq - cpca_info$P_hat
      contribs_design <- lapply(contribs, function(S) {
        M <- cpca_info$P_hat %*% S %*% cpca_info$P_hat
        (M + t(M)) / 2
      })
      contribs_resid <- lapply(contribs, function(S) {
        M <- P_resid %*% S %*% P_resid
        (M + t(M)) / 2
      })
    }
  }

  if (ridge > 0) Chat <- Chat + ridge * diag(q)
  Chat <- (Chat + t(Chat)) / 2

  if (solver == "pooled") {
    eigChat <- eigen(Chat, symmetric = TRUE)
    eig_vectors_full <- eigChat$vectors
    eig_values_full <- eigChat$values

    # Retain only eigenpairs that are positive relative to the observed
    # spectral scale. There is deliberately no absolute floor: uniformly
    # rescaling a valid moment must not change its effective rank.
    eig_tol <- .dkge_positive_eigen_tolerance(eig_values_full)
    effective_rank <- sum(eig_values_full > eig_tol)
    rank_reduced <- FALSE

    # Warn if requested rank exceeds effective rank
    if (rank > effective_rank) {
      warning(sprintf(
        "Requested rank %d exceeds effective rank %d. Reducing to %d components.",
        rank, effective_rank, effective_rank
      ), call. = FALSE)
      rank <- effective_rank
      rank_reduced <- TRUE
    }

    eig_vectors <- eig_vectors_full[, seq_len(rank), drop = FALSE]
    eig_values <- eig_values_full[seq_len(rank)]
    pos_idx <- eig_values > eig_tol
    if (!all(pos_idx)) {
      eig_vectors <- eig_vectors[, pos_idx, drop = FALSE]
      eig_values <- eig_values[pos_idx]
      rank <- length(eig_values)
      rank_reduced <- TRUE
    }

    sdev <- sqrt(pmax(eig_values, 0))
    U_hat <- eig_vectors
    U <- prepped$kernels$Kihalf %*% U_hat

    if (!is.null(cpca_info)) {
      if (cpca_part %in% c("design", "both")) {
        cpca_info$U_design <- U
        cpca_info$evals_design <- eigChat$values
      }
      if (cpca_part == "resid") {
        cpca_info$U_resid <- U
        cpca_info$evals_resid <- eigChat$values
      } else if (cpca_part == "both") {
        eg_resid <- eigen(cpca_info$Chat_resid, symmetric = TRUE)
        cpca_info$evals_resid <- eg_resid$values
        cpca_info$U_resid <- prepped$kernels$Kihalf %*%
          eg_resid$vectors[, seq_len(rank), drop = FALSE]
      }
    }

    return(list(
      Chat = Chat,
      eig = eigChat,
      eig_vectors_full = eig_vectors_full,
      eig_values_full = eig_values_full,
      U_hat = U_hat,
      U = U,
      sdev = sdev,
      rank = rank,
      effective_rank = effective_rank,
      rank_reduced = rank_reduced,
      cpca_info = cpca_info,
      solver = solver,
      jd = NULL
    ))
  }

  # JD branch --------------------------------------------------------------

  A_list_solver <- switch(
    cpca_part,
    design = contribs_design,
    resid = contribs_resid,
    both = contribs_design,
    contribs
  )
  if (is.null(A_list_solver)) A_list_solver <- contribs

  mask_list <- NULL
  if (is.null(jd_mask)) {
    mask_list <- replicate(length(A_list_solver), NULL, simplify = FALSE)
  } else if (is.matrix(jd_mask)) {
    mask_list <- replicate(length(A_list_solver), jd_mask, simplify = FALSE)
  } else if (is.list(jd_mask)) {
    if (length(jd_mask) != length(A_list_solver)) {
      stop("jd_mask list must match the number of subject contributions.")
    }
    mask_list <- jd_mask
  } else {
    stop("jd_mask must be NULL, a matrix, or a list of matrices.")
  }

  Q_init <- NULL
  if (!is.null(jd_init)) {
    stopifnot(is.matrix(jd_init), nrow(jd_init) == q, ncol(jd_init) == q)
    Q_init <- .dkge_jd_retract(jd_init)
  }

  jd_res <- dkge_jd_solve(
    A_list = A_list_solver,
    weights = weights,
    rank = rank,
    Q_init = Q_init,
    mask_list = mask_list,
    Chat = Chat,
    control = jd_control
  )

  eig_vectors_full <- jd_res$Q
  eig_values_full <- jd_res$diag_vals

  # Scale-relative positivity tolerance; see the pooled branch above.
  eig_tol <- .dkge_positive_eigen_tolerance(eig_values_full)
  effective_rank <- sum(eig_values_full > eig_tol)
  rank_reduced <- FALSE

  # Warn if requested rank exceeds effective rank
  if (rank > effective_rank) {
    warning(sprintf(
      "Requested rank %d exceeds effective rank %d. Reducing to %d components.",
      rank, effective_rank, effective_rank
    ), call. = FALSE)
    rank <- effective_rank
    rank_reduced <- TRUE
  }

  eig_vectors <- eig_vectors_full[, seq_len(rank), drop = FALSE]
  eig_values <- eig_values_full[seq_len(rank)]
  pos_idx <- eig_values > eig_tol
  if (!all(pos_idx)) {
    eig_vectors <- eig_vectors[, pos_idx, drop = FALSE]
    eig_values <- eig_values[pos_idx]
    rank <- length(eig_values)
    rank_reduced <- TRUE
  }

  sdev <- sqrt(pmax(eig_values, 0))
  U_hat <- eig_vectors
  U <- prepped$kernels$Kihalf %*% U_hat

  if (!is.null(cpca_info)) {
    if (cpca_part %in% c("design", "both")) {
      cpca_info$U_design <- U
      cpca_info$evals_design <- eig_values_full
      cpca_info$jd_design <- jd_res
    }
    if (cpca_part == "resid") {
      cpca_info$U_resid <- U
      cpca_info$evals_resid <- eig_values_full
      cpca_info$jd_resid <- jd_res
    } else if (cpca_part == "both") {
      if (!is.null(contribs_resid)) {
        jd_resid <- dkge_jd_solve(
          A_list = contribs_resid,
          weights = weights,
          rank = rank,
          Q_init = NULL,
          mask_list = mask_list,
          Chat = cpca_info$Chat_resid,
          control = jd_control
        )
        cpca_info$evals_resid <- jd_resid$diag_vals
        cpca_info$U_resid <- prepped$kernels$Kihalf %*%
          jd_resid$Q[, seq_len(rank), drop = FALSE]
        cpca_info$jd_resid <- jd_resid
      } else {
        if (is.null(contribs_resid)) {
          warning(
            "cpca_part='both' requested but contribs_resid is NULL; U_resid will equal U_design.",
            call. = FALSE
          )
        }
        cpca_info$evals_resid <- eig_values_full
        cpca_info$U_resid <- U
      }
    }
  }

  list(
    Chat = Chat,
    eig = list(values = eig_values_full, vectors = eig_vectors_full),
    eig_vectors_full = eig_vectors_full,
    eig_values_full = eig_values_full,
    U_hat = U_hat,
    U = U,
    sdev = sdev,
    rank = rank,
    effective_rank = effective_rank,
    rank_reduced = rank_reduced,
    cpca_info = cpca_info,
    solver = solver,
    jd = jd_res
  )
}

#' Classify the algebra represented by a fitted moment
#'
#' @keywords internal
#' @noRd
.dkge_fit_representation <- function(prepped, accum, solved, ridge,
                                     force_qspace = FALSE) {
  reasons <- character(0)
  if (isTRUE(force_qspace)) {
    reasons <- c(reasons, "q-space representation was explicitly requested")
  }
  if (!identical(solved$solver, "pooled")) {
    reasons <- c(reasons, "solver is not the pooled symmetric eigensolver")
  }
  effect_method <- accum$effect_weight_spec$method %||% "none"
  if (!identical(effect_method, "none")) {
    reasons <- c(reasons, "effect-pair reliability normalization is active")
  }
  if (!identical(accum$missingness %||% "none", "none")) {
    reasons <- c(reasons, "a missingness transformation is active")
  }
  if (!identical(accum$debias %||% "none", "none")) {
    reasons <- c(reasons, "a debiased or cross-half moment is active")
  }
  if (!is.null(solved$cpca_info)) {
    reasons <- c(reasons, "CPCA filtering is active")
  }
  if (!isTRUE(all.equal(as.numeric(ridge), 0))) {
    reasons <- c(reasons, "a ridge term is present in Chat")
  }
  list(
    kind = if (length(reasons)) "qspace_moment" else "block_biprojector",
    reasons = reasons
  )
}

#' Build the exact physical block matrix for a factorizable DKGE moment
#'
#' Missing effect rows are zeroed in the raw effect coordinate before the
#' pooled ruler or kernel can mix rows.
#'
#' @keywords internal
#' @noRd
.dkge_build_xstar <- function(prepped, accum) {
  dataset <- prepped$dataset
  S <- dataset$n_subjects
  q <- prepped$q
  total_clusters <- 0L
  block_indices <- vector("list", S)
  X_blocks <- vector("list", S)

  for (s in seq_len(S)) {
    B <- as.matrix(dataset$betas[[s]])
    mask_s <- if (is.null(accum$obs_masks) || length(accum$obs_masks) < s) {
      NULL
    } else {
      accum$obs_masks[[s]]
    }
    P_s <- ncol(B)
    block_indices[[s]] <- seq_len(P_s) + total_clusters
    total_clusters <- total_clusters + P_s

    w_s <- if (is.list(accum$voxel_weights_subject)) {
      accum$voxel_weights_subject[[s]]
    } else {
      accum$voxel_weights_subject
    }
    if (is.null(w_s)) w_s <- accum$voxel_weights
    Omega <- dataset$omega[[s]]
    B <- .dkge_right_weighted_beta(B, Omega, w_s, mask_s)

    block <- prepped$kernels$Khalf %*% t(prepped$ruler$R) %*% B
    block <- sqrt(max(accum$subject_weights[s], 0)) * block
    X_blocks[[s]] <- block
  }

  list(
    Xstar = if (length(X_blocks)) do.call(cbind, X_blocks) else NULL,
    block_indices = block_indices,
    total_clusters = total_clusters
  )
}

#' Validate fit-object algebra before exposing it to downstream consumers
#'
#' @keywords internal
#' @noRd
.dkge_validate_fit_algebra <- function(fit, Xstar = NULL, tolerance = 1e-10) {
  r <- fit$rank
  if (r > 0L) {
    metric_error <- max(abs(crossprod(fit$U, fit$K %*% fit$U) - diag(r)))
    if (!is.finite(metric_error) || metric_error > tolerance) {
      stop(sprintf(
        "Internal DKGE algebra error: U is not K-orthonormal (max error %.3e).",
        metric_error
      ), call. = FALSE)
    }
    if (identical(fit$solver, "pooled")) {
      spectral_target <- fit$eig_vectors_full[, seq_len(r), drop = FALSE] %*%
        diag(fit$eig_values_full[seq_len(r)], r) %*%
        t(fit$eig_vectors_full[, seq_len(r), drop = FALSE])
      spectral_error <- norm(tcrossprod(fit$s) - spectral_target, "F") /
        max(1, norm(spectral_target, "F"))
      if (!is.finite(spectral_error) || spectral_error > tolerance) {
        stop(sprintf(
          "Internal DKGE algebra error: q-space scores do not reconstruct the retained spectrum (relative error %.3e).",
          spectral_error
        ), call. = FALSE)
      }
    }
  }

  if (identical(fit$representation, "block_biprojector")) {
    if (is.null(Xstar)) {
      stop("Internal DKGE algebra error: block representation has no Xstar factor.",
           call. = FALSE)
    }
    chat_error <- norm(fit$Chat - tcrossprod(Xstar), "F") /
      max(1, norm(fit$Chat, "F"))
    if (!is.finite(chat_error) || chat_error > tolerance) {
      stop(sprintf(
        "Internal DKGE algebra error: Chat is not Xstar Xstar' (relative error %.3e).",
        chat_error
      ), call. = FALSE)
    }
    if (r > 0L) {
      v_error <- max(abs(crossprod(fit$v) - diag(r)))
      if (!is.finite(v_error) || v_error > tolerance) {
        stop(sprintf(
          "Internal DKGE algebra error: advertised block loadings are not orthonormal (max error %.3e).",
          v_error
        ), call. = FALSE)
      }
    }
  }
  invisible(fit)
}

#' Assemble the final dkge fit object
#'
#' Combines prepared payload, accumulation results, and eigen solution into the
#' multiblock object returned by `dkge_fit()`.
#'
#' @keywords internal
#' @noRd
.dkge_fit_assemble <- function(prepped,
                               accum,
                               solved,
                               keep_X,
                               w_method,
                               w_tau,
                               ridge,
                               force_qspace = FALSE,
                               missingness = c("none", "rescale", "mask", "shrink"),
                               miss_args = list()) {
  missingness <- match.arg(missingness)
  dataset <- prepped$dataset
  S <- dataset$n_subjects
  Btil <- prepped$Btil
  kernels <- prepped$kernels
  Omega_list <- dataset$omega
  rank <- solved$rank
  q <- prepped$q

  representation <- .dkge_fit_representation(
    prepped, accum, solved, ridge, force_qspace = force_qspace
  )
  if (keep_X && !identical(representation$kind, "block_biprojector")) {
    stop(sprintf(
      paste0("`keep_X = TRUE` is unavailable for representation='qspace_moment' (%s). ",
             "The fitted moment has no exact subject-by-voxel block factor."),
      paste(representation$reasons, collapse = "; ")
    ), call. = FALSE)
  }

  if (identical(representation$kind, "block_biprojector")) {
    block_layout <- .dkge_build_xstar(prepped, accum)
    X_concat <- block_layout$Xstar
    block_indices <- block_layout$block_indices
    total_clusters <- block_layout$total_clusters
  } else {
    widths <- vapply(Btil, ncol, integer(1))
    ends <- cumsum(widths)
    starts <- ends - widths + 1L
    block_indices <- Map(
      function(from, to, width) if (width > 0L) seq.int(from, to) else integer(0),
      starts, ends, widths
    )
    total_clusters <- sum(widths)
    X_concat <- NULL
  }

  scores <- if (rank > 0L) {
    solved$U_hat %*% diag(solved$sdev, nrow = rank)
  } else {
    matrix(0, nrow = q, ncol = 0)
  }
  V <- if (identical(representation$kind, "block_biprojector")) {
    if (rank > 0L) {
      t(X_concat) %*% solved$U_hat %*%
        diag(1 / solved$sdev, nrow = rank)
    } else {
      matrix(0, nrow = total_clusters, ncol = 0)
    }
  } else {
    NULL
  }

  multivarious_obj <- NULL
  if (identical(representation$kind, "block_biprojector")) {
    preproc_obj <- multivarious::fit(multivarious::pass(), X_concat)
    multivarious_obj <- multivarious::multiblock_biprojector(
      v = V,
      s = scores,
      sdev = solved$sdev,
      preproc = preproc_obj,
      block_indices = block_indices,
      classes = "dkge_core"
    )
  }

  X_store <- if (keep_X && identical(representation$kind, "block_biprojector")) {
    X_concat
  } else {
    NULL
  }

  fit_provenance <- prepped$provenance
  error_models <- lapply(dataset$subjects, `[[`, "error_model")
  names(error_models) <- prepped$subject_ids
  fit_provenance$error_models <- error_models

  fit <- list(
    v = V,
    s = scores,
    sdev = solved$sdev,
    U = solved$U,
    evals = solved$eig$values,
    R = prepped$ruler$R,
    K = prepped$K,
    Khalf = kernels$Khalf,
    Kihalf = kernels$Kihalf,
    Chat = solved$Chat,
    contribs = accum$contribs,
    effect_moment = accum$effect_moment,
    effect_moments = accum$effect_moments,
    effect_moments_raw = accum$effect_moments_raw,
    noise_moments = accum$noise_moments,
    moment_diagnostics = accum$moment_diagnostics,
    weights = accum$subject_weights,
    Braw = dataset$betas,
    Btil = Btil,
    Omega = Omega_list,
    subject_ids = prepped$subject_ids,
    effects = prepped$effects,
    provenance = fit_provenance,
    error_models = error_models,
    subjects = dataset$subjects,
    effect_scaling = prepped$effect_scaling,
    ruler_source = prepped$ruler$source %||% "identity",
    kernel_info = prepped$kernel_info,
    block_indices = block_indices,
    X_concat = X_store,
    eig_vectors_full = solved$eig_vectors_full,
    eig_values_full = solved$eig_values_full,
    rank = rank,
    cpca = solved$cpca_info,
    solver = solved$solver,
    jd = solved$jd,
    weight_spec = prepped$weight_spec,
    effect_weight_spec = accum$effect_weight_spec,
    effect_precision = accum$effect_precision,
    effect_precision_diagnostics = accum$effect_precision_diagnostics,
    voxel_weights = accum$voxel_weights,
    voxel_weights_subject = accum$voxel_weights_subject,
    voxel_weights_prior = prepped$weight_eval$prior,
    voxel_weights_adapt = prepped$weight_eval$adapt,
    w_method = w_method,
    w_tau = w_tau,
    pair_counts = accum$pair_counts,
    pair_weight = accum$pair_weight,
    pair_ess = accum$pair_ess,
    debias = accum$debias,
    missingness = missingness,
    miss_args = miss_args,
    ridge_input = ridge,
    rank_requested = prepped$rank_requested,
    effective_rank = solved$effective_rank,
    rank_reduced = solved$rank_reduced
  )

  fit$representation <- representation$kind
  fit$representation_reasons <- representation$reasons
  fit$moment_estimator <- .dkge_moment_estimator_name(fit)
  fit$reconstruction <- if (identical(representation$kind, "block_biprojector")) {
    list(domain = "block_matrix", target = "rank_r_svd_of_Xstar")
  } else if (identical(solved$solver, "pooled")) {
    list(domain = "q_space", target = "retained_positive_spectrum_of_Chat")
  } else {
    list(domain = "q_space", target = "solver_specific_joint_diagonalization")
  }

  fit$Chat_sym <- solved$Chat
  fit$KU <- fit$K %*% fit$U

  fit$scores_matrix <- fit$s

  .dkge_validate_fit_algebra(fit, Xstar = X_concat)

  if (identical(representation$kind, "block_biprojector")) {
    for (nm in names(fit)) {
      multivarious_obj[[nm]] <- fit[[nm]]
    }
    class(multivarious_obj) <- unique(c("dkge", class(multivarious_obj)))
    return(multivarious_obj)
  }
  class(fit) <- c("dkge_qspace", "dkge", "list")
  fit
}
