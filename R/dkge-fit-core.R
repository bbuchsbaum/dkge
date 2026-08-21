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
#   scores_matrix      S x r    subject scores matrix (alias of s; used by tests)
#   v                  P x r    physical block loadings (NULL unless
#                               representation = block_biprojector)
#   s                  q x r    retained positive q-space factor
#   sdev               r        standard deviations per dimension
#   KU                 q x r    K %*% U (precomputed product)
#   representation     character  block_biprojector or qspace_moment
#   representation_reasons  character  why a fit is q-space-only (empty otherwise)
#
# INTERNAL (may change between versions):
#   Khalf              q x q    K^{1/2} factor
#   Kihalf             q x q    K^{-1/2} factor
#   Chat               q x q    raw compressed covariance accumulator
#   contribs           list     per-subject contribution matrices in the K-metric,
#                               unweighted (contribs[[s]] = A' M_s A). The subject
#                               weights are applied when pooling, so the identity
#                               is sum_s weights[s] * contribs[[s]] == Chat, and
#                               only when pooling is linear: missingness="none"
#                               and no effect weighting.
#   weights            numeric  per-subject weights
#   Braw               list[S]  input betas (q x P_s), zero-filled on unobserved rows
#   Omega              list[S]  per-subject AR/noise covariance structures
#   subjects           list[S]  slim subject records (debiasing sufficient
#                               statistics; beta/design/omega stripped)
#   provenance         list     data provenance metadata
#   kernel_info        list     kernel construction metadata
#   block_indices      list     multivarious block index mapping
#   X_concat           matrix   concatenated *training* blocks in K^{1/2} space
#                               (NULL unless keep_X=TRUE and the fit is an
#                               exact block biprojector). Chat equals
#                               X_concat %*% t(X_concat) only for that
#                               representation.
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
#   effect_scaling     character  "pooled_design" or "none"
#   effect_moment      q x q      pooled raw-effect-space moment (pre-K, pre-R)
#   effect_moments     list[S]    per-subject raw-effect-space moments (post-debias)
#   effect_moments_raw list[S]    moments before analytic noise subtraction. Only
#                                 debias="analytic" produces a distinct copy; for
#                                 every other setting this aliases effect_moments.
#                                 Under debias="split_half" that alias is the
#                                 debiased split cross-moment itself -- there is no
#                                 separate "raw" second moment on that path.
#   noise_moments      list[S]    subtracted noise moments (NULL unless debias="analytic")
#   effect_weight_spec list       dkge_effect_weights() specification used
#   effect_precision   list[S]    resolved per-effect precision vectors
#   effect_precision_diagnostics list  method-specific precision diagnostics
#                                 (DL tau2 / cap shares for random_effects)
#   pool_cache         list       sample-independent pair matrices reused when re-pooling
#   moment_diagnostics list       spectral diagnostics of the effect-space and
#                                 transformed moments
#   debias             character  finite-trial noise treatment used
#   pair_counts        matrix     sample-weighted effect-pair coverage mass used by
#                                 missingness (a plain count only when the sample
#                                 weights are all 1, as in a batch fit)
#   pair_weight        matrix     summed precision weight per effect pair
#   pair_ess           matrix     Kish effective sample size per effect pair
#   missingness        character  partial-effect missingness policy used at fit time
#   miss_args          list       missingness policy arguments
#   ridge_input        numeric    ridge penalty applied to input
#   rank_requested     integer    rank requested by caller
#   effective_rank     integer    effective rank after regularisation
#   rank_reduced       logical    whether rank was reduced from requested

#' Reconcile design-kernel labels with the dataset effect order
#'
#' `dkge_data()` fixes a global effect order; a kernel carrying its own
#' dimnames (as [design_kernel()] returns) must index the same effects. A set
#' mismatch is an error, and a permutation is repaired rather than silently
#' reindexing the kernel.
#'
#' @return A list with the reconciled kernel (`K`), the effect labels the fit
#'   should carry (`effects`), and `kernel_info` permuted to the same order
#'   whenever `K` is. The labels only change on the placeholder path, where
#'   the kernel's own (unique) names are more informative than
#'   `effect1...effectq`.
#' @keywords internal
#' @noRd
.dkge_align_kernel_effects <- function(K, effects, kernel_info = NULL) {
  out <- function(K, effects, kernel_info) {
    list(K = K, effects = effects, kernel_info = kernel_info)
  }
  if (is.null(effects)) return(out(K, effects, kernel_info))
  rn <- rownames(K)
  cn <- colnames(K)
  # Checked before any early return: a kernel whose rows and columns carry
  # different labels cannot be reconciled with anything, and letting it through
  # on the placeholder / duplicate-label paths allowed `fit$effects` to disagree
  # with `rownames(fit$K)`.
  if (!is.null(rn) && !is.null(cn) && !identical(as.character(rn), as.character(cn))) {
    stop("Design kernel row and column names must be identical.", call. = FALSE)
  }
  labels <- rn %||% cn
  if (is.null(labels)) return(out(K, effects, kernel_info))
  labels <- as.character(labels)
  effects <- as.character(effects)
  # `effect1...effectq` are the placeholders `dkge_subject()` invents when the
  # designs carry no column names: there is nothing to reconcile against, so the
  # kernel's own order stands. Its labels are the only real effect names in play,
  # so the fit adopts them -- unless they repeat, in which case they cannot
  # identify effects and are dropped rather than left to contradict
  # `fit$effects`.
  if (identical(effects, .default_effect_names(length(effects)))) {
    if (anyDuplicated(labels)) {
      dimnames(K) <- NULL
      return(out(K, effects, kernel_info))
    }
    dimnames(K) <- list(labels, labels)
    return(out(K, labels, kernel_info))
  }
  # Duplicated kernel labels cannot be matched by name; keep the kernel's own
  # order but say so, since a permuted kernel would then be used as-is.
  if (anyDuplicated(labels)) {
    if (!identical(labels, effects)) {
      warning(
        "Design kernel labels contain duplicates and cannot be matched to the ",
        "data effect names; the kernel is used in its own order.",
        call. = FALSE
      )
    }
    return(out(K, effects, kernel_info))
  }
  if (!setequal(labels, effects)) {
    missing_eff <- setdiff(effects, labels)
    extra_eff <- setdiff(labels, effects)
    stop(sprintf(
      paste0("Design kernel labels do not match the data effects.%s%s ",
             "Give the kernel dimnames matching the effect labels, ",
             "or drop them with `dimnames(K) <- NULL` if they are incidental."),
      if (length(missing_eff)) sprintf(" Missing from kernel: %s.", paste(missing_eff, collapse = ", ")) else "",
      if (length(extra_eff)) sprintf(" Unknown in kernel: %s.", paste(extra_eff, collapse = ", ")) else ""
    ), call. = FALSE)
  }
  if (!identical(labels, effects)) {
    idx <- match(effects, labels)
    K <- K[idx, idx, drop = FALSE]
    dimnames(K) <- list(effects, effects)
    kernel_info <- .dkge_permute_kernel_info(
      kernel_info,
      idx = idx,
      from_labels = labels,
      to_labels = effects
    )
    message("Design kernel reordered to match the data effect order.")
  } else {
    dimnames(K) <- list(effects, effects)
  }
  out(K, effects, kernel_info)
}

.dkge_kernel_info_space <- function(info, kernel_labels = NULL) {
  if (!is.list(info)) return(NULL)
  declared <- NULL
  if (is.list(info$coordinate_space)) {
    declared <- info$coordinate_space$kernel
  }
  declared <- declared %||% info$basis
  if (!is.null(declared)) {
    declared <- as.character(declared)
    if (length(declared) != 1L || !declared %in% c("cell", "effect")) {
      stop("kernel_info must declare its kernel coordinate space as 'cell' or 'effect'.",
           call. = FALSE)
    }
    return(declared)
  }

  # Backward-compatible derivation is name-based. Equal lengths alone are
  # deliberately not evidence of a shared coordinate system.
  kernel_labels <- as.character(kernel_labels %||% character())
  name_matches <- function(x) {
    x <- as.character(x %||% character())
    length(x) && !anyDuplicated(x) && !anyDuplicated(kernel_labels) &&
      setequal(x, kernel_labels)
  }
  cell_match <- name_matches(info$cell_labels)
  effect_match <- name_matches(info$effect_labels)
  if (cell_match && effect_match) {
    stop(
      "kernel_info coordinate space is ambiguous; declare `basis` or `coordinate_space$kernel`.",
      call. = FALSE
    )
  }
  if (cell_match) return("cell")
  if (effect_match) return("effect")
  NULL
}

.dkge_permute_kernel_info <- function(info, idx, from_labels = NULL,
                                      to_labels = NULL) {
  if (!is.list(info)) {
    return(info)
  }
  q <- length(idx)
  space <- .dkge_kernel_info_space(info, from_labels)
  declared <- if (is.list(info$coordinate_space)) info$coordinate_space else list()

  if (!is.null(info$kernel_labels)) {
    info$kernel_labels <- info$kernel_labels[idx]
  }
  if (identical(space, "cell") && !is.null(info$cell_labels)) {
    info$cell_labels <- info$cell_labels[idx]
  }
  if (identical(space, "effect") && !is.null(info$effect_labels)) {
    info$effect_labels <- info$effect_labels[idx]
  }
  if (identical(space, "cell") && !is.null(info$cells)) {
    cells <- as.data.frame(info$cells, stringsAsFactors = FALSE)
    info$cells <- cells[idx, , drop = FALSE]
    rownames(info$cells) <- NULL
  }
  if (!is.null(info$blocks)) {
    info$blocks <- lapply(info$blocks, function(b) {
      if (!is.numeric(b)) return(b)
      match(as.integer(b), idx)
    })
  }
  if (!is.null(info$map) && is.matrix(info$map)) {
    if (identical(space, "effect")) {
      info$map <- info$map[, idx, drop = FALSE]
    } else if (identical(space, "cell")) {
      info$map <- info$map[idx, , drop = FALSE]
    }
  }
  known <- c("cell_labels", "effect_labels", "kernel_labels", "cells",
             "blocks", "map", "levels", "coordinate_space",
             "factor_names", "term_names", "factor_scope", "term_scope",
             "block_factors", "basis", "dims")
  for (nm in setdiff(names(info), known)) {
    val <- info[[nm]]
    if (!is.atomic(val) || !is.null(dim(val))) next
    val_space <- declared[[nm]] %||% NULL
    named_on_kernel <- !is.null(names(val)) &&
      !anyDuplicated(names(val)) &&
      !is.null(from_labels) &&
      setequal(names(val), from_labels)
    if (named_on_kernel) {
      info[[nm]] <- val[match(to_labels, names(val))]
    } else if (!is.null(space) && identical(val_space, space) && length(val) == q) {
      info[[nm]] <- val[idx]
    }
  }
  info
}

.dkge_check_kernel_info <- function(fit) {
  info <- fit$kernel_info
  if (!is.list(info) || is.null(fit$effects)) {
    return(invisible(fit))
  }
  effects <- as.character(fit$effects)
  kernel_labels <- rownames(fit$K) %||% colnames(fit$K) %||% effects
  space <- .dkge_kernel_info_space(info, kernel_labels)
  if (is.null(space)) return(invisible(fit))

  if (!is.null(info$kernel_labels) &&
      !identical(as.character(info$kernel_labels), effects)) {
    stop("kernel_info$kernel_labels must match fit$effects after kernel alignment.",
         call. = FALSE)
  }

  cell_labels <- as.character(info$cell_labels %||% character())
  if (length(cell_labels) && anyDuplicated(cell_labels)) {
    stop("kernel_info$cell_labels must be unique.", call. = FALSE)
  }
  if (!is.null(info$cells) && nrow(as.data.frame(info$cells)) != length(cell_labels)) {
    stop("kernel_info$cells must have one row per cell label.", call. = FALSE)
  }

  if (identical(space, "cell") &&
      !identical(cell_labels, effects)) {
    stop("kernel_info$cell_labels must match fit$effects for a cell-basis kernel.",
         call. = FALSE)
  }

  if (identical(space, "effect")) {
    effect_labels <- as.character(info$effect_labels %||% character())
    if (!identical(effect_labels, effects)) {
      stop("kernel_info$effect_labels must match fit$effects for an effect-basis kernel.",
           call. = FALSE)
    }
    if (!is.null(info$map)) {
      map <- info$map
      if (!is.matrix(map) || nrow(map) != length(cell_labels) ||
          ncol(map) != length(effect_labels)) {
        stop("kernel_info$map must have cell rows and effect columns.", call. = FALSE)
      }
      if (!identical(as.character(rownames(map)), cell_labels)) {
        stop("kernel_info$map row names must match cell_labels.", call. = FALSE)
      }
      if (!identical(as.character(colnames(map)), effect_labels)) {
        stop("kernel_info$map column names must match effect_labels.", call. = FALSE)
      }
    }
  }

  if (!is.null(info$blocks)) {
    block_idx <- unlist(info$blocks, use.names = FALSE)
    if (length(block_idx) &&
        (!is.numeric(block_idx) || anyNA(block_idx) ||
         any(block_idx < 1L | block_idx > length(effects)))) {
      stop("kernel_info$blocks must index the declared kernel coordinate space.",
           call. = FALSE)
    }
  }
  invisible(fit)
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
  aligned <- .dkge_align_kernel_effects(K, effects, kernel_info)
  K <- aligned$K
  kernel_info <- aligned$kernel_info
  if (!identical(aligned$effects, effects)) {
    # The kernel supplied real effect labels where the data only had
    # `effect1...effectq` placeholders; carry them onto the fit and the dataset
    # so `fit$effects` and `rownames(fit$K)` stay the same vector.
    effects <- aligned$effects
    dataset$effects <- effects
  }
  K <- .dkge_validate_kernel(K)

  rank_requested <- if (is.null(rank)) {
    q
  } else {
    .dkge_validate_positive_integer(rank, "rank")
  }

  if (is.null(Omega_list)) {
    Omega_list <- vector("list", S)
  }
  Omega_list <- Map(function(om, B) .validate_omega(om, ncol(B)),
                    Omega_list, betas)
  dataset$omega <- Omega_list

  if (identical(effect_scaling, "pooled_design")) {
    ruler <- .dkge_compute_shared_ruler(designs)
    Btil <- .dkge_row_standardize(betas, ruler$R)
  } else {
    R_identity <- diag(1, q)
    if (!is.null(effects)) {
      dimnames(R_identity) <- list(effects, effects)
    }
    ruler <- list(R = R_identity, G_pool = R_identity)
    Btil <- betas
  }
  kernels <- .dkge_kernel_roots(K)
  if (kernels$rank == 0L) {
    .dkge_abort(
      "Kernel `K` has numerical rank zero; no DKGE component is estimable.",
      "dkge_zero_rank_error"
    )
  }
  rank <- min(rank_requested, q, kernels$rank)
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
  obs_masks <- .dkge_obs_masks_from_provenance(prepped$provenance,
                                               prepped$subject_ids,
                                               prepped$q)

  subject_weights <- .dkge_subject_weights(Btil, Omega_list, kernels$Khalf,
                                           w_method, w_tau,
                                           obs_masks = obs_masks)

  voxel_weights <- prepped$weight_eval$total
  voxel_weights_subject <- prepped$weight_eval$total_subject
  voxel_payload <- voxel_weights_subject %||% voxel_weights

  if (is.null(obs_masks)) {
    obs_masks <- replicate(prepped$S, rep(TRUE, prepped$q), simplify = FALSE)
  }
  effect_precision <- .dkge_resolve_effect_precision(
    prepped$dataset,
    prepped$effect_weight_spec,
    obs_masks = obs_masks
  )
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
    effect_precision_diagnostics = attr(effect_precision, "diagnostics"),
    effect_weight_spec = prepped$effect_weight_spec,
    effect_moment = accum$pooled,
    effect_moments = accum$moments,
    # Only analytic debiasing makes the raw moments differ from `moments`;
    # otherwise this aliases the same matrices instead of copying them.
    effect_moments_raw = accum$moments_raw %||% accum$moments,
    noise_moments = accum$noise_moments,
    pool_cache = accum$pool_cache,
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
    basis <- .dkge_select_k_basis(
      Chat,
      K = prepped$K,
      roots = prepped$kernels,
      requested_rank = rank,
      label = "Full-data DKGE"
    )
    eigChat <- basis$decomposition
    eig_vectors_full <- basis$vectors
    eig_values_full <- basis$values
    effective_rank <- basis$effective_rank
    rank_reduced <- basis$rank_reduced

    if (rank > basis$rank) {
      message(sprintf(
        "Requested rank %d exceeds effective rank %d. Reducing to %d components.",
        rank, effective_rank, basis$rank
      ))
    }
    rank <- basis$rank
    sdev <- sqrt(pmax(basis$kept_values, 0))
    U_hat <- basis$U_hat
    U <- basis$U

    if (!is.null(cpca_info)) {
      if (cpca_part %in% c("design", "both")) {
        cpca_info$U_design <- U
        cpca_info$evals_design <- eigChat$values
      }
      if (cpca_part == "resid") {
        cpca_info$U_resid <- U
        cpca_info$evals_resid <- eigChat$values
      } else if (cpca_part == "both") {
        residual_basis <- .dkge_select_k_basis(
          cpca_info$Chat_resid,
          K = prepped$K,
          roots = prepped$kernels,
          requested_rank = rank,
          label = "Residual CPCA"
        )
        cpca_info$evals_resid <- residual_basis$values
        cpca_info$U_resid <- residual_basis$U
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
      moment_rank = basis$moment_rank,
      effective_rank = effective_rank,
      rank_reduced = rank_reduced,
      moment_tolerance = basis$moment_tolerance,
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

  basis <- .dkge_select_k_basis(
    Chat,
    K = prepped$K,
    roots = prepped$kernels,
    requested_rank = rank,
    vectors = jd_res$Q,
    values = jd_res$diag_vals,
    label = "Full-data DKGE joint-diagonalization"
  )
  eig_vectors_full <- basis$vectors
  eig_values_full <- basis$values
  effective_rank <- basis$effective_rank
  rank_reduced <- basis$rank_reduced

  if (rank > basis$rank) {
    message(sprintf(
      "Requested rank %d exceeds effective rank %d. Reducing to %d components.",
      rank, effective_rank, basis$rank
    ))
  }
  rank <- basis$rank
  sdev <- sqrt(pmax(basis$kept_values, 0))
  U_hat <- basis$U_hat
  U <- basis$U

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
    moment_rank = basis$moment_rank,
    effective_rank = effective_rank,
    rank_reduced = rank_reduced,
    moment_tolerance = basis$moment_tolerance,
    cpca_info = cpca_info,
    solver = solver,
    jd = jd_res
  )
}

#' Classify the algebra represented by a fitted moment
#'
#' An exact `block_biprojector` is advertised only when Chat is the Gram of a
#' subject-by-voxel factor. Pair-normalized, missingness-transformed, debiased,
#' ridged, CPCA, and JD moments are `qspace_moment`.
#'
#' @keywords internal
#' @noRd
.dkge_fit_representation <- function(accum, solved, ridge,
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

#' Runtime check that an advertised block factor really factors Chat
#'
#' @keywords internal
#' @noRd
.dkge_block_factor_diagnostics <- function(fit, Xstar, tolerance = 1e-10) {
  r <- fit$rank
  chat_residual <- fit$Chat - tcrossprod(Xstar)
  chat_scale <- max(1, norm(fit$Chat, "F"))
  chat_error <- norm(chat_residual, "F") / chat_scale
  symmetry_error <- norm(fit$Chat - t(fit$Chat), "F") / chat_scale

  out <- list(
    chat_scale = chat_scale,
    factor_scale = norm(Xstar, "F"),
    chat_relative_error = chat_error,
    symmetry_relative_error = symmetry_error,
    orthogonality_error = 0,
    orthogonality_tolerance = tolerance,
    backward_error = 0,
    condition_proxy = 1
  )
  if (r == 0L) {
    return(out)
  }

  Q <- fit$eig_vectors_full[, seq_len(r), drop = FALSE]
  component_variances <- fit$sdev^2
  projected_residual <-
    crossprod(Q, fit$Chat %*% Q) - diag(component_variances, nrow = r)

  # With D = diag(sdev), V = Xstar' Q D^-1 and therefore
  #
  #   V'V - I = D^-1 Q'(Xstar Xstar' - Chat)Q D^-1
  #            + D^-1(Q' Chat Q - D^2)D^-1.
  #
  # The measured backward error below is consequently an upper bound, up to
  # ordinary rounding in the matrix products, on loss of orthogonality caused
  # by a small retained component. This is more informative than applying the
  # same absolute threshold to well- and ill-conditioned fits.
  min_variance <- min(component_variances)
  backward_error <-
    (norm(chat_residual, "2") + norm(projected_residual, "2")) /
    min_variance
  condition_proxy <- max(component_variances) / min_variance

  # Inverse square-root normalization can lose at most the square root of the
  # primary factorization tolerance. The error must also be explained by the
  # measured backward residual; poor conditioning alone never excuses it.
  orthogonality_tolerance <- max(
    tolerance,
    min(
      sqrt(tolerance),
      2 * backward_error + 32 * .Machine$double.eps
    )
  )

  out$orthogonality_error <-
    norm(crossprod(fit$v) - diag(r), "2")
  out$orthogonality_tolerance <- orthogonality_tolerance
  out$backward_error <- backward_error
  out$condition_proxy <- condition_proxy
  out
}

.dkge_validate_block_factor <- function(fit, Xstar, tolerance = 1e-10) {
  if (!identical(fit$representation, "block_biprojector")) {
    return(invisible(fit))
  }
  if (length(tolerance) != 1L || !is.numeric(tolerance) ||
      !is.finite(tolerance) || tolerance <= 0 || tolerance >= 1) {
    stop("Internal DKGE algebra error: block-factor tolerance must be a finite scalar in (0, 1).",
         call. = FALSE)
  }
  if (is.null(Xstar)) {
    stop("Internal DKGE algebra error: block representation has no Xstar factor.",
         call. = FALSE)
  }
  if (!is.matrix(fit$Chat) || !is.numeric(fit$Chat) ||
      nrow(fit$Chat) != ncol(fit$Chat) || any(!is.finite(fit$Chat))) {
    stop("Internal DKGE algebra error: Chat must be a finite square numeric matrix.",
         call. = FALSE)
  }
  if (!is.matrix(Xstar) || !is.numeric(Xstar) ||
      nrow(Xstar) != nrow(fit$Chat) || any(!is.finite(Xstar))) {
    stop("Internal DKGE algebra error: Xstar must be a finite numeric matrix with one row per Chat row.",
         call. = FALSE)
  }

  r <- fit$rank
  if (length(r) != 1L || !is.numeric(r) || !is.finite(r) ||
      r < 0 || r != as.integer(r)) {
    stop("Internal DKGE algebra error: block-factor rank must be one non-negative integer.",
         call. = FALSE)
  }
  r <- as.integer(r)
  if (!is.matrix(fit$v) || !is.numeric(fit$v) ||
      nrow(fit$v) != ncol(Xstar) || ncol(fit$v) != r ||
      any(!is.finite(fit$v))) {
    stop("Internal DKGE algebra error: block loadings must be a finite ncol(Xstar)-by-rank matrix.",
         call. = FALSE)
  }
  if (r > 0L) {
    if (!is.matrix(fit$eig_vectors_full) ||
        !is.numeric(fit$eig_vectors_full) ||
        nrow(fit$eig_vectors_full) != nrow(fit$Chat) ||
        ncol(fit$eig_vectors_full) < r ||
        any(!is.finite(fit$eig_vectors_full[, seq_len(r), drop = FALSE]))) {
      stop("Internal DKGE algebra error: fitted eigenvectors do not match Chat and rank.",
           call. = FALSE)
    }
    if (!is.numeric(fit$sdev) || length(fit$sdev) != r ||
        any(!is.finite(fit$sdev)) || any(fit$sdev <= 0)) {
      stop("Internal DKGE algebra error: block representation requires one positive finite sdev per retained component.",
           call. = FALSE)
    }
    component_variances <- fit$sdev^2
    if (any(!is.finite(component_variances)) ||
        any(component_variances <= 0)) {
      stop("Internal DKGE algebra error: block representation requires one positive finite sdev per retained component.",
           call. = FALSE)
    }
  }

  diagnostics <- .dkge_block_factor_diagnostics(fit, Xstar, tolerance)
  if (!is.finite(diagnostics$symmetry_relative_error) ||
      diagnostics$symmetry_relative_error > tolerance) {
    stop(sprintf(
      "Internal DKGE algebra error: Chat is not symmetric (relative error %.3e).",
      diagnostics$symmetry_relative_error
    ), call. = FALSE)
  }
  if (!is.finite(diagnostics$chat_relative_error) ||
      diagnostics$chat_relative_error > tolerance) {
    stop(sprintf(
      "Internal DKGE algebra error: Chat is not Xstar Xstar' (relative error %.3e).",
      diagnostics$chat_relative_error
    ), call. = FALSE)
  }
  if (!is.finite(diagnostics$orthogonality_error) ||
      !is.finite(diagnostics$orthogonality_tolerance) ||
      diagnostics$orthogonality_error > diagnostics$orthogonality_tolerance) {
    stop(sprintf(
      paste0(
        "Internal DKGE algebra error: advertised block loadings are not ",
        "orthonormal (spectral error %.3e; scale-aware bound %.3e; ",
        "backward error %.3e; condition proxy %.3e)."
      ),
      diagnostics$orthogonality_error,
      diagnostics$orthogonality_tolerance,
      diagnostics$backward_error,
      diagnostics$condition_proxy
    ), call. = FALSE)
  }
  invisible(fit)
}

#' Assemble the final dkge fit object
#'
#' Combines prepared payload, accumulation results, and eigen solution into the
#' object returned by `dkge_fit()`. Exact unregularized pooled moments inherit
#' from `multiblock_biprojector`; every other moment is a `dkge_qspace` object.
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

  representation <- .dkge_fit_representation(accum, solved, ridge)
  if (keep_X && !identical(representation$kind, "block_biprojector")) {
    stop(sprintf(
      paste0("`keep_X = TRUE` is unavailable for representation='qspace_moment' (%s). ",
             "The fitted moment has no exact subject-by-voxel block factor."),
      paste(representation$reasons, collapse = "; ")
    ), call. = FALSE)
  }

  total_clusters <- 0L
  block_indices <- vector("list", S)
  X_blocks <- vector("list", S)
  build_blocks <- identical(representation$kind, "block_biprojector")
  for (s in seq_len(S)) {
    Bts <- Btil[[s]]
    P_s <- ncol(Bts)
    idx <- seq_len(P_s) + total_clusters  # empty when P_s == 0 (avoids reversed `:`)
    block_indices[[s]] <- idx
    total_clusters <- total_clusters + P_s
    if (!build_blocks) {
      next
    }

    # Resolve to this subject's width exactly as `.dkge_build_moment_pool()`
    # does, so the training blocks stay consistent with Chat.
    w_s <- .dkge_subject_voxel_weights(
      accum$voxel_weights_subject %||% accum$voxel_weights,
      s, P_s, prepped$subject_ids[[s]]
    )
    Bw <- .dkge_scale_effect_columns(Bts, w_s)
    Omega <- Omega_list[[s]]
    block <- if (is.null(Omega)) {
      Bw
    } else if (is.vector(Omega)) {
      stopifnot(length(Omega) == ncol(Bw))
      sweep(Bw, 2L, sqrt(pmax(Omega, 0)), "*")
    } else {
      Omega <- as.matrix(Omega)
      stopifnot(nrow(Omega) == ncol(Bw), ncol(Omega) == ncol(Bw))
      Bw %*% sqrtm_sym(Omega)
    }
    # Rows this subject does not observe are already zero in *raw* effect
    # space, so `Btil = R' B` carries the masking through the pooled ruler and
    # the full `Khalf` congruence applies. Restricting to `Khalf[obs, obs]`
    # would break `Chat == X X'` and desync these blocks from
    # `dkge_transform_block()`.
    block <- kernels$Khalf %*% block
    # Rows are dropped from the effect labelling on purpose: after the `Khalf`
    # congruence each row is a K^{1/2}-weighted mixture of effects, not the
    # effect it would be named after. `.dkge_kernel_roots()` leaves the roots
    # unnamed for the same reason, so `dkge_transform_block()` -- which produces
    # this very matrix for new data -- returns unnamed rows; labelling the
    # training blocks here would desync the two. Column (parcel) labels are
    # genuine and are preserved through the Omega multiplication.
    rownames(block) <- NULL
    colnames(block) <- colnames(Bts)
    block <- sqrt(max(accum$subject_weights[s], 0)) * block
    X_blocks[[s]] <- block
  }

  scores <- if (rank > 0) {
    solved$U_hat %*% diag(solved$sdev, nrow = rank)
  } else {
    matrix(0, nrow = q, ncol = 0)
  }

  X_concat <- NULL
  V <- NULL
  multivarious_obj <- NULL
  if (build_blocks) {
    X_concat <- if (length(X_blocks)) do.call(cbind, X_blocks) else NULL
    if (rank > 0) {
      safe_sdev <- ifelse(solved$sdev > 0, solved$sdev, 1)
      V <- t(X_concat) %*% solved$U_hat %*% diag(1 / safe_sdev, nrow = rank)
      zero_cols <- which(solved$sdev <= 0)
      if (length(zero_cols) > 0) {
        V[, zero_cols] <- 0
      }
    } else {
      V <- matrix(0, nrow = total_clusters, ncol = 0)
    }
    x_for_preproc <- X_concat
    if (is.null(x_for_preproc)) {
      x_for_preproc <- matrix(numeric(0), nrow = q, ncol = 0)
    }
    preproc_obj <- multivarious::fit(multivarious::pass(), x_for_preproc)
    multivarious_obj <- multivarious::multiblock_biprojector(
      v = V,
      s = scores,
      sdev = solved$sdev,
      preproc = preproc_obj,
      block_indices = block_indices,
      classes = "dkge_core"
    )
  }

  X_store <- if (keep_X && build_blocks) X_concat else NULL

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
    pool_cache = accum$pool_cache,
    moment_diagnostics = accum$moment_diagnostics,
    weights = accum$subject_weights,
    Braw = dataset$betas,
    Btil = Btil,
    Omega = Omega_list,
    subject_ids = prepped$subject_ids,
    effects = prepped$effects,
    provenance = prepped$provenance,
    # Slim subject records: downstream re-pooling only needs the debiasing
    # sufficient statistics. `beta` lives in `$Braw` and the T_s x q designs are
    # not used after the pooled ruler is built.
    subjects = lapply(dataset$subjects, function(s) {
      s[setdiff(names(s), c("beta", "design", "omega"))]
    }),
    effect_scaling = prepped$effect_scaling,
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
    kernel_rank = kernels$rank,
    kernel_nullity = kernels$nullity,
    moment_rank = solved$moment_rank,
    effective_rank = solved$effective_rank,
    rank_reduced = rank < prepped$rank_requested,
    spectral_diagnostics = list(
      absolute_tolerance = kernels$absolute_tolerance,
      relative_tolerance = kernels$relative_tolerance,
      kernel_tolerance = kernels$tolerance,
      moment_tolerance = solved$moment_tolerance,
      kernel_rank = kernels$rank,
      kernel_nullity = kernels$nullity,
      moment_rank = solved$moment_rank,
      effective_rank = solved$effective_rank
    )
  )

  fit$representation <- representation$kind
  fit$representation_reasons <- representation$reasons
  fit$Chat_sym <- accum$Chat_sym
  fit$KU <- fit$K %*% fit$U
  if (fit$rank > 0L) {
    gram_error <- max(abs(crossprod(fit$U, fit$KU) - diag(fit$rank)))
    if (!is.finite(gram_error) || gram_error > 1e-7) {
      .dkge_abort(
        sprintf(
          "Full-data DKGE basis failed the K-orthonormal postcondition (maximum Gram error %.3e).",
          gram_error
        ),
        "dkge_kernel_geometry_error"
      )
    }
  }
  fit$scores_matrix <- fit$s
  .dkge_validate_block_factor(fit, X_concat)

  if (identical(representation$kind, "block_biprojector")) {
    for (nm in names(fit)) {
      multivarious_obj[[nm]] <- fit[[nm]]
    }
    class(multivarious_obj) <- unique(c("dkge", class(multivarious_obj)))
    .dkge_check_kernel_info(multivarious_obj)
    return(multivarious_obj)
  }

  class(fit) <- c("dkge_qspace", "dkge", "list")
  .dkge_check_kernel_info(fit)
  fit
}
