# dkge-fit-core.R
# Internal staging functions that modularise the dkge_fit lifecycle.

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
                              rank = NULL) {
  if (inherits(data, "dkge_data")) {
    dataset <- data
    if (!is.null(Omega_list)) {
      stopifnot(length(Omega_list) == dataset$n_subjects)
      dataset$omega <- Omega_list
    }
  } else {
    dataset <- dkge_data(data, designs = designs, omega = Omega_list)
  }

  betas <- dataset$betas
  designs <- dataset$designs
  Omega_list <- dataset$omega
  subject_ids <- dataset$subject_ids
  effects <- dataset$effects
  q <- dataset$q
  S <- dataset$n_subjects

  kernel_info <- NULL
  if (is.list(K) && !is.null(K$K)) {
    kernel_info <- if (!is.null(K$info)) K$info else NULL
    K <- K$K
  }

  stopifnot(is.matrix(K), nrow(K) == q, ncol(K) == q)

  rank_requested <- if (is.null(rank)) min(q, 10L) else rank
  rank <- max(1L, min(rank_requested, q))

  if (is.null(Omega_list)) {
    Omega_list <- vector("list", S)
  }
  dataset$omega <- Omega_list

  ruler <- .dkge_compute_shared_ruler(designs)
  Btil <- .dkge_row_standardize(betas, ruler$R)
  kernels <- .dkge_kernel_roots(K)
  weight_spec <- if (is.null(weights)) dkge_weights(adapt = "none") else weights
  stopifnot(inherits(weight_spec, "dkge_weights"))

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
    weight_eval = weight_eval,
    subject_ids = subject_ids,
    effects = effects,
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
                                 w_tau) {
  Btil <- prepped$Btil
  Omega_list <- prepped$dataset$omega
  kernels <- prepped$kernels

  subject_weights <- .dkge_subject_weights(Btil, Omega_list, kernels$Khalf,
                                           w_method, w_tau)

  voxel_weights <- prepped$weight_eval$total
  voxel_weights_subject <- prepped$weight_eval$total_subject
  voxel_payload <- voxel_weights_subject %||% voxel_weights

  accum <- .dkge_accumulate_chat(Btil, Omega_list, kernels$Khalf,
                                 subject_weights,
                                 voxel_weights = voxel_payload)

  list(
    Chat = accum$Chat,
    Chat_sym = accum$Chat,
    contribs = accum$contribs,
    subject_weights = subject_weights,
    voxel_weights = voxel_weights,
    voxel_weights_subject = voxel_weights_subject
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
                            ridge) {
  q <- prepped$q
  K <- prepped$K
  Chat <- accum$Chat

  cpca_info <- NULL
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
    } else if (cpca_part == "resid") {
      Chat <- cpca_info$Chat_resid
    } else if (cpca_part == "both") {
      Chat <- cpca_info$Chat_design
    }
  }

  if (ridge > 0) Chat <- Chat + ridge * diag(q)

  eigChat <- eigen((Chat + t(Chat)) / 2, symmetric = TRUE)
  eig_vectors_full <- eigChat$vectors
  eig_values_full <- eigChat$values

  eig_vectors <- eig_vectors_full[, seq_len(rank), drop = FALSE]
  eig_values <- eig_values_full[seq_len(rank)]
  pos_idx <- eig_values > 1e-12
  if (!all(pos_idx)) {
    eig_vectors <- eig_vectors[, pos_idx, drop = FALSE]
    eig_values <- eig_values[pos_idx]
    rank <- length(eig_values)
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
      eg_resid <- eigen((cpca_info$Chat_resid + t(cpca_info$Chat_resid)) / 2,
                        symmetric = TRUE)
      cpca_info$evals_resid <- eg_resid$values
      cpca_info$U_resid <- prepped$kernels$Kihalf %*%
        eg_resid$vectors[, seq_len(rank), drop = FALSE]
    }
  }

  list(
    Chat = Chat,
    eig = eigChat,
    eig_vectors_full = eig_vectors_full,
    eig_values_full = eig_values_full,
    U_hat = U_hat,
    U = U,
    sdev = sdev,
    rank = rank,
    cpca_info = cpca_info
  )
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
                               ridge) {
  dataset <- prepped$dataset
  S <- dataset$n_subjects
  Btil <- prepped$Btil
  kernels <- prepped$kernels
  Omega_list <- dataset$omega
  rank <- solved$rank
  q <- prepped$q

  total_clusters <- 0L
  block_indices <- vector("list", S)
  X_blocks <- vector("list", S)
  for (s in seq_len(S)) {
    Bts <- Btil[[s]]
    P_s <- ncol(Bts)
    idx <- (total_clusters + 1L):(total_clusters + P_s)
    block_indices[[s]] <- idx
    total_clusters <- total_clusters + P_s

    w_s <- if (is.list(accum$voxel_weights_subject)) {
      accum$voxel_weights_subject[[s]]
    } else {
      accum$voxel_weights_subject
    }
    if (is.null(w_s)) w_s <- accum$voxel_weights
    if (!is.null(w_s) && length(w_s) != ncol(Bts)) {
      w_s <- rep(w_s, length.out = ncol(Bts))
    }
    Bw <- if (is.null(w_s) || length(w_s) == 0L) {
      Bts
    } else {
      sweep(Bts, 2L, sqrt(pmax(w_s, 0)), "*")
    }
    Omega <- Omega_list[[s]]
    if (is.null(Omega)) {
      block <- Bw
    } else if (is.vector(Omega)) {
      stopifnot(length(Omega) == ncol(Bw))
      block <- sweep(Bw, 2L, sqrt(pmax(Omega, 0)), "*")
    } else {
      Omega <- as.matrix(Omega)
      stopifnot(nrow(Omega) == ncol(Bw), ncol(Omega) == ncol(Bw))
      block <- Bw %*% sqrtm_sym(Omega)
    }
    block <- kernels$Khalf %*% block
    block <- sqrt(max(accum$subject_weights[s], 0)) * block
    X_blocks[[s]] <- block
  }

  X_concat <- if (length(X_blocks)) do.call(cbind, X_blocks) else NULL
  if (rank > 0) {
    safe_sdev <- ifelse(solved$sdev > 0, solved$sdev, 1)
    V <- t(X_concat) %*% solved$U_hat %*% diag(1 / safe_sdev, nrow = rank)
    zero_cols <- which(solved$sdev <= 0)
    if (length(zero_cols) > 0) {
      V[, zero_cols] <- 0
    }
    scores <- solved$U_hat %*% diag(solved$sdev, nrow = rank)
  } else {
    V <- matrix(0, nrow = total_clusters, ncol = 0)
    scores <- matrix(0, nrow = q, ncol = 0)
  }

  preproc_obj <- multivarious::prep(multivarious::pass())
  multivarious_obj <- multivarious::multiblock_biprojector(
    v = V,
    s = scores,
    sdev = solved$sdev,
    preproc = preproc_obj,
    block_indices = block_indices,
    classes = "dkge_core"
  )

  X_store <- if (keep_X) X_concat else NULL

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
    weights = accum$subject_weights,
    Btil = Btil,
    Omega = Omega_list,
    subject_ids = prepped$subject_ids,
    effects = prepped$effects,
    kernel_info = prepped$kernel_info,
    block_indices = block_indices,
    X_concat = X_store,
    V_full = solved$eig_vectors_full,
    eig_vectors_full = solved$eig_vectors_full,
    eig_values_full = solved$eig_values_full,
    rank = rank,
    cpca = solved$cpca_info,
    weight_spec = prepped$weight_spec,
    voxel_weights = accum$voxel_weights,
    voxel_weights_subject = accum$voxel_weights_subject,
    voxel_weights_prior = prepped$weight_eval$prior,
    voxel_weights_adapt = prepped$weight_eval$adapt,
    w_method = w_method,
    w_tau = w_tau,
    ridge_input = ridge,
    rank_requested = prepped$rank_requested
  )

  fit$Chat_sym <- accum$Chat_sym
  fit$KU <- fit$K %*% fit$U

  fit$variables <- fit$v
  fit$scores_matrix <- fit$s

  for (nm in names(fit)) {
    multivarious_obj[[nm]] <- fit[[nm]]
  }
  class(multivarious_obj) <- unique(c("dkge", class(multivarious_obj)))
  multivarious_obj
}
