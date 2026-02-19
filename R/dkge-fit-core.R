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
      dataset$omega <- Map(function(om, B) .validate_omega(om, ncol(B)),
                           Omega_list, dataset$betas)
    }
  } else {
    dataset <- dkge_data(data, designs = designs, omega = Omega_list)
  }

  betas <- dataset$betas
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
  K <- .dkge_validate_kernel(K)

  rank_requested <- if (is.null(rank)) q else rank
  rank <- max(1L, min(rank_requested, q))

  if (is.null(Omega_list)) {
    Omega_list <- vector("list", S)
  }
  Omega_list <- Map(function(om, B) .validate_omega(om, ncol(B)),
                    Omega_list, betas)
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
    provenance = provenance,
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

    # Track effective rank (eigenvalues > 1e-12)
    effective_rank <- sum(eig_values_full > 1e-12)
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
    pos_idx <- eig_values > 1e-12
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

  # Track effective rank (eigenvalues > 1e-12)
  effective_rank <- sum(eig_values_full > 1e-12)
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
  pos_idx <- eig_values > 1e-12
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
    provenance = prepped$provenance,
    kernel_info = prepped$kernel_info,
    block_indices = block_indices,
    X_concat = X_store,
    V_full = solved$eig_vectors_full,
    eig_vectors_full = solved$eig_vectors_full,
    eig_values_full = solved$eig_values_full,
    rank = rank,
    cpca = solved$cpca_info,
    solver = solved$solver,
    jd = solved$jd,
    weight_spec = prepped$weight_spec,
    voxel_weights = accum$voxel_weights,
    voxel_weights_subject = accum$voxel_weights_subject,
    voxel_weights_prior = prepped$weight_eval$prior,
    voxel_weights_adapt = prepped$weight_eval$adapt,
    w_method = w_method,
    w_tau = w_tau,
    ridge_input = ridge,
    rank_requested = prepped$rank_requested,
    effective_rank = solved$effective_rank,
    rank_reduced = solved$rank_reduced
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
