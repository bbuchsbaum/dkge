#' DKGE classifier information maps
#'
#' Functions for translating latent-space classifier weights into subject- and
#' group-level spatial maps through DKGE's transport and rendering stack.
#'
#' @name dkge-info-maps
NULL

#' Decoder-style information map from latent classifier weights
#'
#' Pulls a cross-fitted latent weight vector \eqn{\beta^{(-s)}} back to each
#' subject's cluster space via \eqn{m_s = A_s \beta^{(-s)}}, transports those
#' maps to anchors, aggregates across subjects, and optionally performs
#' inference.
#'
#' @param fit Fitted `dkge` object.
#' @param betas Either a numeric vector of length `r` or a list of length `S`
#'   (number of subjects) containing one latent weight vector per subject.
#' @param renderer Renderer produced by [dkge_build_renderer()].
#' @param lambda Non-negative smoothing parameter applied via
#'   [dkge_anchor_aggregate()].
#' @param to_vox Logical; when `TRUE` and a decoder is available in `renderer`,
#'   a dense voxel map is produced.
#' @param inference One of `"none"`, `"signflip"`, or `"parametric"` specifying the
#'   group inference routine applied to anchor maps.
#' @return A list containing the mean anchor field, optional voxel map,
#'   per-subject anchor maps, and inference outputs (t- and p-values when
#'   requested).
#' @export
#' @examples
#' \dontrun{
#' info <- dkge_info_map_from_classifier(fit, clf$beta_by_subject, renderer)
#' }
dkge_info_map_from_classifier <- function(fit,
                                           betas,
                                           renderer,
                                           lambda = 0,
                                           to_vox = TRUE,
                                           inference = c("none", "signflip", "parametric")) {
  inference <- match.arg(inference)
  stopifnot(inherits(fit, "dkge"))
  stopifnot(is.list(renderer), !is.null(renderer$mapper_fits))

  A_list <- dkge_cluster_loadings(fit)
  S <- length(A_list)
  r <- ncol(A_list[[1]])

  beta_list <- .dkge_normalise_betas(betas, S, r)

  subj_cluster <- lapply(seq_len(S), function(s) {
    as.numeric(A_list[[s]] %*% beta_list[[s]])
  })

  subj_anchor <- .dkge_map_subject_clusters(renderer, subj_cluster)
  agg <- .dkge_aggregate_anchor_maps(renderer, subj_anchor, lambda)

  vox <- NULL
  if (to_vox && !is.null(renderer$decoder)) {
    vox <- dkge_anchor_to_voxel_apply(renderer$decoder, agg$field)
  }

  infer <- .dkge_anchor_inference(subj_anchor, method = inference)

  list(
    mean_anchor = agg$field,
    mean_voxel = vox,
    subj_anchor = subj_anchor,
    t_anchor = infer$stat,
    p_anchor = infer$p,
    aggregate = agg$details,
    renderer = renderer,
    meta = list(kind = "decoder", lambda = lambda, inference = inference)
  )
}

#' Haufe-style encoding maps from latent classifiers
#'
#' Converts discriminative weights into encoding (activation) patterns in latent
#' space using the Haufe transform \eqn{a_z = \Sigma_{zz} \beta}, then projects
#' them to subject cluster space and through the renderer pipeline.
#'
#' @param clf Cross-fitted classifier returned by
#'   [dkge_cv_train_latent_classifier()].
#' @param Z_by_subject Optional list of latent cluster features used to estimate
#'   fold covariances when they are not stored in `clf`.
#' @inheritParams dkge_info_map_from_classifier
#' @return List mirroring [dkge_info_map_from_classifier()] with
#'   `meta$kind = "haufe"`.
#' @export
#' @examples
#' \dontrun{
#' enc <- dkge_info_map_haufe(fit, clf, renderer, inference = "signflip")
#' }
dkge_info_map_haufe <- function(fit,
                                clf,
                                renderer,
                                Z_by_subject = NULL,
                                lambda = 0,
                                to_vox = TRUE,
                                inference = c("none", "signflip", "parametric")) {
  inference <- match.arg(inference)
  stopifnot(inherits(clf, "dkge_clf"))

  if (is.null(Z_by_subject)) {
    Z_by_subject <- dkge_project_clusters_to_latent(fit)
  }
  stopifnot(is.list(Z_by_subject), length(Z_by_subject) == length(fit$Btil))
  Z_by_subject <- lapply(Z_by_subject, as.matrix)

  A_list <- dkge_cluster_loadings(fit)
  S <- length(A_list)
  r <- ncol(A_list[[1]])

  beta_list <- .dkge_normalise_betas(clf$beta_by_subject, S, r)
  Sigma_by_fold <- .dkge_haufe_fold_covariances(fit, clf, Z_by_subject)

  subj_cluster <- vector("list", S)
  for (s in seq_len(S)) {
    fold_id <- clf$fold_assignment[s]
    beta <- beta_list[[s]]
    Sigma <- Sigma_by_fold[[fold_id]]
    a_z <- as.numeric(Sigma %*% beta)
    subj_cluster[[s]] <- as.numeric(A_list[[s]] %*% a_z)
  }

  subj_anchor <- .dkge_map_subject_clusters(renderer, subj_cluster)
  agg <- .dkge_aggregate_anchor_maps(renderer, subj_anchor, lambda)

  vox <- NULL
  if (to_vox && !is.null(renderer$decoder)) {
    vox <- dkge_anchor_to_voxel_apply(renderer$decoder, agg$field)
  }
  infer <- .dkge_anchor_inference(subj_anchor, method = inference)

  list(
    mean_anchor = agg$field,
    mean_voxel = vox,
    subj_anchor = subj_anchor,
    t_anchor = infer$stat,
    p_anchor = infer$p,
    aggregate = agg$details,
    renderer = renderer,
    meta = list(kind = "haufe", lambda = lambda, inference = inference)
  )
}

#' Group-LOCO anchor importance (zeroing proxy)
#'
#' Approximates the influence of each anchor by zeroing its neighbourhood in the
#' subject maps derived from classifier weights and accumulating the lost margin
#' magnitude across subjects.
#'
#' @param clf Cross-fitted classifier from [dkge_cv_train_latent_classifier()].
#' @param neighborhoods Optional list of integer vectors defining anchor
#'   neighbourhoods. When `NULL`, the function derives them from the renderer's
#'   anchor graph.
#' @param k_nn When neighbourhoods are not supplied, limits the automatically
#'   derived neighbourhood size.
#' @param aggregate Aggregation method across subjects: `"mean"` or `"sum"`.
#' @inheritParams dkge_info_map_from_classifier
#' @return A list with LOCO scores per anchor, per-subject anchor maps, and the
#'   neighbourhood definition.
#' @export
dkge_info_map_loco <- function(fit,
                               clf,
                               renderer,
                               neighborhoods = NULL,
                               k_nn = 32,
                               aggregate = c("mean", "sum")) {
  aggregate <- match.arg(aggregate)
  stopifnot(inherits(clf, "dkge_clf"))

  A_list <- dkge_cluster_loadings(fit)
  S <- length(A_list)
  r <- ncol(A_list[[1]])

  beta_list <- .dkge_normalise_betas(clf$beta_by_subject, S, r)
  subj_cluster <- lapply(seq_len(S), function(s) {
    as.numeric(A_list[[s]] %*% beta_list[[s]])
  })
  subj_anchor <- .dkge_map_subject_clusters(renderer, subj_cluster)
  anchor_mat <- do.call(rbind, subj_anchor)

  neighborhoods <- neighborhoods %||% .dkge_default_anchor_neighborhoods(renderer, k_nn)

  loco <- numeric(length(neighborhoods))
  for (j in seq_along(neighborhoods)) {
    idx <- neighborhoods[[j]]
    removed <- rowSums(anchor_mat[, idx, drop = FALSE])
    loco[j] <- if (aggregate == "mean") mean(abs(removed)) else sum(abs(removed))
  }

  list(
    loco_anchor = loco,
    subj_anchor = subj_anchor,
    neighborhoods = neighborhoods,
    renderer = renderer,
    meta = list(kind = "loco_zero", aggregate = aggregate, k_nn = k_nn)
  )
}

# -----------------------------------------------------------------------------
# Internal helpers ------------------------------------------------------------

#' @keywords internal
#' @noRd
.dkge_normalise_betas <- function(betas, S, r) {
  if (is.list(betas)) {
    stopifnot(length(betas) == S)
    beta_list <- lapply(betas, function(b) {
      stopifnot(is.numeric(b), length(b) == r)
      as.numeric(b)
    })
  } else if (is.numeric(betas) && length(betas) == r) {
    beta_list <- replicate(S, as.numeric(betas), simplify = FALSE)
  } else {
    stop("'betas' must be a length-r numeric vector or a list of length S of such vectors.")
  }
  beta_list
}

#' @keywords internal
#' @noRd
.dkge_map_subject_clusters <- function(renderer, subj_cluster) {
  mapper_fits <- renderer$mapper_fits
  stopifnot(length(mapper_fits) == length(subj_cluster))
  lapply(seq_along(subj_cluster), function(s) {
    apply_mapper(mapper_fits[[s]], subj_cluster[[s]])
  })
}

#' @keywords internal
#' @noRd
.dkge_aggregate_anchor_maps <- function(renderer, subj_anchor, lambda) {
  weights <- renderer$weights %||% rep(1, length(subj_anchor))
  graph <- renderer$graph
  L <- if (!is.null(graph)) graph$L else NULL
  agg <- dkge_anchor_aggregate(subj_anchor,
                               subj_weights = weights,
                               L = L,
                               lambda = lambda)
  list(field = agg$y, details = agg)
}

#' @keywords internal
#' @noRd
.dkge_anchor_inference <- function(subj_anchor, method = c("none", "signflip", "parametric"), B = 2000) {
  method <- match.arg(method)
  if (identical(method, "none")) {
    return(list(stat = NULL, p = NULL))
  }
  Y <- do.call(rbind, subj_anchor)
  if (identical(method, "signflip")) {
    res <- dkge_signflip_maxT(Y, B = B)
    return(list(stat = res$stat, p = res$p))
  }
  mu <- colMeans(Y)
  se <- apply(Y, 2, stats::sd) / sqrt(nrow(Y))
  stat <- mu / (se + 1e-12)
  p <- 2 * stats::pt(-abs(stat), df = nrow(Y) - 1)
  list(stat = stat, p = p)
}

#' @keywords internal
#' @noRd
.dkge_haufe_fold_covariances <- function(fit, clf, Z_by_subject) {
  S <- length(Z_by_subject)
  r <- ncol(Z_by_subject[[1]])
  Sigma_by_fold <- vector("list", clf$folds$k)

  for (k in seq_len(clf$folds$k)) {
    mdl <- clf$models_by_fold[[k]]
    if (!is.null(mdl$Sigma)) {
      Sigma <- mdl$Sigma
      std <- mdl$standardize
      if (!is.null(std)) {
        D <- diag(std$sd, nrow = r)
        Sigma <- D %*% Sigma %*% D
      }
    } else {
      train_idx <- setdiff(seq_len(S), clf$folds$assignments[[k]])
      Z_train <- do.call(rbind, Z_by_subject[train_idx])
      std <- mdl$standardize
      if (!is.null(std)) {
        Z_proc <- sweep(sweep(Z_train, 2, std$mu, FUN = "-"), 2, std$sd, FUN = "/")
        Sigma_proc <- stats::cov(Z_proc)
        D <- diag(std$sd, nrow = r)
        Sigma <- D %*% Sigma_proc %*% D
      } else {
        Sigma <- stats::cov(Z_train)
      }
    }
    Sigma <- (Sigma + t(Sigma)) / 2
    diag(Sigma) <- diag(Sigma) + 1e-6
    Sigma_by_fold[[k]] <- as.matrix(Sigma)
  }

  Sigma_by_fold
}

#' @keywords internal
#' @noRd
.dkge_default_anchor_neighborhoods <- function(renderer, k_nn) {
  graph <- renderer$graph
  if (is.null(graph) || is.null(graph$W)) {
    stop("Renderer does not contain an anchor graph; supply 'neighborhoods' explicitly.")
  }
  W <- graph$W
  Q <- nrow(W)
  neighborhoods <- vector("list", Q)
  for (j in seq_len(Q)) {
    nz <- which(W[j, ] != 0)
    if (length(nz) == 0) {
      neighborhoods[[j]] <- j
      next
    }
    if (!is.null(k_nn) && length(nz) > k_nn) {
      vals <- as.numeric(W[j, nz])
      ord <- order(-vals)
      nz <- nz[ord[seq_len(k_nn)]]
    }
    neighborhoods[[j]] <- unique(c(j, nz))
  }
  neighborhoods
}
