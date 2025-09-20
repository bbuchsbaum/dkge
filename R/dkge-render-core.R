# dkge-render-core.R
# Dense anchor rendering pipeline: anchors, graph, decoder, aggregation, renderer helpers.

#' Build or validate anchor coordinates in MNI space
#'
#' @param xyz Optional \eqn{N \times 3} matrix of grey-matter coordinates (in mm).
#'   Used when anchors need to be derived from voxel locations.
#' @param anchors Optional precomputed \eqn{Q \times 3} anchor matrix. When
#'   supplied it is validated and returned unchanged.
#' @param n_anchor Target number of anchors when deriving them from `xyz`.
#' @param method Selection strategy when deriving anchors. `"kmeans"`
#'   (default) uses k-means centroids, `"sample"` performs a uniform
#'   subsample.
#' @param seed Optional random seed for reproducibility.
#' @return A numeric \eqn{Q \times 3} matrix of anchor coordinates.
#' @export
dkge_make_anchors <- function(xyz = NULL, anchors = NULL,
                              n_anchor = 20000L,
                              method = c("kmeans", "sample"),
                              seed = NULL) {
  if (!is.null(anchors)) {
    stopifnot(is.matrix(anchors), ncol(anchors) == 3)
    storage.mode(anchors) <- "double"
    return(anchors)
  }

  stopifnot(is.matrix(xyz), ncol(xyz) == 3)
  method <- match.arg(method)
  if (!is.null(seed)) set.seed(seed)

  N <- nrow(xyz)
  if (method == "sample") {
    idx <- if (N <= n_anchor) seq_len(N) else sample.int(N, n_anchor)
    A <- xyz[idx, , drop = FALSE]
  } else {
    centers <- min(n_anchor, N)
    km <- stats::kmeans(xyz, centers = centers, iter.max = 25, nstart = 1)
    A <- km$centers
  }
  storage.mode(A) <- "double"
  A
}

#' Construct a kNN anchor graph and Laplacian
#'
#' @param anchors \eqn{Q \times 3} matrix of anchor coordinates.
#' @param k Number of neighbours used in the graph (default 10).
#' @param sigma Optional Gaussian length-scale (mm). If `NULL`, the
#'   `neighborweights` defaults are used.
#' @param weight_mode Edge weighting scheme passed to
#'   [neighborweights::graph_weights()]. Defaults to `"heat"`.
#' @param type Graph symmetrisation strategy (see
#'   [neighborweights::graph_weights()]); `"mutual"` helps enforce
#'   symmetry.
#' @return A list containing the neighbour graph, sparse adjacency `W`,
#'   degree matrix `D`, and Laplacian `L`.
#' @export
dkge_anchor_graph <- function(anchors,
                              k = 10,
                              sigma = NULL,
                              weight_mode = c("heat", "normalized", "binary"),
                              type = c("mutual", "normal", "asym")) {
  stopifnot(is.matrix(anchors), ncol(anchors) == 3, k >= 1)
  weight_mode <- match.arg(weight_mode)
  type <- match.arg(type)
  suppressMessages({
    graph <- neighborweights::graph_weights(
      anchors,
      k = k,
      sigma = sigma,
      neighbor_mode = "knn",
      weight_mode = weight_mode,
      type = type
    )
  })
  W <- neighborweights::adjacency(graph)
  # Ensure symmetry in case of floating point mismatches
  W <- (W + Matrix::t(W)) / 2
  D <- Matrix::Diagonal(nrow(W), x = Matrix::rowSums(W))
  L <- neighborweights::laplacian(graph)
  list(graph = graph, W = W, D = D, L = L)
}

#' Fit a sparse anchor-to-voxel decoder
#'
#' Each voxel is represented as a convex combination of its `k` nearest
#' anchors using Gaussian weights. The result can be reused to decode
#' any anchor-level map.
#'
#' @param anchors Anchor coordinate matrix (`Q x 3`).
#' @param vox_xyz Voxel coordinate matrix (`V x 3`).
#' @param k Number of anchors per voxel (default 8).
#' @param sigma Optional Gaussian length-scale. If `NULL`, it is set to the
#'   square root of the median squared distance between voxels and their
#'   nearest anchors.
#' @return A decoder object storing neighbour indices, weights, and a sparse
#'   matrix implementing the transformation.
#' @export
dkge_anchor_to_voxel_fit <- function(anchors, vox_xyz, k = 8, sigma = NULL) {
  stopifnot(is.matrix(anchors), ncol(anchors) == 3,
            is.matrix(vox_xyz), ncol(vox_xyz) == 3,
            k >= 1)
  knn <- FNN::get.knnx(anchors, vox_xyz, k = k)
  dist2 <- knn$nn.dist^2
  if (is.null(sigma)) {
    sigma <- sqrt(stats::median(dist2))
  }
  weights <- exp(-dist2 / (2 * sigma^2))
  weights <- weights / (rowSums(weights) + 1e-12)

  V <- nrow(vox_xyz)
  idx <- knn$nn.index
  i <- rep(seq_len(V), each = k)
  j <- as.vector(t(idx))
  x <- as.vector(t(weights))
  H <- Matrix::sparseMatrix(i = i, j = j, x = x,
                            dims = c(V, nrow(anchors)))

  list(
    idx = idx,
    weights = weights,
    sparse = H,
    n_anchors = nrow(anchors),
    n_vox = V,
    params = list(k = k, sigma = sigma)
  )
}

#' Decode anchor values to voxel space
#'
#' @param decoder Object returned by [dkge_anchor_to_voxel_fit()].
#' @param anchor_values Numeric vector of length `decoder$n_anchors`.
#' @return Numeric vector of voxel values.
#' @export
dkge_anchor_to_voxel_apply <- function(decoder, anchor_values) {
  stopifnot(is.numeric(anchor_values),
            length(anchor_values) == decoder$n_anchors)
  if (!is.null(decoder$sparse)) {
    return(as.numeric(decoder$sparse %*% anchor_values))
  }
  # Fallback for decoders constructed without a sparse matrix
  V <- decoder$n_vox
  out <- numeric(V)
  for (v in seq_len(V)) {
    out[v] <- sum(decoder$weights[v, ] * anchor_values[decoder$idx[v, ]])
  }
  out
}

#' Aggregate anchor fields with optional Laplacian smoothing
#'
#' @param anchor_list List of anchor-valued vectors (length `Q`).
#' @param subj_weights Optional subject weights applied during the average.
#' @param L Optional graph Laplacian for Tikhonov regularisation.
#' @param lambda Non-negative smoothing parameter. `0` disables smoothing.
#' @return A list with the smoothed field `y`, the raw weighted mean `ybar`,
#'   `coverage` (weighted contribution counts), and placeholder `ess` values.
#' @export
dkge_anchor_aggregate <- function(anchor_list,
                                  subj_weights = NULL,
                                  L = NULL,
                                  lambda = 0) {
  stopifnot(length(anchor_list) >= 1L)
  Q <- length(anchor_list[[1]])
  stopifnot(all(vapply(anchor_list, length, integer(1)) == Q))

  S <- length(anchor_list)
  w <- subj_weights %||% rep(1, S)
  stopifnot(length(w) == S)

  weighted <- Map(function(y, ws) ws * y, anchor_list, as.list(w))
  num <- Reduce(`+`, weighted)
  den <- sum(w) + 1e-12
  ybar <- num / den

  y <- ybar
  if (!is.null(L) && lambda > 0) {
    I <- Matrix::Diagonal(Q)
    A <- I + lambda * L
    y <- as.numeric(Matrix::solve(A, ybar))
  }

  coverage <- Reduce(`+`, Map(function(y, ws) abs(ws) * as.numeric(y != 0),
                              anchor_list, as.list(w)))
  ess <- rep(NA_real_, Q)
  list(y = y, ybar = ybar, coverage = coverage, ess = ess)
}

#' Prepare reusable rendering objects for a fitted DKGE model
#'
#' @param fit Fitted `dkge` object.
#' @param centroids List of per-subject centroid matrices (`P_s x 3`).
#'   Must align with `fit$Btil`.
#' @param anchors Optional precomputed anchor coordinate matrix. When `NULL`,
#'   anchors are derived from `anchor_xyz` if provided, otherwise from
#'   `vox_xyz`.
#' @param anchor_xyz Optional matrix of candidate points used to derive anchors
#'   via [dkge_make_anchors()]. Ignored when `anchors` is supplied.
#' @param anchor_n Number of anchors to draw when constructing them from
#'   coordinates.
#' @param anchor_method Method passed to [dkge_make_anchors()] when anchors are
#'   derived. Defaults to `"kmeans"`.
#' @param anchor_seed Optional seed forwarded to [dkge_make_anchors()].
#' @param vox_xyz Optional voxel coordinates for constructing a decoder.
#' @param mapper Mapper specification created with [dkge_mapper()]. Defaults to
#'   the barycentric kNN mapper.
#' @param subject_feats Optional list of matrices supplying latent features per
#'   subject cluster. When provided and the mapper consumes latent information
#'   (e.g., Sinkhorn), they are forwarded via `subj_feats`.
#' @param anchor_feats Optional anchor-level feature matrix aligned with
#'   `anchors`. Derived automatically by pooling subject features when `NULL`
#'   and `subject_feats` are provided.
#' @param feat_lambda Feature cost weight passed to Sinkhorn mappers. Ignored by
#'   kNN.
#' @param feat_sigma Feature bandwidth used when computing feature costs.
#' @param reliabilities Optional list of per-subject reliability vectors passed
#'   to the mapper during fitting.
#' @param graph_k Optional integer; when provided, an anchor graph of this
#'   neighbourhood size is constructed for subsequent smoothing.
#' @param decoder_k Number of anchors per voxel when building the decoder.
#' @return A list bundling anchors, optional graph/decoder, fitted per-subject
#'   mappers, and subject weights.
#' @export
dkge_build_renderer <- function(fit,
                                centroids,
                                anchors = NULL,
                                anchor_xyz = NULL,
                                anchor_n = 20000L,
                                anchor_method = c("kmeans", "sample"),
                                anchor_seed = NULL,
                                vox_xyz = NULL,
                                mapper = dkge_mapper("knn", k = 8, sigx = 3.0),
                                graph_k = NULL,
                                decoder_k = 8,
                                reliabilities = NULL,
                                subject_feats = NULL,
                                anchor_feats = NULL,
                                feat_lambda = NULL,
                                feat_sigma = NULL) {
  stopifnot(inherits(fit, "dkge"))
  S <- length(fit$Btil)
  stopifnot(length(centroids) == S)

  anchor_method <- match.arg(anchor_method)

  anchor_mat <- if (!is.null(anchors)) {
    dkge_make_anchors(anchors = anchors)
  } else {
    xyz <- anchor_xyz %||% vox_xyz
    if (is.null(xyz)) {
      stop("Provide either 'anchors', 'anchor_xyz', or 'vox_xyz' to construct anchors.")
    }
    dkge_make_anchors(xyz = xyz,
                      n_anchor = anchor_n,
                      method = anchor_method,
                      seed = anchor_seed)
  }

  graph <- if (!is.null(graph_k)) {
    dkge_anchor_graph(anchor_mat, k = graph_k)
  } else {
    NULL
  }

  decoder <- if (!is.null(vox_xyz)) {
    dkge_anchor_to_voxel_fit(anchor_mat, vox_xyz, k = decoder_k)
  } else {
    NULL
  }

  if (!is.null(subject_feats)) {
    stopifnot(length(subject_feats) == S)
    feat_cols <- ncol(subject_feats[[1]])
    stopifnot(all(vapply(subject_feats, ncol, integer(1)) == feat_cols))
  }

  if (is.null(anchor_feats) && !is.null(subject_feats)) {
    # Pool subject features barycentrically to initialise anchor feature vectors
    feat_dim <- ncol(subject_feats[[1]])
    anchor_feats <- matrix(0, nrow = nrow(anchor_mat), ncol = feat_dim)
    weight_tot <- numeric(nrow(anchor_mat))
    mapper_for_feats <- dkge_mapper("knn", k = 8, sigx = 5)
    for (s in seq_len(S)) {
      feats <- subject_feats[[s]]
      stopifnot(nrow(feats) == nrow(centroids[[s]]))
      fit_feat <- fit_mapper(mapper_for_feats,
                             subj_points = centroids[[s]],
                             anchor_points = anchor_mat)
      contrib <- apply_mapper(fit_feat, feats, normalize_by_reliab = FALSE)
      anchor_feats <- anchor_feats + contrib
      weight_tot <- weight_tot + apply_mapper(fit_feat,
                                              rep(1, nrow(centroids[[s]])),
                                              normalize_by_reliab = FALSE)
    }
    scale_idx <- weight_tot > 1e-8
    anchor_feats[scale_idx, ] <- anchor_feats[scale_idx, , drop = FALSE] /
      weight_tot[scale_idx]
  }

  mapper_fits <- vector("list", S)
  mapper_stats <- vector("list", S)
  for (s in seq_len(S)) {
    subj_feat_s <- if (!is.null(subject_feats)) subject_feats[[s]] else NULL
    rel_s <- if (!is.null(reliabilities)) reliabilities[[s]] else NULL
    mapper_obj <- mapper
    if (inherits(mapper_obj, "dkge_mapper_sinkhorn")) {
      mapper_obj <- mapper
      mapper_obj$pars <- mapper_obj$pars
      mapper_obj$pars$lambda_feat <- feat_lambda %||% mapper_obj$pars$lambda_feat
      mapper_obj$pars$sigz <- feat_sigma %||% mapper_obj$pars$sigz
    }
    fit_s <- fit_mapper(mapper_obj,
                        subj_points = centroids[[s]],
                        anchor_points = anchor_mat,
                        subj_feats = subj_feat_s,
                        anchor_feats = anchor_feats,
                        reliab = rel_s)
    mapper_fits[[s]] <- fit_s

    stats <- NULL
    if (inherits(fit_s, "dkge_mapper_fit_sinkhorn")) {
      plan <- fit_s$plan
      plan_vals <- as.numeric(plan)
      plan_vals <- plan_vals[plan_vals > 0]
      entropy <- if (length(plan_vals)) -sum(plan_vals * log(plan_vals)) else 0
      support <- sum(plan_vals > 1e-6)
      cost <- if (!is.null(fit_s$stats)) fit_s$stats$transport_cost else NA_real_
      eps <- if (!is.null(fit_s$epsilon)) fit_s$epsilon else NA_real_
      stats <- list(transport_cost = cost,
                    plan_entropy = entropy,
                    effective_support = support,
                    epsilon = eps)
    }
    mapper_stats[[s]] <- stats
  }

  list(
    anchors = anchor_mat,
    graph = graph,
    decoder = decoder,
    mapper = mapper,
    mapper_fits = mapper_fits,
    weights = fit$weights %||% rep(1, S),
    anchor_feats = anchor_feats,
    mapper_stats = mapper_stats
  )
}

#' Render per-subject values to anchors and voxels
#'
#' @param renderer Object produced by [dkge_build_renderer()].
#' @param values_list List of per-subject value vectors (aligned with centroids).
#' @param lambda Optional Laplacian smoothing strength applied in anchor space.
#' @param to_vox Logical; when `TRUE` (default) and a decoder is available,
#'   voxel maps are produced.
#' @return A list with `anchor` (dense anchor field), optional `voxel` map,
#'   the aggregation diagnostics, and the intermediate per-subject anchor maps.
#' @export
dkge_render_subject_values <- function(renderer,
                                        values_list,
                                        lambda = 0,
                                        to_vox = TRUE) {
  stopifnot(is.list(renderer$mapper_fits))
  S <- length(renderer$mapper_fits)
  stopifnot(length(values_list) == S)

  anchor_maps <- lapply(seq_len(S), function(s) {
    apply_mapper(renderer$mapper_fits[[s]], values_list[[s]])
  })

  L <- if (!is.null(renderer$graph)) renderer$graph$L else NULL
  agg <- dkge_anchor_aggregate(anchor_maps,
                               subj_weights = renderer$weights,
                               L = L,
                               lambda = lambda)

  stats <- renderer$mapper_stats %||% vector("list", S)
  agg$subject_stats <- stats
  if (any(!vapply(stats, is.null, logical(1)))) {
    costs <- vapply(stats, function(x) x$transport_cost %||% NA_real_, numeric(1))
    entropies <- vapply(stats, function(x) x$plan_entropy %||% NA_real_, numeric(1))
    cost_valid <- costs[is.finite(costs)]
    entropy_valid <- entropies[is.finite(entropies)]
    if (length(cost_valid)) {
      agg$transport_cost_mean <- mean(cost_valid)
    }
    if (length(entropy_valid)) {
      agg$plan_entropy_mean <- mean(entropy_valid)
    }
  }

  vox <- NULL
  if (to_vox && !is.null(renderer$decoder)) {
    vox <- dkge_anchor_to_voxel_apply(renderer$decoder, agg$y)
  }

  list(anchor = agg$y,
       voxel = vox,
       details = agg,
       subject_anchor_maps = anchor_maps)
}
