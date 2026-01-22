# dkge-anchor-fit.R
# User-facing glue between anchor kernels and the DKGE core.

#' Fit DKGE using feature-anchored subject kernels
#'
#' Projects item-level kernels onto a shared anchor basis (via
#' [dkge_build_anchor_kernels()]) and reuses [dkge_fit_from_kernels()] to enter
#' the DKGE pipeline. Anchor provenance, coverage diagnostics, and subject item
#' counts are stored in the resulting `dkge` object.
#'
#' @param features_list List of subject feature matrices (`n_s x d` each).
#' @param K_item_list List of subject item kernels (`n_s x n_s` each).
#' @param folds Optional fold specification passed to
#'   [dkge_build_anchor_kernels()].
#' @param anchors Named list overriding anchor-building defaults (`L`,
#'   `method`, `rho`, `fill`, `seed`, `sigma`, `center`, `whiten`, `eps`,
#'   `unit_trace`, `item_weights`).
#' @param design_kernel Optional design kernel forwarded to
#'   [dkge_fit_from_kernels()]. Defaults to the identity.
#' @param dkge_args Named list of additional arguments forwarded to
#'   [dkge_fit_from_kernels()].
#'
#' @return A `dkge` fit with anchor provenance under `fit$provenance$anchors`.
#' @examples
#' \donttest{
#' set.seed(1)
#' features_list <- list(
#'   s1 = matrix(rnorm(30 * 5), 30, 5),
#'   s2 = matrix(rnorm(40 * 5), 40, 5),
#'   s3 = matrix(rnorm(35 * 5), 35, 5)
#' )
#' K_item_list <- lapply(features_list, function(X) {
#'   Z <- matrix(rnorm(nrow(X) * 4), nrow(X), 4)
#'   tcrossprod(Z)
#' })
#' fit <- dkge_anchor_fit(
#'   features_list, K_item_list,
#'   anchors = list(L = 8, method = "dkpp", seed = 1),
#'   dkge_args = list(rank = 2)
#' )
#' dkge_anchor_diagnostics(fit)$summary
#' }
#' @export
dkge_anchor_fit <- function(features_list,
                            K_item_list,
                            folds = NULL,
                            anchors = list(),
                            design_kernel = NULL,
                            dkge_args = list()) {
  default_anchors <- list(
    L = 128L,
    method = "dkpp",
    rho = 0.5,
    fill = "kcenter",
    seed = 1L,
    sigma = NULL,
    center = TRUE,
    whiten = TRUE,
    eps = 1e-6,
    unit_trace = TRUE,
    item_weights = NULL
  )
  anchor_args <- utils::modifyList(default_anchors, anchors)

  build_args <- c(list(features_list = features_list,
                       K_item_list = K_item_list,
                       folds = folds),
                  anchor_args)
  built <- do.call(dkge_build_anchor_kernels, build_args)
  primary <- built[[1]]

  K_aligned <- primary$K_aligned
  effect_ids <- primary$anchor_ids
  subject_ids <- names(K_aligned)

  fit_args <- c(list(K_list = K_aligned,
                     effect_ids = effect_ids,
                     subject_ids = subject_ids,
                     design_kernel = design_kernel),
                dkge_args)
  fit <- do.call(dkge_fit_from_kernels, fit_args)

  provenance <- fit$provenance %||% list()
  provenance$anchors <- list(
    method = anchor_args$method,
    args = anchor_args,
    sigma = primary$sigma,
    anchors = primary$anchors,
    anchor_ids = effect_ids,
    L = nrow(primary$anchors),
    folds = lapply(built, function(ctx) {
      list(train_idx = ctx$train_idx,
           test_idx = ctx$test_idx,
           sigma = ctx$sigma)
    }),
    item_counts = vapply(features_list, nrow, integer(1)),
    coverage = .dkge_anchor_coverage_summary(features_list, primary$anchors),
    leverage = .dkge_anchor_leverage(K_aligned)
  )
  fit$provenance <- provenance
  fit
}

#' Summarise per-subject anchor coverage
#'
#' @keywords internal
#' @noRd
.dkge_anchor_coverage_summary <- function(features_list, anchors) {
  L <- nrow(anchors)
  subject_ids <- names(features_list)
  if (is.null(subject_ids)) {
    subject_ids <- paste0("subject", seq_len(length(features_list)))
  }
  coverage <- lapply(seq_along(features_list), function(i) {
    Fs <- features_list[[i]]
    if (!is.matrix(Fs) || !nrow(Fs)) {
      return(data.frame(subject = subject_ids[[i]], p50 = NA_real_, p90 = NA_real_, p95 = NA_real_))
    }
    XX <- rowSums(Fs^2)
    ZZ <- rowSums(anchors^2)
    D2 <- outer(XX, ZZ, "+") - 2 * (Fs %*% t(anchors))
    nearest <- sqrt(pmax(apply(D2, 1L, min), 0))
    stats <- stats::quantile(nearest, probs = c(0.5, 0.9, 0.95), names = FALSE, type = 7)
    data.frame(subject = subject_ids[[i]], p50 = stats[[1]], p90 = stats[[2]], p95 = stats[[3]])
  })
  do.call(rbind, coverage)
}

#' Aggregate anchor leverage from aligned kernels
#'
#' @keywords internal
#' @noRd
.dkge_anchor_leverage <- function(K_aligned) {
  if (!length(K_aligned)) {
    return(list(mean_diag = numeric(0), per_anchor = numeric(0)))
  }
  L <- nrow(K_aligned[[1]])
  acc <- Reduce(`+`, K_aligned)
  acc <- acc / length(K_aligned)
  list(mean_diag = sum(diag(acc)) / L, per_anchor = diag(acc))
}

#' Build an anchor contrast from prototype feature sets
#'
#' @param anchors Matrix of anchor coordinates (`L x d`).
#' @param positives Matrix (or vector) of positive prototypes in feature space.
#' @param negatives Optional matrix (or vector) of negative prototypes.
#' @param sigma Optional bandwidth for the prototype kernel. Defaults to the
#'   median distance between anchors and prototypes.
#' @param normalize Logical; L2-normalise the resulting contrast.
#'
#' @return Numeric vector of length `L` suitable for [dkge_contrast()].
#' @examples
#' anchors <- matrix(rnorm(10 * 4), 10, 4)
#' pos <- anchors[1:2, , drop = FALSE]
#' neg <- anchors[3:4, , drop = FALSE]
#' w <- dkge_anchor_contrast_from_prototypes(anchors, positives = pos, negatives = neg)
#' length(w)
#' @export
dkge_anchor_contrast_from_prototypes <- function(anchors,
                                                  positives,
                                                  negatives = NULL,
                                                  sigma = NULL,
                                                  normalize = TRUE) {
  stopifnot(is.matrix(anchors))
  if (is.vector(positives)) positives <- matrix(positives, nrow = 1)
  stopifnot(is.matrix(positives), ncol(positives) == ncol(anchors))
  if (!is.null(negatives)) {
    if (is.vector(negatives)) negatives <- matrix(negatives, nrow = 1)
    stopifnot(is.matrix(negatives), ncol(negatives) == ncol(anchors))
  }
  if (is.null(sigma)) {
    sigma <- .dkge_anchor_bandwidth(rbind(anchors, positives, negatives %||% positives), anchors)
  }
  pos_score <- rowMeans(.dKGE_anchor_rbf_kernel(anchors, positives, sigma = sigma))
  neg_score <- if (is.null(negatives)) 0 else rowMeans(.dKGE_anchor_rbf_kernel(anchors, negatives, sigma = sigma))
  contrast <- pos_score - neg_score
  if (normalize) {
    norm <- sqrt(sum(contrast^2))
    if (is.finite(norm) && norm > 0) {
      contrast <- contrast / norm
    }
  }
  contrast
}

#' Build an anchor contrast from a feature-space direction
#'
#' @param anchors Matrix of anchor coordinates (`L x d`).
#' @param direction Numeric vector (length `d`) describing a linear probe in the
#'   shared feature space.
#' @param sigma Optional bandwidth; defaults to the anchor median heuristic.
#' @param normalize Logical; L2-normalise the resulting contrast.
#'
#' @return Numeric vector of length `L` suitable for [dkge_contrast()].
#' @examples
#' anchors <- matrix(rnorm(10 * 4), 10, 4)
#' direction <- rnorm(4)
#' w <- dkge_anchor_contrast_from_direction(anchors, direction)
#' length(w)
#' @export
dkge_anchor_contrast_from_direction <- function(anchors,
                                                 direction,
                                                 sigma = NULL,
                                                 normalize = TRUE) {
  stopifnot(is.matrix(anchors), is.numeric(direction), length(direction) == ncol(anchors))
  if (is.null(sigma)) {
    sigma <- .dkge_anchor_bandwidth(rbind(anchors, matrix(direction, nrow = 1)), anchors)
  }
  scores <- as.numeric(.dKGE_anchor_rbf_kernel(anchors, matrix(direction, nrow = 1), sigma = sigma))
  if (normalize) {
    norm <- sqrt(sum(scores^2))
    if (is.finite(norm) && norm > 0) scores <- scores / norm
  }
  scores
}

#' Extract anchor diagnostics from a DKGE fit
#'
#' @param fit Object returned by [dkge_anchor_fit()] or [dkge_fit_from_kernels()].
#'
#' @return List containing per-subject coverage quantiles and per-anchor leverage
#'   estimates.
#' @examples
#' \donttest{
#' set.seed(1)
#' features_list <- list(
#'   s1 = matrix(rnorm(20 * 5), 20, 5),
#'   s2 = matrix(rnorm(25 * 5), 25, 5),
#'   s3 = matrix(rnorm(22 * 5), 22, 5)
#' )
#' K_item_list <- lapply(features_list, function(X) tcrossprod(matrix(rnorm(nrow(X) * 3), nrow(X), 3)))
#' fit <- dkge_anchor_fit(features_list, K_item_list,
#'                        anchors = list(L = 6, method = "dkpp", seed = 1),
#'                        dkge_args = list(rank = 2))
#' diag <- dkge_anchor_diagnostics(fit)
#' names(diag)
#' }
#' @export
dkge_anchor_diagnostics <- function(fit) {
  stopifnot(inherits(fit, "dkge"))
  info <- fit$provenance$anchors
  if (is.null(info)) {
    stop("Fit does not contain anchor provenance.", call. = FALSE)
  }
  leverage <- info$leverage$per_anchor %||% numeric(0)
  list(
    summary = list(
      method = info$method,
      sigma = info$sigma,
      L = info$L,
      mean_item_count = mean(info$item_counts)
    ),
    coverage = info$coverage,
    leverage = data.frame(anchor = info$anchor_ids,
                          leverage = leverage)
  )
}
