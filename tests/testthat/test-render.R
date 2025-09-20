# test-render.R
# Basic validation for the dense rendering core.

library(testthat)

set.seed(1)


test_that("kNN mapper preserves constant fields on covered anchors", {
  anchors <- matrix(runif(30, -10, 10), ncol = 3)
  subj_points <- anchors[1:6, , drop = FALSE]
  mapper <- dkge_mapper("knn", k = 3, sigx = 5)
  fit <- fit_mapper(mapper,
                    subj_points = subj_points,
                    anchor_points = anchors)
  const_vals <- rep(2.5, nrow(subj_points))
  mapped <- apply_mapper(fit, const_vals)
  covered <- sort(unique(as.vector(fit$idx)))
  expect_equal(mapped[covered], rep(2.5, length(covered)), tolerance = 1e-6)
})


test_that("kNN mapper reduces to identity when anchors equal centroids", {
  anchors <- matrix(seq_len(9), ncol = 3)
  mapper <- dkge_mapper("knn", k = 1, sigx = 1)
  fit <- fit_mapper(mapper,
                    subj_points = anchors,
                    anchor_points = anchors)
  vals <- rnorm(nrow(anchors))
  mapped <- apply_mapper(fit, vals)
  expect_equal(mapped, vals, tolerance = 1e-6)
})


test_that("anchor smoothing reduces Laplacian energy", {
  anchors <- cbind(seq(0, 1, length.out = 5), 0, 0)
  graph <- dkge_anchor_graph(anchors, k = 2)
  y1 <- c(-1, -0.5, 0, 0.5, 1)
  y2 <- rev(y1)
  agg_raw <- dkge_anchor_aggregate(list(y1, y2), lambda = 0)
  agg_smooth <- dkge_anchor_aggregate(list(y1, y2), L = graph$L, lambda = 0.5)
  energy_raw <- as.numeric(t(agg_raw$y) %*% (graph$L %*% agg_raw$y))
  energy_smooth <- as.numeric(t(agg_smooth$y) %*% (graph$L %*% agg_smooth$y))
  expect_lt(energy_smooth, energy_raw + 1e-8)
})

test_that("renderer builds anchors from voxels when not supplied", {
  fit_stub <- structure(list(Btil = list(matrix(0, 3, 1)),
                             weights = 1),
                        class = "dkge")
  centroids <- list(matrix(runif(9, -10, 10), ncol = 3))
  vox_xyz <- matrix(runif(300, -40, 40), ncol = 3)

  renderer <- dkge_build_renderer(fit_stub,
                                 centroids = centroids,
                                 anchors = NULL,
                                 vox_xyz = vox_xyz,
                                 anchor_n = 15L,
                                 anchor_method = "sample",
                                 anchor_seed = 1L,
                                 mapper = dkge_mapper("knn", k = 3, sigx = 5))

  expect_true(is.matrix(renderer$anchors))
  expect_equal(nrow(renderer$anchors), min(15L, nrow(vox_xyz)))
  expect_false(is.null(renderer$decoder))
})

test_that("barycentric aggregation recovers smooth analytic field", {
  anchors <- cbind(seq(0, 4), 0, 0)
  f <- function(x) 2 * x - 1
  mapper <- dkge_mapper("knn", k = 3, sigx = 0.3)

  jitter1 <- c(-0.05, 0.02, 0.04, -0.03, 0.01)
  jitter2 <- c(0.06, -0.04, 0.03, 0.05, -0.06)
  subj_points1 <- anchors
  subj_points1[, 1] <- subj_points1[, 1] + jitter1
  subj_points2 <- anchors
  subj_points2[, 1] <- subj_points2[, 1] + jitter2

  vals <- f(anchors[, 1])
  fit1 <- fit_mapper(mapper, subj_points = subj_points1, anchor_points = anchors)
  fit2 <- fit_mapper(mapper, subj_points = subj_points2, anchor_points = anchors)
  y1 <- apply_mapper(fit1, vals)
  y2 <- apply_mapper(fit2, vals)

  agg <- dkge_anchor_aggregate(list(y1, y2))
  expect_equal(agg$y, vals, tolerance = 1e-2)
})


test_that("reliability weights dampen uncertain contributions", {
  anchors <- rbind(c(0, 0, 0), c(2, 0, 0))
  subj_points <- rbind(c(0.1, 0, 0), c(1.9, 0, 0))
  mapper <- dkge_mapper("knn", k = 2, sigx = 5)
  fit <- fit_mapper(mapper, subj_points = subj_points, anchor_points = anchors)

  values <- c(4, -2)
  reliab <- c(1, 0.1)
  mapped_unweighted <- apply_mapper(fit, values, normalize_by_reliab = FALSE)
  mapped_weighted <- apply_mapper(fit, values, reliab = reliab, normalize_by_reliab = TRUE)

  W <- fit$weights
  idx <- fit$idx
  numer <- numeric(fit$Q)
  denom <- numeric(fit$Q)
  rv <- reliab * values
  for (t in seq_len(ncol(idx))) {
    js <- idx[, t]
    wt <- W[, t]
    numer[js] <- numer[js] + wt * rv
    denom[js] <- denom[js] + wt * reliab
  }
  expected_weighted <- numeric(fit$Q)
  keep <- denom > 1e-12
  expected_weighted[keep] <- numer[keep] / denom[keep]

  expect_false(isTRUE(all.equal(mapped_weighted, mapped_unweighted)))
  expect_equal(mapped_weighted, expected_weighted, tolerance = 1e-8)
})


test_that("anchor-to-voxel decoder preserves aligned values", {
  anchors <- rbind(c(0, 0, 0), c(1, 0, 0), c(2, 0, 0))
  vox_xyz <- anchors
  decoder <- dkge_anchor_to_voxel_fit(anchors, vox_xyz, k = 1, sigma = 1)
  anchor_values <- c(5, -1.5, 3.2)
  vox_values <- dkge_anchor_to_voxel_apply(decoder, anchor_values)
  expect_equal(vox_values, anchor_values, tolerance = 1e-10)
})

test_that("sinkhorn renderer uses latent features and reports diagnostics", {
  fit_stub <- structure(list(Btil = list(matrix(0, 2, 1), matrix(0, 2, 1)),
                             weights = c(1, 1)),
                        class = "dkge")
  centroids <- list(matrix(0, nrow = 2, ncol = 3), matrix(0, nrow = 2, ncol = 3))
  anchors <- matrix(0, nrow = 2, ncol = 3)
  subj_feats <- list(rbind(c(1, 0), c(0, 1)),
                     rbind(c(1, 0), c(0, 1)))
  anchor_feats <- rbind(c(1, 0), c(0, 1))

  mapper <- dkge_mapper("sinkhorn",
                        epsilon = 1e-3,
                        lambda_xyz = 0,
                        lambda_feat = 1,
                        sigz = 1)

  renderer <- dkge_build_renderer(fit_stub,
                                 centroids = centroids,
                                 anchors = anchors,
                                 mapper = mapper,
                                 subject_feats = subj_feats,
                                 anchor_feats = anchor_feats,
                                 feat_lambda = 1,
                                 feat_sigma = 1)

  expect_equal(renderer$anchor_feats, anchor_feats, tolerance = 1e-8)
  expect_false(all(vapply(renderer$mapper_stats, is.null, logical(1))))

  values_list <- list(c(1, -1), c(-1, 1))
  rendered <- dkge_render_subject_values(renderer, values_list, lambda = 0, to_vox = FALSE)
  stats <- rendered$details$subject_stats
  expect_true(all(vapply(stats, function(x) !is.null(x$plan_entropy), logical(1))))
  expect_true(rendered$details$plan_entropy_mean > 0)
})
