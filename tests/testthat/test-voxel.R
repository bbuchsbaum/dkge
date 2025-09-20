# test-voxel.R

library(testthat)

make_voxel_fixture <- function(S = 3, q = 2, P = 3, V = 4, seed = 303) {
  set.seed(seed)
  betas <- replicate(S, matrix(rnorm(q * P), q, P), simplify = FALSE)
  designs <- replicate(S, {
    X <- matrix(rnorm(40 * q), 40, q)
    qr.Q(qr(X))
  }, simplify = FALSE)
  centroids <- replicate(S, matrix(runif(P * 3), P, 3), simplify = FALSE)
  vox_feat <- replicate(S, matrix(runif(V * q), V, q), simplify = FALSE)
  fit <- dkge_fit(betas, designs, K = diag(q), rank = q)
  fit$centroids <- centroids
  list(fit = fit,
       values = lapply(betas, function(B) B[1, ]),
       voxels = vox_feat)
}

test_that("dkge_transport_to_voxels maps values", {
  fx <- make_voxel_fixture()
  res <- dkge_transport_to_voxels(fx$fit,
                                  values = fx$values,
                                  voxels = fx$voxels,
                                  mapper = "ridge")
  expect_equal(nrow(res$subj_values), length(fx$values))
  expect_equal(length(res$value), ncol(res$subj_values))
})
