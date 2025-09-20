# test-backproj.R
# Diagnostics for dkge-backproj helpers

library(testthat)
library(dkge)

make_backproj_fixture <- function(S = 3, q = 3, P = 4, T = 24, seed = 4242) {
  set.seed(seed)
  effects <- paste0("eff", seq_len(q))
  betas <- replicate(S, {
    mat <- matrix(rnorm(q * P), q, P)
    rownames(mat) <- effects
    mat
  }, simplify = FALSE)
  designs <- replicate(S, {
    X <- matrix(rnorm(T * q), T, q)
    X <- qr.Q(qr(X))
    colnames(X) <- effects
    X
  }, simplify = FALSE)
  centroids <- replicate(S, matrix(runif(P * 3), P, 3), simplify = FALSE)

  fit <- dkge_fit(dkge_data(betas, designs = designs),
                  K = diag(q),
                  rank = 2,
                  keep_X = TRUE)
  loadings <- dkge_predict_loadings(fit, betas)

  list(fit = fit,
       betas = betas,
       centroids = centroids,
       loadings = loadings,
       S = S,
       P = P)
}

test_that("dkge_transport_loadings_to_medoid keeps medoid subject unchanged", {
  fixture <- make_backproj_fixture()
  res <- dkge_transport_loadings_to_medoid(
    fixture$fit,
    medoid = 1,
    centroids = fixture$centroids,
    loadings = fixture$loadings,
    method = "sinkhorn"
  )

  expect_length(res$subjects, fixture$fit$rank)
  expect_equal(dim(res$subjects[[1]]), c(fixture$S, fixture$P))
  expect_equal(res$subjects[[1]][1, ], fixture$loadings[[1]][, 1],
               tolerance = 1e-6)
  expect_equal(res$group[[1]], apply(res$subjects[[1]], 2, stats::median),
               tolerance = 1e-12)
})

test_that("dkge_transport_contrasts_to_medoid returns aligned subject maps", {
  fixture <- make_backproj_fixture(seed = 2024)
  contrast_obj <- dkge_contrast(fixture$fit, c(1, -1, 0), method = "loso")

  res <- dkge_transport_contrasts_to_medoid(
    fixture$fit,
    contrast_obj,
    medoid = 1,
    centroids = fixture$centroids,
    betas = fixture$betas,
    method = "sinkhorn"
  )

  expect_named(res, names(contrast_obj$values))
  first <- res[[1]]
  expect_equal(dim(first$subj_values), c(fixture$S, fixture$P))
  expect_equal(first$subj_values[1, ], contrast_obj$values[[1]][[1]],
               tolerance = 1e-6)
  expect_equal(first$plans[[1]], diag(1, fixture$P), tolerance = 1e-6)
})

test_that("dkge_transport_loadings_to_medoid falls back to stored loadings", {
  fixture <- make_backproj_fixture()
  fixture$fit$input$betas <- NULL
  res <- dkge_transport_loadings_to_medoid(
    fixture$fit,
    medoid = 1,
    centroids = fixture$centroids,
    method = "sinkhorn"
  )

  expect_equal(res$subjects[[1]][1, ], fixture$loadings[[1]][, 1], tolerance = 1e-6)
})
