# test-project.R
# Diagnostics for projection and prediction helpers

library(testthat)
library(dkge)

make_projection_fixture <- function(S = 3, q = 3, P = 4, Tn = 40, seed = 123) {
  set.seed(seed)
  effects <- paste0("eff", seq_len(q))
  betas <- replicate(S, {
    mat <- matrix(rnorm(q * P), q, P)
    rownames(mat) <- effects
    mat
  }, simplify = FALSE)
  designs <- replicate(S, {
    X <- matrix(rnorm(Tn * q), Tn, q)
    X <- qr.Q(qr(X))
    colnames(X) <- effects
    X
  }, simplify = FALSE)
  omega <- replicate(S, runif(P, min = 0.5, max = 1.5), simplify = FALSE)
  centroids <- replicate(S, matrix(runif(P * 3), P, 3), simplify = FALSE)
  list(betas = betas, designs = designs, omega = omega, effects = effects,
       centroids = centroids)
}

context_fixture <- make_projection_fixture()

fit_fixture <- dkge_fit(dkge_data(context_fixture$betas,
                                  designs = context_fixture$designs,
                                  omega = context_fixture$omega),
                        K = diag(length(context_fixture$effects)),
                        rank = 2,
                        keep_X = TRUE)

test_that("subject weights shrink correctly", {
  Btil <- list(matrix(1, 2, 3), matrix(2, 2, 3))
  Khalf <- diag(2)
  w_raw_expected <- c(1 / sum(Btil[[1]]^2), 1 / sum(Btil[[2]]^2))
  w_norm_expected <- w_raw_expected / mean(w_raw_expected)
  w0 <- dkge:::.dkge_subject_weights(Btil, list(NULL, NULL), Khalf,
                                     w_method = "energy", w_tau = 0)
  w1 <- dkge:::.dkge_subject_weights(Btil, list(NULL, NULL), Khalf,
                                     w_method = "energy", w_tau = 1)
  expect_equal(w0, w_norm_expected, tolerance = 1e-10)
  expect_equal(w1, rep(1, 2))
})

test_that("Omega vector accumulation matches sqrt weighting", {
  Btil <- list(matrix(c(1, 2, 3, 4), 2, 2))
  Omega <- list(c(1.5, 0.5))
  Khalf <- diag(2)
  res <- dkge:::.dkge_accumulate_chat(Btil, Omega, Khalf, weights = 1)
  W <- Btil[[1]] * rep(sqrt(Omega[[1]]), each = 2)
  expected <- Khalf %*% (W %*% t(W)) %*% Khalf
  expect_equal(res$contribs[[1]], expected)
  expect_equal(res$Chat, expected)
})

test_that("preprocess recreates training blocks", {
  X_new <- dkge_preprocess_blocks(fit_fixture,
                                  context_fixture$betas,
                                  fit_fixture$Omega,
                                  fit_fixture$weights)
  expect_equal(X_new, fit_fixture$X_concat, tolerance = 1e-10, check.attributes = FALSE)
})

test_that("project blocks reproduces training scores", {
  scores <- dkge_project_blocks(fit_fixture,
                                context_fixture$betas,
                                fit_fixture$Omega,
                                fit_fixture$weights)
  expect_equal(scores, fit_fixture$s, tolerance = 1e-10, check.attributes = FALSE)
})

test_that("cluster projection matches training loadings", {
  for (s in seq_along(context_fixture$betas)) {
    osm <- context_fixture$omega[[s]]
    for (j in seq_len(ncol(context_fixture$betas[[s]]))) {
      omega_j <- if (is.null(osm)) 1 else osm[j]
      proj <- dkge_project_cluster(fit_fixture,
                                   context_fixture$betas[[s]][, j],
                                   omega = omega_j,
                                   w = fit_fixture$weights[s])
  expect_equal(proj, fit_fixture$v[fit_fixture$block_indices[[s]][j], ],
               tolerance = 1e-8, check.attributes = FALSE)
    }
  }
})

test_that("predict loadings reproduce training loadings", {
  preds <- dkge_predict_loadings(fit_fixture, context_fixture$betas)
  expected <- lapply(fit_fixture$Btil, function(Bts) t(Bts) %*% fit_fixture$K %*% fit_fixture$U)
  expect_equal(preds, expected, tolerance = 1e-10, check.attributes = FALSE)
})

test_that("predict respects effect reordering", {
  shuffled <- lapply(context_fixture$betas, function(mat) mat[rev(seq_len(nrow(mat))), , drop = FALSE])
  preds_shuffled <- dkge_predict_loadings(fit_fixture, shuffled)
  preds <- dkge_predict_loadings(fit_fixture, context_fixture$betas)
  expect_equal(preds_shuffled, preds, tolerance = 1e-10, check.attributes = FALSE)
})

test_that("dkge_predict preserves subject identifiers", {
  subjects <- lapply(seq_along(context_fixture$betas), function(s) {
    dkge_subject(context_fixture$betas[[s]],
                 design = context_fixture$designs[[s]],
                 id = paste0("sub", s))
  })
  contrast <- diag(length(context_fixture$effects))[, 1]
  res <- dkge_predict(fit_fixture, subjects, list(c1 = contrast))
  expect_equal(names(res$values), paste0("sub", seq_along(subjects)))
  expect_equal(names(res$A_list), paste0("sub", seq_along(subjects)))
})


test_that("latent loadings transport to medoid", {
  centroids <- context_fixture$centroids
  loadings <- dkge_predict_loadings(fit_fixture, context_fixture$betas)
  transport <- dkge_transport_loadings_to_medoid(fit_fixture, medoid = 1,
                                                 centroids = centroids,
                                                 loadings = loadings)
  first_component <- transport$subjects[[1]][1, ]
  expect_equal(first_component, loadings[[1]][, 1], tolerance = 1e-6)
})

test_that("contrast transport to medoid matches loadings", {
  centroids <- context_fixture$centroids
  contrast <- dkge_contrast(fit_fixture, c(1, -1, 0), method = "loso")
  res <- dkge_transport_contrasts_to_medoid(fit_fixture, contrast, medoid = 1,
                                            centroids = centroids,
                                            betas = context_fixture$betas)
  first_subject <- res[[1]]$subj_values[1, ]
  base_vals <- contrast$values[[1]][[1]]
  expect_equal(first_subject, base_vals, tolerance = 1e-6)
})
