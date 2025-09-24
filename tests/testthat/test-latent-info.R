# test-latent-info.R

library(testthat)

make_demo_fit <- function() {
  q <- 3
  r <- 2
  B_list <- list(
    matrix(c(1.0, 0.0, 0.0, 0.2, 1.0, 0.1), nrow = q, ncol = 2),
    matrix(c(0.9, 0.1, 0.0, -0.1, 0.8, 0.2), nrow = q, ncol = 2),
    matrix(c(-0.8, 0.1, 0.2, 0.3, -0.6, 0.0), nrow = q, ncol = 2),
    matrix(c(-0.9, -0.2, 0.1, 0.2, -0.5, 0.1), nrow = q, ncol = 2)
  )
  structure(list(
    Btil = B_list,
    K = diag(q),
    U = cbind(c(1, 0, 0), c(0, 1, 0)),
    weights = rep(1, length(B_list)),
    subject_ids = paste0("s", seq_along(B_list))
  ), class = "dkge")
}

make_demo_renderer <- function(fit) {
  anchors <- rbind(
    c(0, 0, 0),
    c(1, 0, 0),
    c(0, 1, 0),
    c(0, 0, 1)
  )
  centroids <- list(
    rbind(c(0, 0, 0), c(1, 0, 0)),
    rbind(c(0, 1, 0), c(1, 1, 0)),
    rbind(c(0, 0, 1), c(1, 0, 1)),
    rbind(c(0, 1, 1), c(1, 1, 1))
  )
  dkge_build_renderer(fit,
                      centroids = centroids,
                      anchors = anchors,
                      mapper = dkge_mapper("knn", k = 2, sigx = 2.0),
                      graph_k = 2)
}

make_demo_classifier <- function(betas) {
  S <- length(betas)
  folds <- structure(list(k = 1, assignments = list(seq_len(S))),
                     class = "dkge_folds")
  mdl <- list(
    Sigma = diag(length(betas[[1]])),
    standardize = NULL,
    model = "lda",
    classes = c("a", "b")
  )
  structure(list(
    model = "lda",
    level = "subject",
    folds = folds,
    models_by_fold = list(mdl),
    beta_by_subject = betas,
    fold_assignment = rep(1L, S),
    r = length(betas[[1]])
  ), class = "dkge_clf")
}

set.seed(123)


test_that("latent projections have expected shapes", {
  fit <- make_demo_fit()
  Z_list <- dkge_project_clusters_to_latent(fit)
  A_list <- dkge_cluster_loadings(fit)
  expect_length(Z_list, 4)
  expect_length(A_list, 4)
  expect_true(all(vapply(Z_list, ncol, integer(1)) == 2))
  expect_true(all(vapply(A_list, ncol, integer(1)) == 2))
  # Projection and loadings are consistent for identity K
  beta <- c(0.5, -1.2)
  decoded <- as.numeric(A_list[[1]] %*% beta)
  manual <- as.numeric(crossprod(fit$Btil[[1]], fit$K %*% fit$U %*% beta))
  expect_equal(decoded, manual, tolerance = 1e-10)
})


test_that("cross-fitted latent classifier produces per-subject betas", {
  fit <- make_demo_fit()
  y <- factor(c("case", "case", "control", "control"))
  folds <- structure(list(k = 2, assignments = list(c(1L, 3L), c(2L, 4L))),
                     class = "dkge_folds")
  clf <- dkge_cv_train_latent_classifier(fit, y,
                                         folds = folds,
                                         model = "lda",
                                         standardize = TRUE)
  expect_s3_class(clf, "dkge_clf")
  expect_length(clf$beta_by_subject, 4)
  expect_true(all(vapply(clf$beta_by_subject, length, integer(1)) == 2))
  expect_equal(clf$fold_assignment, c(1L, 2L, 1L, 2L))
  expect_true(all(vapply(clf$models_by_fold, function(m) m$model, character(1)) == "lda"))
})


test_that("decoder map matches manual aggregation", {
  fit <- make_demo_fit()
  renderer <- make_demo_renderer(fit)
  betas <- list(
    c(1.0, 0.4),
    c(0.6, 0.2),
    c(-0.5, -0.1),
    c(-0.7, -0.2)
  )
  info <- dkge_info_map_from_classifier(fit, betas, renderer,
                                        lambda = 0,
                                        to_vox = FALSE,
                                        inference = "none")
  A_list <- dkge_cluster_loadings(fit)
  subj_cluster <- lapply(seq_along(betas), function(s) as.numeric(A_list[[s]] %*% betas[[s]]))
  subj_anchor <- lapply(seq_along(subj_cluster), function(s) {
    apply_mapper(renderer$mapper_fits[[s]], subj_cluster[[s]])
  })
  manual <- dkge_anchor_aggregate(subj_anchor, subj_weights = renderer$weights)
  expect_equal(info$mean_anchor, manual$y, tolerance = 1e-10)
  expect_null(info$t_anchor)
  expect_equal(info$meta$kind, "decoder")
})


test_that("Haufe map collapses to decoder when covariance is identity", {
  fit <- make_demo_fit()
  renderer <- make_demo_renderer(fit)
  betas <- list(
    c(0.7, -0.2),
    c(0.5, -0.1),
    c(-0.6, 0.1),
    c(-0.4, 0.2)
  )
  clf <- make_demo_classifier(betas)
  dec <- dkge_info_map_from_classifier(fit, betas, renderer,
                                       lambda = 0,
                                       to_vox = FALSE,
                                       inference = "none")
  haufe <- dkge_info_map_haufe(fit, clf, renderer,
                               lambda = 0,
                               to_vox = FALSE,
                               inference = "none")
  expect_equal(haufe$mean_anchor, dec$mean_anchor, tolerance = 1e-6)
  expect_equal(haufe$meta$kind, "haufe")
})


test_that("LOCO proxy returns non-negative scores of correct length", {
  fit <- make_demo_fit()
  renderer <- make_demo_renderer(fit)
  betas <- list(
    c(1.0, 0.0),
    c(0.8, 0.1),
    c(-0.6, -0.1),
    c(-0.5, -0.2)
  )
  clf <- make_demo_classifier(betas)
  loco <- dkge_info_map_loco(fit, clf, renderer,
                             aggregate = "mean")
  expect_length(loco$loco_anchor, nrow(renderer$anchors))
  expect_true(all(loco$loco_anchor >= 0))
  expect_equal(loco$meta$kind, "loco_zero")
})


test_that("decoder accepts a single shared beta vector", {
  fit <- make_demo_fit()
  renderer <- make_demo_renderer(fit)
  beta_vec <- c(0.3, -0.7)
  info <- dkge_info_map_from_classifier(fit, beta_vec, renderer,
                                        lambda = 0,
                                        to_vox = FALSE,
                                        inference = "none")
  expect_length(info$subj_anchor, 4)
  expect_equal(info$meta$kind, "decoder")
})
