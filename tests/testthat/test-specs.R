testthat::local_edition(3)
library(dkge)

test_that("dkge_transport_spec validates medoid and shape", {
  set.seed(1)
  centroids <- list(matrix(runif(12), 4, 3), matrix(runif(15), 5, 3))
  spec <- dkge_transport_spec(centroids = centroids, medoid = 2, lambda_size = 0.1)
  expect_s3_class(spec, "dkge_transport_spec")
  expect_equal(spec$medoid, 2L)
  expect_equal(spec$method, "sinkhorn")
  expect_error(dkge_transport_spec(list(), medoid = 1), "length")
})

test_that("dkge_inference_spec records arguments", {
  spec <- dkge_inference_spec(B = 500, tail = "left", center = "median")
  expect_s3_class(spec, "dkge_inference_spec")
  expect_equal(spec$B, 500L)
  expect_equal(spec$tail, "left")
  expect_equal(spec$center, "median")
})

test_that("dkge_classification_spec stores metadata", {
  spec <- expect_no_warning(dkge_classification_spec(targets = ~ condition,
                                                     method = "logit",
                                                     metric = c("accuracy", "logloss"),
                                                     mode = "cell",
                                                     n_perm = 10))
  expect_s3_class(spec, "dkge_classification_spec")
  expect_equal(spec$method, "logit")
  expect_true("n_perm" %in% names(spec))
  expect_equal(spec$mode, "cell")
})

test_that("dkge_predict_subjects handles list inputs", {
  fixtures <- make_small_fit(S = 3, q = 3, P = 4, T = 10, rank = 2, seed = 7)
  contrasts <- list(c1 = rep(0, fixtures$q))
  contrasts$c1[1] <- 1
  pred <- dkge_predict_subjects(fixtures$fit,
                                 betas = fixtures$betas,
                                 contrasts = contrasts,
                                 ids = paste0("s", seq_along(fixtures$betas)))
  expect_equal(length(pred$values), length(fixtures$betas))
  expect_true(all(names(pred$values) == paste0("s", seq_along(fixtures$betas))))
})

test_that("dkge_predict_subjects accepts dkge_data", {
  fixtures <- make_small_fit(S = 2, q = 4, P = 3, T = 8, rank = 2, seed = 11)
  subjects <- lapply(seq_along(fixtures$betas), function(i) {
    dkge_subject(fixtures$betas[[i]], design = diag(fixtures$q), id = paste0("d", i))
  })
  data_bundle <- dkge_data(subjects)
  contrasts <- matrix(c(1, rep(0, fixtures$q - 1)), fixtures$q, 1)
  pred <- dkge_predict_subjects(fixtures$fit,
                                 betas = data_bundle,
                                 contrasts = contrasts,
                                 return_loadings = FALSE)
  expect_equal(length(pred$values), length(fixtures$betas))
  expect_true(all(names(pred$values) == data_bundle$subject_ids))
  expect_null(pred$A_list)
})
