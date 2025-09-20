# test-mapper.R
# Direct tests for mapper specification, fitting, and prediction.

library(testthat)

test_that("sinkhorn mapper produces balanced transport plan", {
  source_feat <- matrix(c(0, 0,
                          1, 0), nrow = 2, byrow = TRUE)
  target_feat <- matrix(c(0, 0,
                          2, 0), nrow = 2, byrow = TRUE)

  spec <- dkge_mapper_spec("sinkhorn", epsilon = 0.1, lambda_spa = 0)
  mapping <- fit_mapper(spec, source_feat = source_feat, target_feat = target_feat)

  plan <- mapping$operator
  expect_equal(rowSums(plan), rep(0.5, 2), tolerance = 1e-3)
  expect_equal(colSums(plan), rep(0.5, 2), tolerance = 1e-3)

  values <- c(3, 7)
  mapped <- predict_mapper(mapping, values)
  expect_equal(mapped, as.numeric(t(plan) %*% values), tolerance = 1e-8)
})

test_that("ridge mapper matches closed-form solution", {
  source_feat <- diag(3)
  target_feat <- matrix(c(1, 2, 0,
                          0, 1, 3,
                          2, 0, 1), nrow = 3, byrow = TRUE)
  spec <- dkge_mapper_spec("ridge", lambda = 0)
  mapping <- fit_mapper(spec, source_feat = source_feat, target_feat = target_feat)

  expect_equal(mapping$operator, t(target_feat), tolerance = 1e-10)

  vals <- c(5, 4, -1)
  expect_equal(predict_mapper(mapping, vals), as.numeric(target_feat %*% vals), tolerance = 1e-10)
})

test_that("ridge mapper applies regularisation", {
  source_feat <- diag(2)
  target_feat <- matrix(c(2, 0,
                          0, 2), nrow = 2, byrow = TRUE)
  spec <- dkge_mapper_spec("ridge", lambda = 1)
  mapping <- fit_mapper(spec, source_feat = source_feat, target_feat = target_feat)
  expect_equal(mapping$operator, diag(2 / (1 + 1), 2), tolerance = 1e-10)
})

test_that("sinkhorn mapper solves near-identity transport when supports align", {
  subj_points <- matrix(c(0, 0, 0,
                          1, 0, 0), ncol = 3, byrow = TRUE)
  anchor_points <- subj_points
  mapper <- dkge_mapper("sinkhorn",
                        epsilon = 1e-3,
                        sigx = 1,
                        lambda_xyz = 1)
  fit <- fit_mapper(mapper,
                    subj_points = subj_points,
                    anchor_points = anchor_points)
  plan <- fit$plan
  expect_equal(Matrix::rowSums(plan), rep(0.5, 2), tolerance = 1e-6)
  expect_equal(Matrix::colSums(plan), rep(0.5, 2), tolerance = 1e-6)

  vals <- c(-2, 4)
  mapped <- apply_mapper(fit, vals)
  expect_equal(mapped, vals, tolerance = 1e-2)
})


test_that("sinkhorn mapper honours reliability weights", {
  subj_points <- matrix(c(0, 0, 0,
                          2, 0, 0,
                          4, 0, 0), ncol = 3, byrow = TRUE)
  anchor_points <- cbind(seq(0, 4, length.out = 5), 0, 0)
  mapper <- dkge_mapper("sinkhorn",
                        epsilon = 5e-3,
                        sigx = 2,
                        lambda_xyz = 1)
  reliab <- c(1, 0.5, 0.25)
  fit <- fit_mapper(mapper,
                    subj_points = subj_points,
                    anchor_points = anchor_points,
                    reliab = reliab)

  values <- c(3, 6, 9)
  mapped <- apply_mapper(fit, values, normalize_by_reliab = TRUE)
  expect_true(all(is.finite(mapped)))
  expect_gt(max(mapped), min(mapped))
})

test_that("sinkhorn mapper leverages features when geometry is ambiguous", {
  subj_points <- matrix(0, nrow = 2, ncol = 3)
  anchor_points <- matrix(0, nrow = 2, ncol = 3)
  subj_feats <- rbind(c(1, 0), c(0, 1))
  anchor_feats <- rbind(c(1, 0), c(0, 1))

  mapper <- dkge_mapper("sinkhorn",
                        epsilon = 1e-3,
                        lambda_xyz = 0,
                        lambda_feat = 1,
                        sigz = 1)
  fit <- fit_mapper(mapper,
                    subj_points = subj_points,
                    anchor_points = anchor_points,
                    subj_feats = subj_feats,
                    anchor_feats = anchor_feats)

  plan <- as.matrix(fit$plan)
  expect_gt(plan[1, 1], 0.45)
  expect_gt(plan[2, 2], 0.45)
  expect_lt(plan[1, 2], 0.05)
  expect_lt(plan[2, 1], 0.05)
})
