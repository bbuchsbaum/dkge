test_that("fold builder uses training-only adaptive weights and weighted G", {
  fit   <- toy_fold_fit(nsub = 3L, Q = 12L, V = 10L)
  folds <- fit$folds_index

  wts <- dkge_weights(
    prior   = NULL,
    adapt   = "kenergy",
    combine = "product",
    mix     = 1,
    shrink  = list(alpha = 1, winsor = 0.9999, normalize = "mean"),
    scope   = "fold"
  )

  fold_info <- dkge:::.dkge_build_fold_bases(fit, folds, weights = wts, ridge = 1e-6)
  f1 <- fold_info$folds[[1L]]

  expect_true(is.matrix(f1$U_minus))
  expect_true(is.numeric(f1$D_minus))
  expect_length(f1$weights$w_total, ncol(fit$Btil[[1L]]))

  holdout <- folds[[1L]]
  train_ids <- setdiff(seq_len(length(fit$Btil)), holdout)
  B_train <- fit$Btil[train_ids]

  exp_w_adapt <- Reduce("+", lapply(B_train, function(B) colSums(B^2))) / length(B_train)
  exp_w_adapt <- dkge:::.dkge_winsor(exp_w_adapt, upper = wts$shrink$winsor)
  exp_w_adapt <- dkge:::.dkge_norm_vec(exp_w_adapt, "mean")
  expect_equal(f1$weights$w_total, exp_w_adapt, tolerance = 1e-4)

  w <- f1$weights$w_total
  accum <- dkge:::.dkge_accumulate_chat(B_train, vector("list", length(B_train)),
                                        diag(nrow(B_train[[1L]])),
                                        rep(1, length(B_train)),
                                        voxel_weights = w)
  G_ref <- accum$Chat
  G_ref <- G_ref + (1e-6 * sum(diag(G_ref)) / nrow(G_ref)) * diag(nrow(G_ref))
  G_ref <- 0.5 * (G_ref + t(G_ref))
  G_rec <- f1$U_minus %*% diag(f1$D_minus, ncol(f1$U_minus)) %*% t(f1$U_minus)
  expect_equal(G_rec, G_ref, tolerance = 1e-6)

  expect_gte(min(f1$D_minus), -1e-10)
})

test_that("time collapse in weights integrates with fold builder", {
  fit   <- toy_fold_fit(nsub = 3L, Q = 12L, V = 10L)
  folds <- fit$folds_index

  wts <- dkge_weights(
    prior   = NULL,
    adapt   = "kenergy",
    combine = "product",
    mix     = 1,
    shrink  = list(alpha = 1, winsor = 0.9999, normalize = "mean"),
    scope   = "fold",
    collapse = list(time = list(method = "mean", window = 2:3))
  )

  fold_info <- dkge:::.dkge_build_fold_bases(fit, folds, weights = wts, ridge = 1e-6)
  expect_true(is.list(fold_info$folds[[1]]$weights))
  expect_equal(mean(fold_info$folds[[1]]$weights$w_total), 1, tolerance = 1e-12)
})
