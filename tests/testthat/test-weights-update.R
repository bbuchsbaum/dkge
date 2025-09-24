test_that("dkge_update_weights rebuilds bases with new weights", {
  wts1 <- dkge_weights(
    adapt = "kenergy", combine = "product", mix = 1,
    shrink = list(alpha = 1, winsor = 0.9999, normalize = "mean")
  )

  fit <- toy_real_fit(weights = wts1)
  expect_equal(fit$weight_spec$adapt, "kenergy")

  wts2 <- dkge_weights(
    adapt = "precision", combine = "product", mix = 1,
    shrink = list(alpha = 1, winsor = 0.9999, normalize = "mean")
  )

  fit2 <- dkge_update_weights(fit, wts2)
  expect_s3_class(fit2$weight_spec, "dkge_weights")
  expect_equal(fit2$weight_spec$adapt, "precision")

  w1 <- fit$voxel_weights
  w2 <- fit2$voxel_weights
  expect_equal(mean(w2), 1, tolerance = 1e-12)
  expect_gt(mean(abs(w1 - w2)), 1e-4)
})
