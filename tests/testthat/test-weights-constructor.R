test_that("dkge_weights constructor/print behave", {
  wts <- dkge_weights_auto()
  expect_s3_class(wts, "dkge_weights")
  expect_equal(wts$adapt, "kenergy_prec")
  expect_equal(wts$combine, "product")
  expect_true(wts$mix >= 0 && wts$mix <= 1)

  expect_output(print(wts), "dkge weight specification")
  expect_output(print(wts), "adapt")
  expect_output(print(wts), "combine")
  expect_output(print(wts), "shrink")
})

test_that("prior sugar helpers build numeric vectors", {
  V <- 10L
  mask <- c(rep(TRUE, 3), rep(FALSE, V - 3))
  p1 <- dkge_weights_prior_mask(mask, value_in = 2, value_out = 0)
  expect_type(p1, "double")
  expect_equal(sum(p1 > 0), 3)

  labs <- factor(rep(1:2, each = V / 2))
  p2 <- dkge_weights_prior_roi(labs, roi_values = c("1" = 2, "2" = 0.5))
  expect_length(p2, V)
  expect_true(all(p2[labs == 1] == 2))
  expect_true(all(p2[labs == 2] == 0.5))
})
