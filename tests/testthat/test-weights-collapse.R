test_that("k-energy kernel resolves with time collapse correctly", {
  kinf <- toy_kernel_info()

  wts <- dkge_weights(
    adapt   = "kenergy",
    combine = "product",
    mix     = 0.6,
    collapse = list(time = list(method = "mean", window = 2:3))
  )

  Keff <- dkge:::.dkge_weights_resolve_k(wts, kinf)

  IA <- diag(2); IB <- diag(2)
  w  <- c(0, 0.5, 0.5)
  Mtime <- outer(w, w, "*")
  Kcell <- kronecker(kronecker(IA, IB), Mtime)

  expect_equal(Keff, 0.5 * (Keff + t(Keff)), tolerance = 1e-12)
  expect_equal(Keff, Kcell, tolerance = 1e-12)
})
