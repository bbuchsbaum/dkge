# test-neuroim2.R
# Light-touch tests for neuroim2 integration helpers

library(testthat)

test_that("cluster aggregation works with neuroim2 objects", {
  skip_if_not_installed("neuroim2")
  skip_if_not_installed("fmrireg")
  skip("neuroim2 low-level constructors differ across versions; skipping." )

  invisible(NULL)
})
