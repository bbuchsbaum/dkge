library(testthat)

test_that("dkge_kernel_prescreen orders kernels by alignment", {
  skip_on_cran()
  K_true <- diag(3)
  K_shift <- diag(c(0.8, 0.7, 0.6))
  K_bad <- matrix(c(1, 0.8, 0.8,
                    0.8, 1, -0.2,
                    0.8, -0.2, 1), nrow = 3)
  K_list <- list(true = K_true, shift = K_shift, bad = K_bad)

  C <- diag(3)
  tab <- dkge_kernel_prescreen(K_list,
                               C = C,
                               normalize_k = FALSE,
                               top_k = 3)

  expect_equal(tab$kernel[1], "true")
  expect_setequal(tab$kernel, names(K_list))
  expect_true(all(diff(tab$align) <= 0))
})
