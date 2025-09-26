library(testthat)
library(dkge)

make_plot_fit <- function(S = 3, P = 4, seed = 1) {
  set.seed(seed)
  factors <- list(A = list(L = 2), B = list(L = 2), time = list(L = 4))
  dk <- design_kernel(factors, basis = "effect")
  q <- nrow(dk$K)
  betas <- replicate(S, matrix(rnorm(q * P), q, P), simplify = FALSE)
  designs <- replicate(S, diag(q), simplify = FALSE)
  fit <- dkge_fit(betas, designs, K = dk, rank = min(3, q))
  list(fit = fit, kernel = dk)
}

test_that("plot suite falls back gracefully when patchwork missing", {
  skip_if_not_installed("ggplot2")
  fixture <- make_plot_fit()
  fit <- fixture$fit
  orig_requireNamespace <- base::requireNamespace
  plot_obj <- with_mocked_bindings(
    dkge_plot_suite(fit, save_path = NULL),
    .package = "base",
    requireNamespace = function(pkg, ...) {
      if (pkg == "patchwork") return(FALSE)
      orig_requireNamespace(pkg, ...)
    }
  )
  expect_s3_class(plot_obj, "ggplot")
})
