# test-as-data-frame.R

library(testthat)

test_that("as.data.frame.dkge_contrasts returns tidy rows", {
  fixture <- make_small_fit()
  contrasts <- diag(fixture$q)[, 1:2, drop = FALSE]
  contrast_obj <- dkge_contrast(fixture$fit, contrasts, method = "analytic",
                                parallel = FALSE, align = FALSE)

  df <- as.data.frame(contrast_obj)

  expected_rows <- length(contrast_obj$values) *
    length(contrast_obj$values[[1]]) *
    length(contrast_obj$values[[1]][[1]])

  expect_equal(nrow(df), expected_rows)
  expect_true(all(c("contrast", "subject", "component", "value", "method") %in%
                  names(df)))
  expected_contrast_names <- colnames(contrasts)
  if (is.null(expected_contrast_names) || any(!nzchar(expected_contrast_names))) {
    expected_contrast_names <- paste0("contrast", seq_len(ncol(contrasts)))
  }
  expect_setequal(unique(df$contrast), expected_contrast_names)
})

test_that("as.data.frame.dkge_inference summarises statistics", {
  fixture <- make_small_fit()
  contrasts <- diag(fixture$q)[, 1:2, drop = FALSE]
  infer_obj <- dkge_infer(fixture$fit, contrasts,
                          method = "analytic",
                          inference = "parametric",
                          correction = "none",
                          n_perm = 10)

  df <- as.data.frame(infer_obj)

  expect_true(all(c("contrast", "component", "statistic",
                    "p_value", "p_adjusted", "significant",
                    "alpha", "method", "inference", "correction") %in%
                  names(df)))
  expect_equal(length(unique(df$contrast)), ncol(contrasts))
  expect_false(any(is.na(df$statistic)))
})
