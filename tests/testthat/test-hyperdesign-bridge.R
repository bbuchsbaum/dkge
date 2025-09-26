test_that("as_dkge_kernel preserves existing matrix behaviour", {
  set.seed(1)
  q <- 3
  P <- 4
  S <- 2

  betas <- replicate(S, matrix(rnorm(q * P), q, P), simplify = FALSE)
  designs <- replicate(S, diag(q), simplify = FALSE)

  kernel_obj <- list(K = diag(q), info = list(tag = "ok"))
  fit_obj <- dkge(betas, designs = designs, kernel = kernel_obj, subject_ids = paste0("s", seq_len(S)))

  expect_equal(fit_obj$kernel_info$tag, "ok")

  fit_default <- dkge(betas, designs = designs, kernel = diag(q), subject_ids = paste0("s", seq_len(S)))
  expect_null(fit_default$kernel_info)
})

test_that("as_dkge_folds accepts subject-fold data.frame", {
  fixtures <- make_small_fit()
  fit <- fixtures$fit

  df <- data.frame(
    subject = fit$subject_ids,
    fold = c("train", "test", "train"),
    stringsAsFactors = FALSE
  )

  folds <- as_dkge_folds(df, fit)
  expect_s3_class(folds, "dkge_folds")
  expect_equal(length(folds$assignments), 2L)
  expect_equal(sort(unname(unlist(folds$assignments))), 1:3)
})

test_that("dkge_contrast accepts folds coercible via as_dkge_folds", {
  fixtures <- make_small_fit()
  fit <- fixtures$fit
  df <- data.frame(
    subject = fit$subject_ids,
    fold = rep(1:3, each = 1L),
    stringsAsFactors = FALSE
  )

  result <- dkge_contrast(fit, contrasts = rep(1 / fixtures$q, fixtures$q), method = "kfold", folds = df)
  expect_s3_class(result, "dkge_contrasts")
  expect_named(result$values)
})
