test_that("canonical default MFA weights match the pre-extension baseline", {
  set.seed(20260820L)
  S <- 6L
  q <- 5L
  P <- 13L
  Tn <- 24L
  effects <- paste0("e", seq_len(q))
  betas <- replicate(S, {
    B <- matrix(stats::rnorm(q * P), q, P)
    rownames(B) <- effects
    B
  }, simplify = FALSE)
  designs <- replicate(S, {
    X <- matrix(stats::rnorm(Tn * q), Tn, q)
    colnames(X) <- effects
    X
  }, simplify = FALSE)
  A <- matrix(stats::rnorm(q * q), q, q)
  K <- crossprod(A) / q + diag(q)
  dimnames(K) <- list(effects, effects)

  set.seed(99)
  caller_rng <- .Random.seed
  fit <- dkge_fit(betas, designs, K = K, rank = 3L)
  expect_identical(.Random.seed, caller_rng)
  expect_identical(
    sprintf("%a", unname(fit$weights)),
    c(
      "0x1.092e59420ccfep+0", "0x1.b1938e2f0683ap-1",
      "0x1.9241e5f5f2ae8p-1", "0x1.07e262b545053p+0",
      "0x1.ff6c5984d258p-1", "0x1.4d4e5d33c865bp+0"
    )
  )

  set.seed(123456)
  fit_other_rng <- dkge_fit(betas, designs, K = K, rank = 3L)
  expect_identical(fit_other_rng$weights, fit$weights)
  expect_identical(fit_other_rng$Chat, fit$Chat)

  fit_none <- dkge_fit(
    betas, designs, K = K, rank = 3L, w_method = "none"
  )
  expect_identical(unname(fit_none$weights), rep(1, S))
})

test_that("stored full-coverage golden evidence satisfies every threshold", {
  evidence <- system.file(
    "extdata", "dkge-full-coverage-golden.csv", package = "dkge"
  )
  metadata <- system.file(
    "extdata", "dkge-full-coverage-golden-metadata.csv", package = "dkge"
  )
  expect_true(file.exists(evidence))
  expect_true(file.exists(metadata))
  if (!file.exists(evidence) || !file.exists(metadata)) return(invisible(NULL))

  results <- utils::read.csv(evidence, stringsAsFactors = FALSE)
  meta <- utils::read.csv(metadata, stringsAsFactors = FALSE)
  expect_setequal(results$path, c("default", "none"))
  expect_true(all(results$pass))
  expect_equal(
    results$threshold[results$path == "default" & results$artifact == "Chat"],
    1e-15
  )
  expect_equal(
    results$threshold[results$path == "default" & results$artifact == "loso"],
    1e-14
  )
  expect_identical(
    meta$source_sha[meta$role == "baseline"],
    "953d2a64ff356d55f4a6fd3278a7f11db018fea6"
  )
  expect_true(all(nzchar(meta$dependencies)))
  expect_true(all(nzchar(meta$command)))
})
