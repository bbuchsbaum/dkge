testthat::local_edition(3)
library(dkge)

test_that('.dkge_leading_sv_squared approximates leading singular value', {
  set.seed(321)
  X <- matrix(rnorm(25), 5, 5)
  approx <- dkge:::.dkge_leading_sv_squared(X, tol = 1e-8, max_iter = 100)
  exact <- svd(X, nu = 0, nv = 0)$d[1]^2
  expect_equal(approx, exact, tolerance = 1e-4)
})

test_that(".dkge_fit_prepare harmonises inputs", {
  S <- 3; q <- 4; P <- 5; T <- 18
  set.seed(123)
  betas <- replicate(S, matrix(rnorm(q * P), q, P), simplify = FALSE)
  designs <- replicate(S, {
    X <- matrix(rnorm(T * q), T, q)
    qr.Q(qr(X))
  }, simplify = FALSE)

  prepped <- dkge:::.dkge_fit_prepare(betas,
                                      designs = designs,
                                      K = diag(q),
                                      Omega_list = NULL,
                                      weights = NULL,
                                      rank = 3)

  expect_equal(prepped$dataset$n_subjects, S)
  expect_equal(length(prepped$Btil), S)
  expect_equal(prepped$rank, 3)
  expect_equal(dim(prepped$K), c(q, q))
  expect_equal(length(prepped$weight_eval), 6)
})

test_that(".dkge_fit_accumulate returns symmetric Chat", {
  S <- 3; q <- 4; P <- 5; T <- 16
  set.seed(99)
  betas <- replicate(S, matrix(rnorm(q * P), q, P), simplify = FALSE)
  designs <- replicate(S, {
    X <- matrix(rnorm(T * q), T, q)
    qr.Q(qr(X))
  }, simplify = FALSE)

  prepped <- dkge:::.dkge_fit_prepare(betas,
                                      designs = designs,
                                      K = diag(q),
                                      Omega_list = NULL,
                                      weights = NULL,
                                      rank = 3)
  accum <- dkge:::.dkge_fit_accumulate(prepped,
                                       w_method = "none",
                                       w_tau = 0)

  expect_true(is.matrix(accum$Chat))
  expect_equal(length(accum$contribs), S)
  expect_true(isSymmetric(unname(accum$Chat)))
  expect_equal(length(accum$subject_weights), S)
})

test_that(".dkge_fit_solve produces orthonormal basis", {
  S <- 3; q <- 4; P <- 5; T <- 20
  set.seed(456)
  betas <- replicate(S, matrix(rnorm(q * P), q, P), simplify = FALSE)
  designs <- replicate(S, {
    X <- matrix(rnorm(T * q), T, q)
    qr.Q(qr(X))
  }, simplify = FALSE)

  prepped <- dkge:::.dkge_fit_prepare(betas,
                                      designs = designs,
                                      K = diag(q),
                                      Omega_list = NULL,
                                      weights = NULL,
                                      rank = 3)
  accum <- dkge:::.dkge_fit_accumulate(prepped,
                                       w_method = "none",
                                       w_tau = 0)
  solved <- dkge:::.dkge_fit_solve(prepped,
                                   accum,
                                   rank = prepped$rank,
                                   cpca_part = "none",
                                   cpca_blocks = NULL,
                                   cpca_T = NULL,
                                   cpca_ridge = 0,
                                   ridge = 0)

  expect_equal(ncol(solved$U), prepped$rank)
  expect_equal(length(solved$sdev), prepped$rank)
  expect_true(all(solved$sdev >= 0))
  expect_equal(nrow(solved$U), q)
})

test_that(".dkge_fit_assemble returns dkge object", {
  S <- 3; q <- 4; P <- 6; T <- 15
  set.seed(2024)
  betas <- replicate(S, matrix(rnorm(q * P), q, P), simplify = FALSE)
  designs <- replicate(S, {
    X <- matrix(rnorm(T * q), T, q)
    qr.Q(qr(X))
  }, simplify = FALSE)

  prepped <- dkge:::.dkge_fit_prepare(betas,
                                      designs = designs,
                                      K = diag(q),
                                      Omega_list = NULL,
                                      weights = NULL,
                                      rank = 3)
  accum <- dkge:::.dkge_fit_accumulate(prepped,
                                       w_method = "none",
                                       w_tau = 0)
  solved <- dkge:::.dkge_fit_solve(prepped,
                                   accum,
                                   rank = prepped$rank,
                                   cpca_part = "none",
                                   cpca_blocks = NULL,
                                   cpca_T = NULL,
                                   cpca_ridge = 0,
                                   ridge = 0)
  fit <- dkge:::.dkge_fit_assemble(prepped,
                                   accum,
                                   solved,
                                   keep_X = TRUE,
                                   w_method = "none",
                                   w_tau = 0,
                                   ridge = 0)

  expect_s3_class(fit, "dkge")
  expect_equal(length(fit$weights), S)
  expect_equal(fit$rank, solved$rank)
  expect_equal(dim(fit$K), c(q, q))
  expect_true(is.matrix(fit$Chat_sym))
  expect_equal(fit$ridge_input, 0)
  expect_true(!is.null(fit$X_concat))
})
