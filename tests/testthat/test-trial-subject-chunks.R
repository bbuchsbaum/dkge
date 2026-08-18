library(testthat)

test_that("chunked trial reduction matches the dense OLS constructor", {
  set.seed(19401)
  cell <- rep(1:3, each = 8)
  X <- model.matrix(~ 0 + factor(cell, levels = 1:3))
  colnames(X) <- paste0("e", 1:3)
  truth <- matrix(rnorm(3 * 9), 3, 9)
  Y <- X %*% truth + matrix(rnorm(nrow(X) * 9, sd = 0.4), nrow(X), 9)
  colnames(Y) <- paste0("v", 1:9)

  dense <- dkge_trial_subject(Y, X)
  chunked <- dkge_trial_subject_chunks(
    list(Y[, 1:2], Y[, 3:6], Y[, 7:9]), X
  )
  expect_equal(chunked$beta, dense$beta, tolerance = 1e-12,
               ignore_attr = TRUE)
  expect_equal(chunked$effect_noise_cov, dense$effect_noise_cov,
               tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(chunked$residual_variance, dense$residual_variance,
               tolerance = 1e-12, ignore_attr = TRUE)
  expect_null(chunked$noise_trace)
  expect_equal(dkge:::.dkge_noise_trace(chunked),
               sum(dense$residual_variance), tolerance = 1e-12)
  expect_equal(chunked$input_type, "trialwise_chunks")
  expect_equal(chunked$error_model$n_chunks, 3L)
  expect_null(chunked$effect_score)
})

test_that("chunked split halves match the dense constructor", {
  set.seed(19402)
  X <- model.matrix(~ 0 + factor(rep(1:2, each = 8)))
  colnames(X) <- c("e1", "e2")
  Y <- matrix(rnorm(nrow(X) * 8), nrow(X), 8)
  colnames(Y) <- paste0("v", 1:8)

  dense <- dkge_trial_subject(Y, X, split = "within_cell")
  chunked <- dkge_trial_subject_chunks(
    list(Y[, 1:3], Y[, 4:8]), X, split = "within_cell"
  )
  expect_equal(chunked$split_betas, dense$split_betas,
               tolerance = 1e-12, ignore_attr = TRUE)
})

test_that("function chunk sources terminate without retaining trialwise Y", {
  set.seed(19403)
  X <- model.matrix(~ 0 + factor(rep(1:2, each = 8)))
  colnames(X) <- c("e1", "e2")
  blocks <- lapply(c(3, 4, 2), function(width) {
    matrix(rnorm(nrow(X) * width), nrow(X), width)
  })
  names_seen <- 0L
  source <- function(i) {
    if (i > length(blocks)) return(NULL)
    block <- blocks[[i]]
    colnames(block) <- paste0("v", names_seen + seq_len(ncol(block)))
    names_seen <<- names_seen + ncol(block)
    block
  }
  subject <- dkge_trial_subject_chunks(source, X)
  expect_equal(dim(subject$beta), c(2L, 9L))
  expect_null(subject$y)
  expect_null(subject$trialwise_y)
  expect_equal(subject$error_model$storage, "feature_chunked")
  expect_equal(object.size(subject$design), object.size(X))
})

test_that("chunk sources fail clearly on malformed input", {
  X <- matrix(1, 4, 1, dimnames = list(NULL, "mean"))
  expect_error(dkge_trial_subject_chunks(list(), X), "no feature blocks")
  expect_error(
    dkge_trial_subject_chunks(function(i) matrix(1, 4, 1), X,
                              max_chunks = 2),
    "max_chunks"
  )
  expect_error(
    dkge_trial_subject_chunks(list(matrix(1, 3, 2)), X),
    "T x P_block"
  )
  duplicate_names <- list(
    matrix(1, 4, 1, dimnames = list(NULL, "v1")),
    matrix(2, 4, 1, dimnames = list(NULL, "v1"))
  )
  expect_error(dkge_trial_subject_chunks(duplicate_names, X), "unique")
})

test_that("chunked analytic debias matches dense under spatial and voxel weights", {
  set.seed(19404)
  cell <- factor(rep(1:2, each = 6), levels = 1:2)
  X <- model.matrix(~ 0 + cell)
  colnames(X) <- c("e1", "e2")
  omega <- c(0.2, 1, 3, 0.5)
  make_Y <- function() {
    truth <- matrix(rnorm(8), 2, 4)
    Y <- X %*% truth + matrix(rnorm(nrow(X) * 4, sd = 0.35), nrow(X), 4)
    colnames(Y) <- paste0("v", 1:4)
    Y
  }
  Ys <- list(make_Y(), make_Y())
  dense <- lapply(seq_along(Ys), function(i) {
    dkge_trial_subject(Ys[[i]], X, id = paste0("s", i), omega = omega)
  })
  chunked <- lapply(seq_along(Ys), function(i) {
    dkge_trial_subject_chunks(
      list(Ys[[i]][, 1:2], Ys[[i]][, 3:4]), X,
      id = paste0("s", i), omega = omega
    )
  })

  for (i in seq_along(Ys)) {
    oracle <- sum(omega * dense[[i]]$residual_variance)
    expect_equal(dkge:::.dkge_noise_trace(dense[[i]], Omega = omega),
                 oracle, tolerance = 1e-12)
    expect_equal(dkge:::.dkge_noise_trace(chunked[[i]], Omega = omega),
                 oracle, tolerance = 1e-12)
    expect_null(chunked[[i]]$noise_trace)
  }

  K <- diag(2)
  dimnames(K) <- list(colnames(X), colnames(X))
  fit_dense <- dkge_fit(dkge_data(dense), K = K, rank = 1, w_method = "none",
                        effect_scaling = "none", debias = "analytic")
  fit_chunked <- dkge_fit(dkge_data(chunked), K = K, rank = 1, w_method = "none",
                          effect_scaling = "none", debias = "analytic")
  expect_equal(fit_chunked$Chat, fit_dense$Chat, tolerance = 1e-12)

  prior <- c(0.5, 2, 0.25, 1.5)
  wspec <- dkge_weights(prior = prior, adapt = "none")
  fit_dense_w <- dkge_fit(dkge_data(dense), K = K, rank = 1, w_method = "none",
                          effect_scaling = "none", debias = "analytic",
                          weights = wspec)
  fit_chunked_w <- dkge_fit(dkge_data(chunked), K = K, rank = 1, w_method = "none",
                            effect_scaling = "none", debias = "analytic",
                            weights = wspec)
  expect_equal(fit_chunked_w$Chat, fit_dense_w$Chat, tolerance = 1e-12)
  expect_false(isTRUE(all.equal(fit_dense_w$Chat, fit_dense$Chat)))
})
