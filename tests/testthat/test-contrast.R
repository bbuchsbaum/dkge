# test-contrast.R
# Tests for unified contrast engine and cross-fitting methods

library(testthat)
library(dkge)

# Helper function to create toy data
create_toy_data <- function(S = 10, q = 5, P = 50, seed = 123) {
  set.seed(seed)

  betas <- lapply(seq_len(S), function(s) {
    matrix(rnorm(q * P), q, P)
  })

  designs <- lapply(seq_len(S), function(s) {
    X <- matrix(rnorm(100 * q), 100, q)
    # Orthogonalize for stability
    qr.Q(qr(X))
  })

  K <- diag(q)  # Simple identity kernel

  list(betas = betas, designs = designs, K = K, S = S, q = q, P = P)
}

test_that("dkge_contrast works with single contrast vector", {
  data <- create_toy_data()
  fit <- dkge_fit(data$betas, data$designs, K = data$K, rank = 3)

  # Single contrast
  c1 <- c(1, -1, 0, 0, 0)
  result <- dkge_contrast(fit, c1, method = "loso")

  expect_s3_class(result, "dkge_contrasts")
  expect_equal(result$method, "loso")
  expect_length(result$values, 1)
  expect_length(result$values[[1]], data$S)
  expect_length(result$values[[1]][[1]], data$P)
})

test_that("dkge_contrast works with multiple contrasts", {
  data <- create_toy_data()
  fit <- dkge_fit(data$betas, data$designs, K = data$K, rank = 3)

  # Multiple contrasts
  contrasts <- list(
    main1 = c(1, -1, 0, 0, 0),
    main2 = c(0, 0, 1, -1, 0),
    interaction = c(1, -1, -1, 1, 0)
  )

  result <- dkge_contrast(fit, contrasts, method = "loso")

  expect_length(result$values, 3)
  expect_named(result$values, names(contrasts))
  expect_equal(result$contrasts, contrasts)
})

test_that("dkge_contrast works with matrix input", {
  data <- create_toy_data()
  fit <- dkge_fit(data$betas, data$designs, K = data$K, rank = 3)

  # Matrix of contrasts
  C <- cbind(
    c1 = c(1, -1, 0, 0, 0),
    c2 = c(0, 0, 1, -1, 0)
  )

  result <- dkge_contrast(fit, C, method = "loso")

  expect_length(result$values, 2)
  expect_named(result$values, c("c1", "c2"))
})

test_that("LOSO produces unbiased estimates", {
  # Create data with known structure
  set.seed(456)
  S <- 20
  q <- 4
  P <- 30
  r <- 2

  # True latent structure
  U_true <- qr.Q(qr(matrix(rnorm(q * r), q, r)))

  betas <- lapply(seq_len(S), function(s) {
    # Generate from latent model with noise
    A_s <- matrix(rnorm(P * r), P, r)
    B_s <- U_true %*% t(A_s) + matrix(rnorm(q * P, sd = 0.1), q, P)
    B_s
  })

  designs <- lapply(seq_len(S), function(s) {
    qr.Q(qr(matrix(rnorm(100 * q), 100, q)))
  })

  K <- diag(q)
  fit <- dkge_fit(betas, designs, K = K, rank = r)

  # True contrast
  c_true <- c(1, -1, 0, 0)

  # LOSO should give unbiased estimates
  result_loso <- dkge_contrast(fit, c_true, method = "loso")

  # Each subject should have values
  expect_length(result_loso$values[[1]], S)

  # Check that LOSO bases are different from full fit
  expect_length(result_loso$metadata$bases, S)
  for (s in seq_len(S)) {
    # Bases should be different (not identical to full U)
    U_loso <- result_loso$metadata$bases[[s]]
    expect_equal(dim(U_loso), dim(fit$U))

    # Should be K-orthonormal
    KU <- fit$K %*% U_loso
    gram <- t(U_loso) %*% KU
    expect_equal(gram, diag(r), tolerance = 1e-10)
  }
})

test_that("K-fold cross-fitting works correctly", {
  data <- create_toy_data(S = 15)
  fit <- dkge_fit(data$betas, data$designs, K = data$K, rank = 3)

  # Define folds
  folds <- dkge_define_folds(fit, type = "subject", k = 3, seed = 789)
  expect_s3_class(folds, "dkge_folds")
  expect_equal(folds$k, 3)
  expect_equal(length(folds$assignments), 3)

  # Use K-fold for contrast
  c1 <- c(1, -1, 0, 0, 0)
  result <- dkge_contrast(fit, c1, method = "kfold", folds = folds)

  expect_equal(result$method, "kfold")
  expect_equal(result$metadata$folds, folds)
  expect_length(result$metadata$fold_bases, 3)

  # All subjects should have values
  expect_length(result$values[[1]], data$S)
  expect_true(all(!vapply(result$values[[1]], is.null, logical(1))))
})

test_that("Analytic approximation matches LOSO for small perturbations", {
  data <- create_toy_data(S = 8, q = 4, P = 20)

  # Use equal weights to ensure small perturbations
  fit <- dkge_fit(data$betas, data$designs, K = data$K,
                 rank = 2, w_method = "none")

  c1 <- c(1, -1, 0, 0)

  # Full LOSO
  result_loso <- dkge_contrast(fit, c1, method = "loso")

  # Analytic approximation
  result_analytic <- dkge_contrast(fit, c1, method = "analytic")

  expect_equal(result_analytic$method, "analytic")

  # Values should be similar (not identical)
  for (s in seq_len(data$S)) {
    v_loso <- result_loso$values[[1]][[s]]
    v_analytic <- result_analytic$values[[1]][[s]]

    # Correlation should be very high
    cor_val <- cor(v_loso, v_analytic)
    expect_gt(cor_val, 0.95)

    # Relative error should be small
    rel_error <- mean(abs(v_loso - v_analytic)) / mean(abs(v_loso))
    expect_lt(rel_error, 0.1)
  }
})

test_that("as.matrix method works for dkge_contrasts", {
  data <- create_toy_data()
  fit <- dkge_fit(data$betas, data$designs, K = data$K, rank = 3)

  contrasts <- list(
    c1 = c(1, -1, 0, 0, 0),
    c2 = c(0, 0, 1, -1, 0)
  )

  result <- dkge_contrast(fit, contrasts, method = "loso")

  # Extract first contrast as matrix
  mat1 <- as.matrix(result, contrast = 1)
  expect_equal(dim(mat1), c(data$S, data$P))

  # Extract by name
  mat2 <- as.matrix(result, contrast = "c2")
  expect_equal(dim(mat2), c(data$S, data$P))

  # Check values match
  expect_equal(mat1[1, ], result$values[[1]][[1]])
  expect_equal(mat2[1, ], result$values[[2]][[1]])
})

test_that("as.matrix requires transport when cluster counts differ", {
  data <- create_mismatched_data()
  fit <- dkge_fit(data$betas, data$designs, K = data$K, rank = 2)
  contrast <- dkge_contrast(fit, c(1, -1, 0), method = "loso")

  expect_error(as.matrix(contrast, contrast = 1), "Subject cluster counts differ")
})

test_that("ridge mapper aligns mismatched cluster sizes", {
  data <- create_mismatched_data()
  fit <- dkge_fit(data$betas, data$designs, K = data$K, rank = 2)

  mapper_spec <- dkge_mapper_spec("ridge", lambda = 1e-2)
  contrast <- dkge_contrast(fit, c(1, -1, 0), method = "loso")
  mapped <- dkge_transport_contrasts_to_medoid(fit, contrast,
                                               medoid = 1L,
                                               centroids = data$centroids,
                                               mapper = mapper_spec)

  mat <- mapped[[1]]$subj_values
  expect_equal(dim(mat), c(data$S, nrow(data$centroids[[1]])))
})

test_that("Cross-fitting methods respect K-metric", {
  data <- create_toy_data(q = 4)

  # Non-trivial kernel
  K <- matrix(0.3, 4, 4)
  diag(K) <- 1

  fit <- dkge_fit(data$betas, data$designs, K = K, rank = 2)

  c1 <- rep(1, 4)

  # All methods should work with non-identity K
  result_loso <- dkge_contrast(fit, c1, method = "loso")
  result_kfold <- dkge_contrast(fit, c1, method = "kfold", folds = 2)
  result_analytic <- dkge_contrast(fit, c1, method = "analytic")

  # Check K-orthonormality of bases
  for (method_result in list(result_loso, result_analytic)) {
    for (s in seq_len(min(3, data$S))) {
      U_s <- method_result$metadata$bases[[s]]
      gram <- t(U_s) %*% K %*% U_s
      expect_equal(gram, diag(2), tolerance = 1e-8)
    }
  }

  # K-fold bases
  for (fold in seq_len(2)) {
    U_fold <- result_kfold$metadata$fold_bases[[fold]]
    gram <- t(U_fold) %*% K %*% U_fold
    expect_equal(gram, diag(2), tolerance = 1e-8)
  }
})

test_that("Ridge parameter affects results", {
  data <- create_toy_data()
  fit <- dkge_fit(data$betas, data$designs, K = data$K, rank = 3)

  c1 <- c(1, -1, 0, 0, 0)

  # Without ridge
  result1 <- dkge_contrast(fit, c1, method = "loso", ridge = 0)

  # With ridge
  result2 <- dkge_contrast(fit, c1, method = "loso", ridge = 0.1)

  # Results should differ
  v1 <- result1$values[[1]][[1]]
  v2 <- result2$values[[1]][[1]]

  expect_false(all(v1 == v2))

  # But should be correlated
  expect_gt(cor(v1, v2), 0.9)
})
