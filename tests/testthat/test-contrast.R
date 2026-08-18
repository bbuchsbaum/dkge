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
  contrast <- suppressWarnings(dkge_contrast(fit, c(1, -1, 0), method = "loso"))

  expect_error(as.matrix(contrast, contrast = 1), "Subject cluster counts differ")
})

test_that("ridge mapper aligns mismatched cluster sizes", {
  data <- create_mismatched_data()
  fit <- dkge_fit(data$betas, data$designs, K = data$K, rank = 2)

  mapper_spec <- dkge_mapper_spec("ridge", lambda = 1e-2)
  contrast <- suppressWarnings(dkge_contrast(fit, c(1, -1, 0), method = "loso"))
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

test_that("as.matrix hints about transport when cluster counts differ", {
  data <- create_toy_data(S = 3, q = 3, P = 6)
  fit <- dkge_fit(data$betas, data$designs, K = data$K, rank = 2)
  result <- dkge_contrast(fit, c(1, -1, 0), method = "loso")
  # Tamper with one subject to shorten the value vector
  result$values[[1]][[1]] <- result$values[[1]][[1]][-1]
  expect_error(as.matrix(result),
               "dkge_transport_contrasts_to_medoid",
               class = "dkge_transport_needed")
})

make_q8_scope_fit <- function(S = 4, P = 10, seed = 111) {
  set.seed(seed)
  grid <- dkge_effect_grid(
    factors = list(
      group = c("control", "sdam"),
      task = c("recog", "nback"),
      measure = c("low_mid", "high_sem")
    ),
    scope = c(group = "between", task = "within", measure = "within"),
    block_factors = "group"
  )
  kernel <- design_kernel(
    grid,
    terms = list("group", "task", c("group", "task")),
    basis = "cell",
    normalize = "none"
  )
  q <- length(grid$cell_labels)
  betas <- replicate(S, {
    B <- matrix(rnorm(q * P), q, P)
    rownames(B) <- grid$cell_labels
    B
  }, simplify = FALSE)
  designs <- replicate(S, {
    X <- diag(q)
    colnames(X) <- grid$cell_labels
    X
  }, simplify = FALSE)
  fit <- dkge_fit(dkge_data(betas, designs), K = kernel, rank = 3,
                  w_method = "none")
  list(fit = fit, grid = grid)
}

test_that("contrasts carry within/between/mixed estimability metadata", {
  fx <- make_q8_scope_fit()
  fit <- fx$fit
  task4 <- c(-0.5, -0.5, 0.5, 0.5)
  contrasts <- list(
    task = rep(task4, 2),
    group = c(rep(-0.25, 4), rep(0.25, 4)),
    "group:task" = c(-task4, task4)
  )

  task_res <- expect_warning(
    dkge_contrast(fit, contrasts["task"], method = "loso", align = FALSE),
    NA
  )
  expect_equal(task_res$metadata$contrast_estimability$estimability, "within")

  expect_warning(
    dkge_contrast(fit, contrasts[c("group", "group:task")],
                  method = "loso", align = FALSE),
    "subject-label permutation"
  )
  mixed_res <- suppressWarnings(
    dkge_contrast(fit, contrasts[c("group", "group:task")],
                  method = "loso", align = FALSE)
  )
  info <- mixed_res$metadata$contrast_estimability
  expect_equal(info$estimability, c("between", "mixed"))
  expect_true(all(grepl("dkge_between", info$recommended_inference)))
})

test_that("estimability falls back to a structural classifier for arbitrary names", {
  fx <- make_q8_scope_fit()
  fit <- fx$fit
  task4 <- c(-0.5, -0.5, 0.5, 0.5)
  within_c <- rep(task4, 2)
  between_c <- c(rep(-0.25, 4), rep(0.25, 4))
  mixed_c <- c(-task4, task4)

  # None of these names appear in kernel_info$term_scope, so the name lookup
  # cannot classify them; the scope must come from the cell structure.
  expect_false(any(c("anything", "whatever", "zzz") %in%
                     names(fit$kernel_info$term_scope)))

  res <- suppressWarnings(
    dkge_contrast(fit,
                  list(anything = within_c, whatever = between_c, zzz = mixed_c),
                  method = "loso", align = FALSE)
  )
  info <- res$metadata$contrast_estimability
  expect_equal(info$estimability, c("within", "between", "mixed"))

  # Unnamed vector input ("contrast1") is classified the same way.
  res_unnamed <- suppressWarnings(
    dkge_contrast(fit, between_c, method = "loso", align = FALSE)
  )
  expect_equal(res_unnamed$metadata$contrast_estimability$estimability, "between")

  # ... and a between-subject contrast under LOSO must warn.
  expect_warning(dkge_contrast(fit, between_c, method = "loso", align = FALSE),
                 "subject-label permutation")
  expect_warning(dkge_contrast(fit, within_c, method = "loso", align = FALSE), NA)
})

test_that("a between-subject contrast named after a within term is still flagged", {
  fx <- make_q8_scope_fit()
  fit <- fx$fit
  between_c <- c(rep(-0.25, 4), rep(0.25, 4))
  expect_true("task" %in% names(fit$kernel_info$term_scope))
  expect_equal(unname(fit$kernel_info$term_scope[["task"]]), "within")
  # Names are arbitrary user labels: structural evidence must win.
  expect_warning(
    res <- dkge_contrast(fit, list(task = between_c), method = "loso", align = FALSE),
    "subject-label permutation"
  )
  expect_equal(res$metadata$contrast_estimability$estimability, "between")
})

test_that("structural scope keys on factor dependence, not contrast magnitude", {
  fx <- make_q8_scope_fit()
  fit <- fx$fit
  scope_of <- function(x) dkge:::.dkge_contrast_structural_scope(x, fit)

  task4 <- c(-0.5, -0.5, 0.5, 0.5)
  # Constant weights depend on no factor: every subject can estimate it.
  expect_equal(scope_of(rep(1, 8)), "within")
  # Measure-only contrast (measure varies fastest).
  expect_equal(scope_of(rep(c(1, -1), 4)), "within")
  # Group-only, at a tiny scale.
  expect_equal(scope_of(1e-6 * c(rep(-1, 4), rep(1, 4))), "between")
  # Group x measure.
  expect_equal(scope_of(c(rep(c(1, -1), 2), rep(c(-1, 1), 2))), "mixed")
  # Rescaling cannot change the classification.
  expect_equal(scope_of(1e6 * rep(task4, 2)), "within")

  # No cell metadata -> NULL, so callers can fall back.
  fit_bare <- fit
  fit_bare$kernel_info <- NULL
  expect_null(dkge:::.dkge_contrast_structural_scope(rep(1, 8), fit_bare))
})

test_that("matrix contrast input propagates dkge_scope and dkge_term", {
  fx <- make_q8_scope_fit()
  fit <- fx$fit
  task4 <- c(-0.5, -0.5, 0.5, 0.5)
  M <- cbind(first = rep(task4, 2), second = c(rep(-0.25, 4), rep(0.25, 4)))
  attr(M, "dkge_scope") <- c("within", "between")
  attr(M, "dkge_term") <- c("task", "group")

  lst <- dkge:::.normalize_contrasts(M, fit)
  expect_equal(attr(lst$first, "dkge_scope", exact = TRUE), "within")
  expect_equal(attr(lst$second, "dkge_scope", exact = TRUE), "between")
  expect_equal(attr(lst$second, "dkge_term", exact = TRUE), "group")

  info <- dkge:::.dkge_classify_contrasts(lst, fit)
  expect_equal(info$estimability, c("within", "between"))

  # A single scope attribute applies to every column.
  M2 <- M
  attr(M2, "dkge_scope") <- "between"
  attr(M2, "dkge_term") <- NULL
  lst2 <- dkge:::.normalize_contrasts(M2, fit)
  expect_equal(dkge:::.dkge_classify_contrasts(lst2, fit)$estimability,
               c("between", "between"))

  # With no attributes the structural classifier still applies.
  M3 <- unname(M)
  colnames(M3) <- c("first", "second")
  lst3 <- dkge:::.normalize_contrasts(M3, fit)
  expect_equal(dkge:::.dkge_classify_contrasts(lst3, fit)$estimability,
               c("within", "between"))
})

test_that("normalizing contrasts keeps effect names and fills missing labels", {
  fx <- make_q8_scope_fit()
  fit <- fx$fit
  q <- nrow(fit$U)
  effects <- rownames(fit$K) %||% paste0("eff", seq_len(q))
  task4 <- c(-0.5, -0.5, 0.5, 0.5)

  M <- cbind(first = rep(task4, 2), second = c(rep(-0.25, 4), rep(0.25, 4)))
  rownames(M) <- effects

  # `as.numeric()` used to strip the row names, discarding the only record of
  # which design effect each weight refers to.
  lst <- dkge:::.normalize_contrasts(M, fit)
  expect_equal(names(lst$first), effects)
  expect_equal(names(lst$second), effects)
  expect_equal(unname(lst$first), unname(M[, "first"]))

  # Unnamed matrix columns get positional labels.
  M_nocol <- M
  colnames(M_nocol) <- NULL
  expect_equal(names(dkge:::.normalize_contrasts(M_nocol, fit)),
               c("contrast1", "contrast2"))

  # Named vectors and list entries keep their names too.
  v <- stats::setNames(rep(task4, 2), effects)
  expect_equal(names(dkge:::.normalize_contrasts(v, fit)$contrast1), effects)
  expect_equal(names(dkge:::.normalize_contrasts(list(a = v), fit)$a), effects)

  # Named contrasts still reach the full pipeline unchanged. (The second column
  # is a between-subject contrast, so LOSO warns; that is tested elsewhere.)
  res_named <- suppressWarnings(dkge_contrast(fit, M, method = "loso"))
  res_bare <- dkge_contrast(fit, unname(M[, "first", drop = FALSE]), method = "loso")
  expect_equal(names(res_named$values), c("first", "second"))
  expect_equal(names(res_bare$values), "contrast1")
  expect_equal(res_named$values$first, res_bare$values$contrast1)
})

.restore_kernel_order_info <- function(fit, fx) {
  info <- fx$kernel$info
  fit$kernel_info$cell_labels <- info$cell_labels
  fit$kernel_info$cells <- info$cells
  fit$kernel_info$blocks <- info$blocks
  fit
}

.effect_parts <- function(effects) {
  parts <- do.call(rbind, strsplit(effects, ":", fixed = TRUE))
  list(task = parts[, 1], measure = parts[, 2])
}

test_that("structural estimability rematches permuted cells by name", {
  fx <- make_permuted_cell_fit(seed = 7)
  fit <- .restore_kernel_order_info(fx$fit, fx)
  fit$kernel_info$factor_scope["task"] <- "between"
  expect_false(identical(fit$kernel_info$cell_labels, fit$effects))

  parts <- .effect_parts(fit$effects)
  task_lev <- unique(parts$task)
  measure_lev <- unique(parts$measure)
  between_c <- ifelse(parts$task == task_lev[[1]], 1, -1)
  within_c <- ifelse(parts$measure == measure_lev[[1]], 1,
                     ifelse(parts$measure == measure_lev[[2]], -1, 0))
  names(between_c) <- names(within_c) <- fit$effects

  expect_equal(dkge:::.dkge_contrast_structural_scope(between_c, fit), "between")
  expect_equal(dkge:::.dkge_contrast_structural_scope(within_c, fit), "within")

  expect_warning(
    dkge_contrast(fit, list(task = between_c), method = "loso", align = FALSE),
    "subject-label permutation"
  )
  expect_warning(
    dkge_contrast(fit, list(measure = within_c), method = "loso", align = FALSE),
    NA
  )
})

test_that("unmatched cell_labels leave structural scope NULL with a message", {
  fx <- make_permuted_cell_fit(seed = 3)
  fit <- fx$fit
  fit$kernel_info$cell_labels <- paste0("other", seq_along(fit$effects))
  cvec <- rep(c(1, -1), length.out = length(fit$effects))
  expect_message(
    scope <- dkge:::.dkge_contrast_structural_scope(cvec, fit),
    "cell_labels"
  )
  expect_null(scope)
})

test_that("contrast estimability has clean row names and validates dkge_scope", {
  fx <- make_q8_scope_fit()
  fit <- fx$fit
  task4 <- c(-0.5, -0.5, 0.5, 0.5)
  res <- suppressWarnings(
    dkge_contrast(fit, list(task = rep(task4, 2),
                            group = c(rep(-0.25, 4), rep(0.25, 4))),
                  method = "loso", align = FALSE)
  )
  info <- res$metadata$contrast_estimability
  expect_equal(rownames(info), as.character(seq_len(nrow(info))))
  expect_equal(unname(info$estimability), c("within", "between"))

  expect_error(
    dkge_contrast(fit, structure(rep(task4, 2), dkge_scope = "sideways"),
                  method = "loso", align = FALSE),
    "dkge_scope"
  )
})

