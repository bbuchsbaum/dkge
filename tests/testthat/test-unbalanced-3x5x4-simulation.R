library(testthat)

test_that("the motivating 3 by 5 by 4 grid remains a named ordinal effect space", {
  grid <- dkge_effect_grid(list(
    A = list(L = 3, levels = paste0("A", 1:3)),
    B = list(L = 5, levels = paste0("B", 1:5)),
    response = list(L = 4, type = "ordinal", levels = as.character(1:4), l = 1)
  ))
  kernel <- design_kernel(
    grid,
    terms = list("A", "B", "response", c("A", "B"),
                 c("A", "response"), c("B", "response"),
                 c("A", "B", "response")),
    basis = "cell", normalize = "max_diag"
  )
  expect_equal(nrow(grid$cells), 60L)
  expect_equal(dim(kernel$K), c(60L, 60L))
  expect_identical(rownames(kernel$K), grid$cell_labels)

  same_factor <- which(grid$cells$A == "A1" & grid$cells$B == "B1")
  response_kernel <- kernel$K[same_factor, same_factor]
  expect_gt(response_kernel[1, 2], response_kernel[1, 4])
})

test_that("partial 60-cell subjects fit without zero-padding provenance loss", {
  set.seed(73100)
  effects <- paste0("cell", seq_len(60))
  local_sets <- list(1:45, 8:55, c(1:20, 31:60))
  subjects <- lapply(seq_along(local_sets), function(s) {
    local <- local_sets[[s]]
    B <- matrix(rnorm(length(local) * 5), length(local), 5,
                dimnames = list(effects[local], NULL))
    X <- diag(length(local))
    colnames(X) <- effects[local]
    counts <- stats::setNames(4 + (local %% 7), effects[local])
    covariance <- diag(1 / counts, length(local))
    dimnames(covariance) <- list(effects[local], effects[local])
    suppressWarnings(dkge_subject(
      B, X, id = paste0("s", s), effect_n = counts,
      effect_noise_cov = covariance,
      residual_variance = rep(0.2, 5), residual_df = 10
    ))
  })
  K <- diag(60)
  dimnames(K) <- list(effects, effects)
  fit <- suppressWarnings(dkge_fit(
    dkge_data(subjects), K = K, rank = 2, w_method = "none",
    effect_weights = dkge_effect_weights("random_effects", within = "count"),
    debias = "analytic", missingness = "rescale"
  ))
  expect_setequal(fit$effects, effects)
  expect_equal(dim(fit$pair_counts), c(60L, 60L))
  expect_true(any(fit$pair_counts == 0))
  expect_true(all(is.finite(fit$Chat)))
  expect_equal(fit$effect_precision_diagnostics$method, "random_effects")
})

test_that("tracked 20-seed simulation artifacts preserve the frozen evidence", {
  replicate_path <- system.file(
    "extdata", "dkge-unbalanced-trialwise-replicates.csv", package = "dkge"
  )
  summary_path <- system.file(
    "extdata", "dkge-unbalanced-trialwise-summary.csv", package = "dkge"
  )
  expect_true(nzchar(replicate_path))
  expect_true(nzchar(summary_path))
  replicates <- utils::read.csv(replicate_path)
  summary <- utils::read.csv(summary_path)

  expect_equal(nrow(replicates), 1400L)
  expect_equal(sort(unique(replicates$seed)), 73101:73120)
  expect_equal(sort(unique(replicates$scenario)), 1:10)
  expect_setequal(unique(replicates$method), c(
    "legacy", "count", "explicit_precision", "random_effects",
    "iid_analytic", "covariance_analytic", "split_half"
  ))
  gates <- summary[summary$row_type == "gate", ]
  expect_equal(nrow(gates), 14L)
  expect_equal(sum(gates$passed), 10L)
  expect_false(gates$passed[gates$gate_id == "G5_covariance_analytic"])
  expect_equal(
    gates$observed[gates$gate_id == "G4_covariance_analytic"],
    0.05,
    tolerance = 1e-12
  )
})
