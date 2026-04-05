# test-cell-vs-cell-cross.R
# Tests that verify the cell vs cell_cross mode distinction is correctly
# implemented — different projection bases and behaviorally distinguishable
# classification results.

# ---------------------------------------------------------------------------
# Structural tests: bases and loaders
# ---------------------------------------------------------------------------

test_that("cell mode loaders use the global fit$U for every fold", {
  toy <- dkge_sim_toy(
    factors = list(A = list(L = 2)), active_terms = "A",
    S = 5, P = 12, snr = 5, seed = 1
  )
  fit  <- dkge(toy$B_list, toy$X_list, K = toy$K, rank = 2)
  fold_bundle <- .dkge_prepare_folds(fit, NULL)

  fi <- .dkge_build_global_fold_loaders(fit, fold_bundle$assignments)

  for (f in fi$folds) {
    expect_equal(f$basis, fit$U,
                 info = sprintf("fold %d basis should equal fit$U", f$index))
  }
})

test_that("cell_cross mode loaders use fold-specific LOSO bases", {
  toy <- dkge_sim_toy(
    factors = list(A = list(L = 2)), active_terms = "A",
    S = 5, P = 12, snr = 5, seed = 1
  )
  fit  <- dkge(toy$B_list, toy$X_list, K = toy$K, rank = 2)
  fold_bundle <- .dkge_prepare_folds(fit, NULL)

  fi <- .dkge_build_fold_bases(fit, fold_bundle$assignments,
                                align = FALSE, loader_scope = "all")

  # At least one fold basis must differ from the global fit$U
  any_different <- any(vapply(fi$folds, function(f) {
    !isTRUE(all.equal(f$basis, fit$U, tolerance = 1e-8))
  }, logical(1)))

  expect_true(any_different)
})

test_that("held-out subject Y matrix differs between cell and cell_cross", {
  toy <- dkge_sim_toy(
    factors = list(A = list(L = 2)), active_terms = "A",
    S = 5, P = 12, snr = 5, seed = 3
  )
  fit <- dkge(toy$B_list, toy$X_list, K = toy$K, rank = 2)
  fold_bundle <- .dkge_prepare_folds(fit, NULL)

  fi_cell  <- .dkge_build_global_fold_loaders(fit, fold_bundle$assignments)
  fi_cross <- .dkge_build_fold_bases(fit, fold_bundle$assignments,
                                     align = FALSE, loader_scope = "all")

  # For each fold, compare the held-out subject's Y matrix
  n_differ <- 0L
  for (fold_idx in seq_along(fi_cell$folds)) {
    s <- fi_cell$folds[[fold_idx]]$subjects[[1]]
    key <- as.character(s)
    Y_cell  <- fi_cell$folds[[fold_idx]]$loaders[[key]]$Y
    Y_cross <- fi_cross$folds[[fold_idx]]$loaders[[key]]$Y
    if (!isTRUE(all.equal(Y_cell, Y_cross, tolerance = 1e-8))) {
      n_differ <- n_differ + 1L
    }
  }
  expect_gt(n_differ, 0L)
})

# ---------------------------------------------------------------------------
# Behavioral tests: do the two modes give different classification results?
# ---------------------------------------------------------------------------

test_that("cell and cell_cross produce different row-level feature values", {
  # Use a small dataset where leaving one subject out changes the basis
  toy <- dkge_sim_toy(
    factors = list(A = list(L = 3)), active_terms = "A",
    S = 6, P = 15, snr = 4, seed = 99
  )
  kern <- design_kernel(factors = list(A = list(L = 3)), basis = "effect")
  fit  <- dkge(toy$B_list, toy$X_list, K = kern, rank = 2)

  clf_cell  <- dkge_classify(fit, targets = ~A, method = "lda", mode = "cell",
                              verbose = FALSE)
  clf_cross <- dkge_classify(fit, targets = ~A, method = "lda", mode = "cell_cross",
                              verbose = FALSE)

  res_cell  <- clf_cell$results[[1]]
  res_cross <- clf_cross$results[[1]]

  # Modes should be recorded correctly
  expect_equal(res_cell$mode,  "cell")
  expect_equal(res_cross$mode, "cell_cross")

  # Feature matrices differ → predictions or probabilities differ for at least
  # one subject (this is nearly certain when bases differ)
  probs_identical <- isTRUE(all.equal(res_cell$probabilities,
                                      res_cross$probabilities,
                                      tolerance = 1e-6))
  expect_false(probs_identical)
})

# ---------------------------------------------------------------------------
# Leakage documentation test (cell mode by design)
# ---------------------------------------------------------------------------

test_that("cell mode basis includes held-out subject (documented leakage)", {
  # cell mode is NOT a fully cross-validated pipeline for the basis step:
  # fit$U is estimated on ALL subjects, then the classifier is cross-validated.
  # This test documents that property explicitly.
  toy <- dkge_sim_toy(
    factors = list(A = list(L = 2)), active_terms = "A",
    S = 4, P = 10, snr = 5, seed = 5
  )
  fit <- dkge(toy$B_list, toy$X_list, K = toy$K, rank = 1)
  fold_bundle <- .dkge_prepare_folds(fit, NULL)

  fi_cell <- .dkge_build_global_fold_loaders(fit, fold_bundle$assignments)

  # The basis for every fold is exactly fit$U, which was estimated on all S subjects
  for (f in fi_cell$folds) {
    expect_equal(f$basis, fit$U)
  }
  # cell_cross avoids this by recomputing the basis without the held-out subject
  fi_cross <- .dkge_build_fold_bases(fit, fold_bundle$assignments,
                                     align = FALSE, loader_scope = "all")
  # At least one fold re-estimated its basis
  refit <- any(vapply(fi_cross$folds, function(f) {
    !isTRUE(all.equal(f$basis, fit$U, tolerance = 1e-8))
  }, logical(1)))
  expect_true(refit)
})

# ---------------------------------------------------------------------------
# auto mode selects cell for within-subject targets
# ---------------------------------------------------------------------------

test_that("auto mode selects cell for within-subject targets", {
  toy <- dkge_sim_toy(
    factors = list(A = list(L = 2)), active_terms = "A",
    S = 4, P = 10, snr = 5, seed = 11
  )
  kern <- design_kernel(factors = list(A = list(L = 2)), basis = "effect")
  fit  <- dkge(toy$B_list, toy$X_list, K = kern, rank = 1)

  clf <- dkge_classify(fit, targets = ~A, method = "lda", mode = "auto",
                       verbose = FALSE)
  # within_subject scope → auto should pick "cell"
  expect_equal(clf$results[[1]]$mode, "cell")
})
