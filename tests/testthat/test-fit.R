# test-fit.R
# Direct diagnostics for dkge_fit helpers

library(testthat)

make_fit_fixture <- function(S = 3, q = 3, P = 4, T = 20, seed = 900) {
  set.seed(seed)
  betas <- replicate(S, matrix(rnorm(q * P), q, P), simplify = FALSE)
  designs <- replicate(S, {
    X <- matrix(rnorm(T * q), T, q)
    qr.Q(qr(X))
  }, simplify = FALSE)
  list(betas = betas, designs = designs, K = diag(q))
}

test_that("pooled design Cholesky matches Gram matrix", {
  fixture <- make_fit_fixture()
  ruler <- dkge:::.dkge_compute_shared_ruler(fixture$designs)
  expect_equal(ruler$G_pool, Reduce(`+`, lapply(fixture$designs, crossprod)))
  expect_equal(t(ruler$R) %*% ruler$R, ruler$G_pool, tolerance = 1e-10)
})

test_that("row standardisation matches manual computation", {
  fixture <- make_fit_fixture()
  ruler <- dkge:::.dkge_compute_shared_ruler(fixture$designs)
  Btil <- dkge:::.dkge_row_standardize(fixture$betas, ruler$R)
  expect_equal(Btil[[1]], t(ruler$R) %*% fixture$betas[[1]])
})

test_that("kernel roots reconstruct original kernel", {
  K <- matrix(c(2, 0.5, 0.5, 1.5), 2, 2)
  roots <- dkge:::.dkge_kernel_roots(K)
  expect_equal(roots$Khalf %*% roots$Khalf, K, tolerance = 1e-10)
  expect_equal(roots$Khalf %*% roots$Kihalf, diag(2), tolerance = 1e-10)
})

test_that("subject weights handle MFA and Omega matrices", {
  Btil <- list(matrix(1:6, 3, 2), matrix(2:7, 3, 2))
  Omega_list <- list(diag(c(1, 2)), matrix(c(2, 1, 1, 3), 2))
  Khalf <- diag(3)

  w_mfa <- dkge:::.dkge_subject_weights(Btil, list(NULL, NULL), Khalf,
                                        w_method = "mfa_sigma1", w_tau = 0)
  expect_true(all(is.finite(w_mfa)))

  w_mat <- dkge:::.dkge_subject_weights(Btil, Omega_list, Khalf,
                                        w_method = "energy", w_tau = 0.5)
  expect_equal(length(w_mat), length(Btil))
  expect_true(all(w_mat > 0))
})

test_that("subject weights apply the full Khalf congruence to zero-filled rows", {
  K <- matrix(0.35, 3, 3)
  diag(K) <- 1
  Khalf <- dkge:::.dkge_kernel_roots(K)$Khalf

  # Rows a subject does not observe are already zero in raw effect space, so
  # masking must be a no-op and the FULL Khalf congruence is what the
  # eigensolve sees (dkge-moments.R fit-path convention).
  B3 <- list(
    matrix(c(1, 2,
             3, 4,
             0, 0), 3, 2, byrow = TRUE),
    matrix(c(5, 1,
             0, 0,
             2, 7), 3, 2, byrow = TRUE)
  )
  masks <- list(c(TRUE, TRUE, FALSE), c(TRUE, FALSE, TRUE))

  w_masked <- dkge:::.dkge_subject_weights(B3, list(NULL, NULL), Khalf,
                                           w_method = "energy", w_tau = 0,
                                           obs_masks = masks)
  w_unmasked <- dkge:::.dkge_subject_weights(B3, list(NULL, NULL), Khalf,
                                             w_method = "energy", w_tau = 0)

  # Passing masks changes nothing: the zero rows carry no energy already.
  expect_equal(w_masked, w_unmasked, tolerance = 1e-12)

  # Golden values, computed once from K = 0.35 off-diagonal / 1 diagonal.
  expect_equal(w_masked, c(1.41368585, 0.58631415), tolerance = 1e-7)

  # The superseded convention (Khalf[obs, obs] on the observed rows only)
  # gives a measurably different answer, so this pins the choice.
  expect_false(isTRUE(all.equal(w_masked, c(1.41529552, 0.58470448),
                                tolerance = 1e-7)))

  # The energy is exactly the trace of the subject's transformed effect moment
  # (R = I here), i.e. the fit-path contribution to Chat.
  energy <- vapply(B3, function(B) {
    M <- dkge:::.dkge_effect_moment(B)
    sum(diag(dkge:::.dkge_transform_effect_moment(M, diag(3), Khalf)))
  }, numeric(1))
  raw <- 1 / (energy + 1e-12)
  expect_equal(w_masked, raw / mean(raw), tolerance = 1e-10)
})

test_that("a subject observing no effect rows gets weight zero", {
  Khalf <- diag(2)
  Btil <- list(matrix(c(1, 2, 3, 4), 2, 2), matrix(0, 2, 2))
  masks <- list(c(TRUE, TRUE), c(FALSE, FALSE))

  expect_warning(
    w <- dkge:::.dkge_subject_weights(Btil, list(NULL, NULL), Khalf,
                                      w_method = "energy", w_tau = 0,
                                      obs_masks = masks),
    "no effect rows"
  )
  expect_equal(w, c(1, 0))

  # Shrinkage toward equal weights must not resurrect the excluded subject.
  expect_warning(
    w_tau <- dkge:::.dkge_subject_weights(Btil, list(NULL, NULL), Khalf,
                                          w_method = "energy", w_tau = 0.5,
                                          obs_masks = masks),
    "no effect rows"
  )
  expect_equal(w_tau, c(1, 0))
})

test_that("missingness helper implements shared coverage policies", {
  Chat <- matrix(c(4, 2, 6,
                   2, 8, 3,
                   6, 3, 9), 3, 3, byrow = TRUE)
  pair_counts <- matrix(c(2L, 1L, 0L,
                          1L, 2L, 1L,
                          0L, 1L, 1L), 3, 3, byrow = TRUE)

  expect_equal(
    dkge:::.dkge_apply_missingness(Chat, pair_counts, "none"),
    Chat
  )

  # rescale: hand-computed entrywise Chat / pair_counts, 0 where uncovered.
  expected_rescale <- matrix(c(2, 2, 0,
                               2, 4, 3,
                               0, 3, 9), 3, 3, byrow = TRUE)
  expect_equal(
    dkge:::.dkge_apply_missingness(Chat, pair_counts, "rescale"),
    expected_rescale
  )

  expected_mask <- matrix(c(4, 0, 0,
                            0, 8, 0,
                            0, 0, 0), 3, 3, byrow = TRUE)
  expect_equal(
    dkge:::.dkge_apply_missingness(Chat, pair_counts, "mask",
                                   list(min_pairs = 2L)),
    expected_mask
  )

  # shrink (gamma = 1): rel = pair_counts / 2, blended toward the diagonal of
  # the RESCALED matrix (2, 4, 9), not of the raw Chat.
  #   [1,1] rel 1   -> 2
  #   [1,2] rel 0.5 -> 0.5 * 2 + 0.5 * 0     = 1
  #   [1,3] rel 0   -> 0 * 0    + 1 * 0      = 0
  #   [2,2] rel 1   -> 4
  #   [2,3] rel 0.5 -> 0.5 * 3 + 0.5 * 0     = 1.5
  #   [3,3] rel 0.5 -> 0.5 * 9 + 0.5 * 9     = 9
  expected_shrink <- matrix(c(2, 1, 0,
                              1, 4, 1.5,
                              0, 1.5, 9), 3, 3, byrow = TRUE)
  expect_equal(
    dkge:::.dkge_apply_missingness(Chat, pair_counts, "shrink",
                                   list(gamma = 1)),
    expected_shrink
  )
})

test_that("missingness arguments are validated the same way on both paths", {
  Chat <- diag(2)
  pair_counts <- matrix(1L, 2, 2)
  for (bad in list(-1, c(1, 2), NA_real_, Inf, "2")) {
    expect_error(
      dkge:::.dkge_apply_missingness(Chat, pair_counts, "mask",
                                     list(min_pairs = bad)),
      "min_pairs"
    )
    expect_error(
      dkge:::.dkge_apply_missingness(Chat, pair_counts, "shrink",
                                     list(gamma = bad)),
      "gamma"
    )
  }
})

test_that("missingness policies read pair_counts as a fractional weighted mass", {
  # `pair_counts` is sum_s a_s * 1[pair observed]; the multiplier bootstrap
  # schemes ("exp", "bayes") make the sample weights a_s fractional, so the
  # policies must gate on positive mass and put `min_pairs` back on the
  # subject-count scale. `pmax(pc, 1L)` turned rescale/shrink into no-ops and
  # `pc < min_pairs` zeroed everything.
  effects <- c("e1", "e2", "e3")
  design <- diag(3)
  colnames(design) <- effects
  set.seed(404)
  rows <- list(c(1L, 2L), c(2L, 3L), 1:3)
  subjects <- lapply(seq_along(rows), function(i) {
    B <- matrix(rnorm(9), 3, 3, dimnames = list(effects, NULL))
    B[-rows[[i]], ] <- 0
    suppressWarnings(dkge_subject(B, design, id = paste0("s", i),
                                  observed_rows = rows[[i]]))
  })
  data <- dkge_data(subjects)
  structural <- lapply(data$provenance$obs_mask,
                       function(m) tcrossprod(as.numeric(m)))
  moments <- lapply(data$betas, tcrossprod)
  xi <- c(0.2, 0.3, 0.25)   # fractional and unequal, as a multiplier draw is

  for (miss in c("rescale", "mask", "shrink")) {
    fit <- dkge_fit(data, K = diag(3), rank = 2, w_method = "none",
                    effect_scaling = "none", missingness = miss,
                    miss_args = list(min_pairs = 2, gamma = 1))

    # Unit sample weights must reproduce the fit's own Chat exactly.
    unit <- dkge:::.dkge_repool_fit(fit, sample_weights = rep(1, 3))
    expect_equal(unname(unit$Chat), unname(fit$Chat), tolerance = 1e-12,
                 info = miss)

    boot <- dkge:::.dkge_repool_fit(fit, sample_weights = xi)
    raw <- dkge:::.dkge_repool_fit(fit, sample_weights = xi,
                                   missingness = "none")
    expect_false(isTRUE(all.equal(boot$Chat, raw$Chat)), info = miss)

    num <- Reduce(`+`, Map(function(w, M) w * M, xi, moments))
    den <- Reduce(`+`, Map(function(w, S) w * S, xi, structural))
    expect_equal(unname(boot$pair_counts), unname(den), tolerance = 1e-12,
                 info = miss)
    per_pair <- matrix(0, 3, 3)
    per_pair[den > 0] <- num[den > 0] / den[den > 0]

    expected <- switch(
      miss,
      rescale = per_pair,
      mask = {
        # den / mean(xi) is the coverage in subjects: 3 for (e2,e2), 2.2 for
        # the pairs s2/s3 share, 1.8 or 1 for the rest.
        keep <- (den / mean(xi)) >= 2
        M <- num
        M[!keep] <- 0
        M
      },
      shrink = {
        rel <- den / max(den)
        rel * per_pair + (1 - rel) * diag(diag(per_pair), 3, 3)
      }
    )
    expected <- (expected + t(expected)) / 2
    expect_equal(unname(boot$pooled), unname(expected), tolerance = 1e-12,
                 info = miss)
  }
})

test_that("training blocks and dkge_transform_block() agree in dimnames", {
  set.seed(7)
  S <- 2; q <- 3; P <- 4; T <- 10
  effects <- c("a", "b", "c")
  betas <- replicate(S, matrix(rnorm(q * P), q, P,
                               dimnames = list(effects, paste0("v", 1:P))),
                     simplify = FALSE)
  designs <- replicate(S, {
    X <- qr.Q(qr(matrix(rnorm(T * q), T, q)))
    colnames(X) <- effects
    X
  }, simplify = FALSE)
  K <- diag(q)
  dimnames(K) <- list(effects, effects)

  fit <- dkge_fit(betas, designs, K = K, rank = 2, keep_X = TRUE,
                  w_method = "none")
  block <- fit$X_concat[, fit$block_indices[[1]], drop = FALSE]
  transformed <- dkge_transform_block(fit, betas[[1]], w_s = fit$weights[1])

  expect_equal(unname(block), unname(transformed), tolerance = 1e-12)
  expect_identical(dimnames(block), dimnames(transformed))
  # Rows are K^{1/2}-whitened mixtures, so they carry no effect labels; parcel
  # labels are genuine and are preserved.
  expect_null(rownames(block))
  expect_equal(colnames(block), colnames(betas[[1]]))
})

test_that("kernel row/column labels must agree before any other reconciliation", {
  K <- diag(3)
  dimnames(K) <- list(c("a", "b", "c"), c("c", "b", "a"))

  # The check has to run before the placeholder and duplicate-label shortcuts,
  # or a kernel whose rows and columns disagree slips through unnoticed.
  expect_error(
    dkge:::.dkge_align_kernel_effects(K, dkge:::.default_effect_names(3)),
    "row and column names must be identical"
  )
  expect_error(
    dkge:::.dkge_align_kernel_effects(K, c("a", "b", "c")),
    "row and column names must be identical"
  )

  K_dup <- diag(3)
  dimnames(K_dup) <- list(c("a", "a", "b"), c("a", "b", "b"))
  expect_error(
    dkge:::.dkge_align_kernel_effects(K_dup, c("a", "a", "b")),
    "row and column names must be identical"
  )

  fixture <- make_fit_fixture(q = 3)
  expect_error(
    dkge_fit(fixture$betas, fixture$designs, K = K, rank = 2),
    "row and column names must be identical"
  )
})

test_that("a labelled kernel names the effects when the data only has placeholders", {
  fixture <- make_fit_fixture(q = 3)
  labels <- c("cue", "load", "cue:load")
  K <- diag(3)
  dimnames(K) <- list(labels, labels)

  fit <- dkge_fit(fixture$betas, fixture$designs, K = K, rank = 2,
                  w_method = "none")
  expect_equal(fit$effects, labels)
  expect_equal(rownames(fit$K), labels)
  expect_equal(colnames(fit$K), labels)

  # Duplicated kernel labels cannot identify effects; they are dropped rather
  # than left to contradict `fit$effects`.
  K_dup <- diag(3)
  dimnames(K_dup) <- list(c("A1", "A1", "A2"), c("A1", "A1", "A2"))
  fit_dup <- dkge_fit(fixture$betas, fixture$designs, K = K_dup, rank = 2,
                      w_method = "none")
  expect_equal(fit_dup$effects, dkge:::.default_effect_names(3))
  expect_null(dimnames(fit_dup$K))
})

test_that("kernel_info is permuted with the kernel to match data effects", {
  fx <- make_permuted_cell_fit(seed = 7)
  fit <- fx$fit
  expect_identical(fit$effects, fx$permuted_labels)
  expect_identical(fit$kernel_info$cell_labels, fit$effects)
  expect_identical(rownames(fit$K), fit$effects)

  cells <- as.data.frame(fit$kernel_info$cells, stringsAsFactors = FALSE)
  split_eff <- do.call(rbind, strsplit(fit$effects, ":", fixed = TRUE))
  expect_equal(as.character(cells$task), split_eff[, 1])
  expect_equal(as.character(cells$measure), split_eff[, 2])
  expect_identical(fit$kernel_info$blocks$cells, match(seq_len(6), fx$perm))

  expect_silent(dkge:::.dkge_check_kernel_info(fit))
  bad <- fit
  bad$kernel_info$cell_labels <- rev(bad$kernel_info$cell_labels)
  expect_error(dkge:::.dkge_check_kernel_info(bad), "cell_labels")

  # Direct aligner: cells, labels, and blocks follow idx.
  labels <- c("a", "b", "c")
  K <- diag(3)
  dimnames(K) <- list(labels, labels)
  info <- list(
    cell_labels = labels,
    cells = data.frame(f = labels, stringsAsFactors = FALSE),
    blocks = list(cells = 1:3)
  )
  aligned <- suppressMessages(
    dkge:::.dkge_align_kernel_effects(K, c("c", "a", "b"), info)
  )
  expect_equal(aligned$kernel_info$cell_labels, c("c", "a", "b"))
  expect_equal(as.character(aligned$kernel_info$cells$f), c("c", "a", "b"))
  expect_equal(aligned$kernel_info$blocks$cells, match(1:3, c(3L, 1L, 2L)))

  # Consumer smoke: once info matches fit$effects, the default design basis
  # is the sum-contrast model matrix of those cells (D1 runs the rest).
  C <- dkge_design_basis(fit, normalize = "none")
  cells_f <- cells
  for (nm in c("task", "measure")) {
    cells_f[[nm]] <- factor(cells_f[[nm]],
                            levels = unique(as.character(cells_f[[nm]])))
  }
  mm <- stats::model.matrix(
    ~ task * measure, data = cells_f,
    contrasts.arg = list(task = stats::contr.sum(nlevels(cells_f$task)),
                         measure = stats::contr.sum(nlevels(cells_f$measure)))
  )
  # Drop intercept name differences; compare the task main-effect column
  # up to a positive scale (unit_K is applied only under the default).
  task_mm <- mm[, attr(mm, "assign") == 1L]
  task_C <- C[, attr(C, "term") == "task"]
  ratio <- as.numeric(task_C / task_mm)
  expect_equal(ratio, rep(ratio[[1]], length(ratio)), tolerance = 1e-10)
})

test_that("min_pairs is compared on the unweighted subject scale", {
  Chat <- matrix(c(4, 2, 2, 8), 2, 2)
  # Half a unit of mass per subject: two subjects share both pairs.
  pair_counts <- matrix(c(1, 1, 1, 1), 2, 2) * 0.5
  kept <- dkge:::.dkge_apply_missingness(Chat, pair_counts, "mask",
                                         list(min_pairs = 2),
                                         weight_scale = 0.25)
  expect_equal(kept, (Chat + t(Chat)) / 2)
  dropped <- dkge:::.dkge_apply_missingness(Chat, pair_counts, "mask",
                                            list(min_pairs = 3),
                                            weight_scale = 0.25)
  expect_equal(dropped, matrix(0, 2, 2))

  # rescale divides by the fractional mass itself, not by pmax(mass, 1).
  expect_equal(
    dkge:::.dkge_apply_missingness(Chat, pair_counts, "rescale"),
    (Chat + t(Chat)) / 2 / 0.5
  )
})

test_that("dkge_fit accepts raw lists and dkge_data", {
  fixture <- make_fit_fixture()
  fit1 <- dkge_fit(fixture$betas, fixture$designs, K = fixture$K, rank = 2,
                   keep_X = TRUE)
  data_bundle <- dkge_data(fixture$betas, designs = fixture$designs)
  fit2 <- dkge_fit(data_bundle, K = fixture$K, rank = 2, keep_X = FALSE)

  expect_equal(fit1$rank, 2)
  expect_equal(fit2$rank, 2)
  expect_true(!is.null(fit1$X_concat))
  expect_null(fit2$X_concat)
  expect_equal(fit1$K, fixture$K)
})

test_that("dkge_fit can leave effect rows on input scale", {
  fixture <- make_fit_fixture(S = 3, q = 4, P = 6, T = 20)
  fit <- dkge_fit(
    fixture$betas,
    fixture$designs,
    K = fixture$K,
    rank = 2,
    w_method = "none",
    effect_scaling = "none"
  )

  expect_equal(fit$effect_scaling, "none")
  expect_equal(fit$R, diag(4), tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(fit$Btil, fixture$betas, tolerance = 1e-12, ignore_attr = TRUE)
})

test_that("dkge_fit default missingness leaves full-coverage fits unchanged", {
  fixture <- make_fit_fixture()
  fit_default <- dkge_fit(fixture$betas, fixture$designs, K = fixture$K,
                          rank = 2, w_method = "none")
  fit_none <- dkge_fit(fixture$betas, fixture$designs, K = fixture$K,
                       rank = 2, w_method = "none", missingness = "none")

  expect_equal(fit_default$Chat, fit_none$Chat)
  expect_equal(fit_default$pair_counts, fit_none$pair_counts)
  expect_equal(unique(as.integer(fit_default$pair_counts)), length(fixture$betas))
})

test_that("dkge_fit applies partial-effect missingness to training Chat", {
  betas <- list(
    matrix(c(1, 2,
             3, 4), 2, 2, byrow = TRUE,
           dimnames = list(c("e1", "e2"), NULL)),
    matrix(c(5, 6,
             7, 8), 2, 2, byrow = TRUE,
           dimnames = list(c("e2", "e3"), NULL))
  )
  designs <- list(
    matrix(c(1, 0,
             0, 1,
             1, 0,
             0, 1), 4, 2, byrow = TRUE,
           dimnames = list(NULL, c("e1", "e2"))),
    matrix(c(1, 0,
             0, 1,
             1, 0,
             0, 1), 4, 2, byrow = TRUE,
           dimnames = list(NULL, c("e2", "e3")))
  )
  data <- dkge_data(betas, designs, subject_ids = c("s1", "s2"))
  K <- matrix(0.25, 3, 3)
  diag(K) <- 1
  dimnames(K) <- list(data$effects, data$effects)

  fit_none <- dkge_fit(data, K = K, rank = 2, w_method = "none",
                       missingness = "none")
  fit_mask <- dkge_fit(data, K = K, rank = 2, w_method = "none",
                       missingness = "mask",
                       miss_args = list(min_pairs = 1L))

  expect_equal(fit_mask$pair_counts["e1", "e3"], 0L)
  # Coverage policies act in raw effect space, before R and K mix rows.
  expect_equal(unname(fit_mask$effect_moment["e1", "e3"]), 0)
  expect_equal(unname(fit_mask$effect_moment["e3", "e1"]), 0)

  ctx <- dkge:::.dkge_fold_weight_context(
    fit_none,
    train_ids = seq_along(fit_none$Btil),
    missingness = "mask",
    miss_args = list(min_pairs = 1L)
  )
  expect_equal(fit_mask$Chat, ctx$Chat, tolerance = 1e-10)
})

test_that("dkge_fit defaults to full rank when rank is NULL", {
  fixture <- make_fit_fixture(S = 4, q = 12, P = 20, T = 35)
  fit <- dkge_fit(fixture$betas, fixture$designs, K = fixture$K, rank = NULL,
                  w_method = "none")

  expect_equal(fit$rank_requested, 12)
  expect_equal(fit$rank, 12)
})

test_that("dkge_fit honours ridge and Omega weighting", {
  fixture <- make_fit_fixture()
  Omega <- lapply(fixture$betas, function(B) runif(ncol(B)))
  fit <- dkge_fit(fixture$betas, fixture$designs, K = fixture$K, rank = 2,
                  ridge = 0.5, Omega_list = Omega)
  expect_true(all(eigen(fit$Chat, symmetric = TRUE)$values >= -1e-8))
  expect_equal(length(fit$weights), length(fixture$betas))
})

test_that("dkge_fit validates kernel PSD and finite values", {
  fixture <- make_fit_fixture(q = 3)

  K_bad_psd <- diag(3)
  K_bad_psd[1, 1] <- -1
  expect_error(
    dkge_fit(fixture$betas, fixture$designs, K = K_bad_psd, rank = 2),
    "positive semidefinite"
  )

  K_bad_finite <- fixture$K
  K_bad_finite[1, 2] <- Inf
  K_bad_finite[2, 1] <- Inf
  expect_error(
    dkge_fit(fixture$betas, fixture$designs, K = K_bad_finite, rank = 2),
    "non-finite"
  )

  K_bad_sym <- fixture$K
  K_bad_sym[1, 2] <- 0.9
  K_bad_sym[2, 1] <- 0
  expect_error(
    dkge_fit(fixture$betas, fixture$designs, K = K_bad_sym, rank = 2),
    "symmetric"
  )
})

test_that("cpca arguments are respected", {
  fixture <- make_fit_fixture(q = 5, P = 8)
  fit <- dkge_fit(fixture$betas, fixture$designs, K = fixture$K, rank = 3,
                  cpca_blocks = 1:2, cpca_part = "design", cpca_ridge = 0.1)
  expect_true(!is.null(fit$cpca))
  expect_equal(fit$cpca$part, "design")
  expect_equal(fit$cpca$blocks, 1:2)
})

test_that("dkge_fit detects dimension mismatches", {
  fixture <- make_fit_fixture()
  bad_designs <- fixture$designs
  bad_designs[[1]] <- bad_designs[[1]][, -1]
  expect_error(dkge_fit(fixture$betas, bad_designs, K = fixture$K, rank = 2))
})

# -------------------------------------------------------------------------
# K-orthonormality property tests
# -------------------------------------------------------------------------

test_that("K-orthonormality holds with identity kernel", {
  withr::local_seed(1001)
  fixture <- make_fit_fixture(S = 3, q = 4, P = 5, T = 20)
  fit <- dkge_fit(fixture$betas, fixture$designs, K = fixture$K, rank = 3,
                  w_method = "none")

  # U^T K U = I_r (identity kernel K=I so this is simply orthonormality)
  UtKU <- t(fit$U) %*% fit$K %*% fit$U
  expect_lt(max(abs(UtKU - diag(fit$rank))), 1e-8)
})

test_that("K-orthonormality holds with RBF (ordinal) kernel", {
  withr::local_seed(1002)
  q <- 5
  # Construct an RBF kernel for ordinal data
  dists <- as.matrix(dist(1:q))
  sigma <- 1.5
  K <- exp(-dists^2 / (2 * sigma^2))
  K <- (K + t(K)) / 2  # ensure symmetric

  fixture <- make_fit_fixture(S = 4, q = q, P = 6, T = 25)
  fit <- dkge_fit(fixture$betas, fixture$designs, K = K, rank = 3,
                  w_method = "none")

  UtKU <- t(fit$U) %*% fit$K %*% fit$U
  expect_lt(max(abs(UtKU - diag(fit$rank))), 1e-8)
})

test_that("K-orthonormality holds with multi-factor kernel", {
  withr::local_seed(1003)
  q <- 6
  # Multi-factor kernel: block diagonal with 2 factors
  K <- matrix(0, q, q)
  K[1:3, 1:3] <- 1
  K[4:6, 4:6] <- 1
  diag(K) <- 1
  K <- K + 0.1 * diag(q)  # ensure PSD

  fixture <- make_fit_fixture(S = 3, q = q, P = 9, T = 20)
  fit <- dkge_fit(fixture$betas, fixture$designs, K = K, rank = 4,
                  w_method = "none")

  UtKU <- t(fit$U) %*% fit$K %*% fit$U
  expect_lt(max(abs(UtKU - diag(fit$rank))), 1e-8)
})

test_that("K-orthonormality holds for rank=1 edge case", {
  withr::local_seed(1004)
  fixture <- make_fit_fixture(S = 3, q = 4, P = 5, T = 20)
  fit <- dkge_fit(fixture$betas, fixture$designs, K = fixture$K, rank = 1,
                  w_method = "none")

  expect_equal(fit$rank, 1)
  UtKU <- t(fit$U) %*% fit$K %*% fit$U
  expect_lt(max(abs(UtKU - diag(fit$rank))), 1e-8)
})

test_that("K-orthonormality holds for full rank (rank=q)", {
  withr::local_seed(1005)
  q <- 4
  fixture <- make_fit_fixture(S = 3, q = q, P = 5, T = 20)
  fit <- dkge_fit(fixture$betas, fixture$designs, K = fixture$K, rank = q,
                  w_method = "none")

  # May be less than q if some eigenvalues are near zero
  expect_lte(fit$rank, q)
  UtKU <- t(fit$U) %*% fit$K %*% fit$U
  expect_lt(max(abs(UtKU - diag(fit$rank))), 1e-8)
})

test_that("K-orthonormality holds for near-full rank (rank=q-1)", {
  withr::local_seed(1006)
  q <- 5
  fixture <- make_fit_fixture(S = 3, q = q, P = 6, T = 20)
  fit <- dkge_fit(fixture$betas, fixture$designs, K = fixture$K, rank = q - 1,
                  w_method = "none")

  expect_lte(fit$rank, q - 1)
  UtKU <- t(fit$U) %*% fit$K %*% fit$U
  expect_lt(max(abs(UtKU - diag(fit$rank))), 1e-8)
})

test_that("K-orthonormality holds with ridge > 0", {
  withr::local_seed(1007)
  fixture <- make_fit_fixture(S = 3, q = 4, P = 5, T = 20)
  fit <- dkge_fit(fixture$betas, fixture$designs, K = fixture$K, rank = 3,
                  ridge = 0.5, w_method = "none")

  UtKU <- t(fit$U) %*% fit$K %*% fit$U
  expect_lt(max(abs(UtKU - diag(fit$rank))), 1e-8)
})

test_that("K-orthonormality holds with MFA subject weights", {
  withr::local_seed(1008)
  fixture <- make_fit_fixture(S = 4, q = 4, P = 5, T = 20)
  fit <- dkge_fit(fixture$betas, fixture$designs, K = fixture$K, rank = 3,
                  w_method = "mfa_sigma1", w_tau = 0.3)

  UtKU <- t(fit$U) %*% fit$K %*% fit$U
  expect_lt(max(abs(UtKU - diag(fit$rank))), 1e-8)
})

test_that("K-orthonormality holds with energy subject weights", {
  withr::local_seed(1009)
  fixture <- make_fit_fixture(S = 4, q = 4, P = 5, T = 20)
  fit <- dkge_fit(fixture$betas, fixture$designs, K = fixture$K, rank = 3,
                  w_method = "energy", w_tau = 0.5)

  UtKU <- t(fit$U) %*% fit$K %*% fit$U
  expect_lt(max(abs(UtKU - diag(fit$rank))), 1e-8)
})

# -------------------------------------------------------------------------
# Pooled design verification tests
# -------------------------------------------------------------------------

test_that("pooled Gram with heterogeneous T_s values", {
  withr::local_seed(2001)
  q <- 4
  S <- 5
  # Different T_s for each subject (heterogeneous time points)
  T_values <- c(15, 20, 25, 18, 22)

  betas <- replicate(S, matrix(rnorm(q * 6), q, 6), simplify = FALSE)
  designs <- lapply(T_values, function(T_s) {
    X <- matrix(rnorm(T_s * q), T_s, q)
    qr.Q(qr(X))
  })

  ruler <- dkge:::.dkge_compute_shared_ruler(designs)

  # G_pool = sum_s X_s^T X_s
  G_manual <- Reduce(`+`, lapply(designs, crossprod))
  expect_equal(ruler$G_pool, G_manual, tolerance = 1e-10)

  # R^T R = G_pool (Cholesky decomposition property)
  expect_equal(t(ruler$R) %*% ruler$R, ruler$G_pool, tolerance = 1e-10)
})

test_that("pooled Gram accumulation with many subjects", {

  withr::local_seed(2002)
  q <- 5
  S <- 10  # Many subjects to verify accumulation doesn't drift
  T_s <- 20
  P <- 8

  betas <- replicate(S, matrix(rnorm(q * P), q, P), simplify = FALSE)
  designs <- replicate(S, {
    X <- matrix(rnorm(T_s * q), T_s, q)
    qr.Q(qr(X))
  }, simplify = FALSE)

  ruler <- dkge:::.dkge_compute_shared_ruler(designs)

  # Manual accumulation of G_pool
  G_manual <- matrix(0, q, q)
  for (s in seq_len(S)) {
    G_manual <- G_manual + crossprod(designs[[s]])
  }

  expect_equal(ruler$G_pool, G_manual, tolerance = 1e-10)
  expect_equal(t(ruler$R) %*% ruler$R, ruler$G_pool, tolerance = 1e-10)

  # Verify PSD property of G_pool
  eigs <- eigen(ruler$G_pool, symmetric = TRUE)$values
  expect_true(all(eigs > 0))
})

test_that("row standardization produces expected scaling", {
  withr::local_seed(2003)
  q <- 4
  S <- 3
  T_s <- 18
  P <- 5

  betas <- replicate(S, matrix(rnorm(q * P), q, P), simplify = FALSE)
  designs <- replicate(S, {
    X <- matrix(rnorm(T_s * q), T_s, q)
    qr.Q(qr(X))
  }, simplify = FALSE)

  ruler <- dkge:::.dkge_compute_shared_ruler(designs)
  Btil <- dkge:::.dkge_row_standardize(betas, ruler$R)

  # B_tilde_s = R^T B_s for each subject
  for (s in seq_len(S)) {
    expected <- t(ruler$R) %*% betas[[s]]
    expect_equal(Btil[[s]], expected, tolerance = 1e-10)
  }

  # Verify dimensions preserved
  expect_equal(nrow(Btil[[1]]), q)
  expect_equal(ncol(Btil[[1]]), P)
})

# -------------------------------------------------------------------------
# Determinism tests
# -------------------------------------------------------------------------

test_that("dkge_fit is deterministic with same seed (w_method='none')", {
  # First run
  fixture1 <- make_fit_fixture(S = 3, q = 4, P = 5, T = 20, seed = 4001)
  fit1 <- dkge_fit(fixture1$betas, fixture1$designs, K = fixture1$K, rank = 2,
                   w_method = "none")

  # Second run with identical inputs

  fixture2 <- make_fit_fixture(S = 3, q = 4, P = 5, T = 20, seed = 4001)
  fit2 <- dkge_fit(fixture2$betas, fixture2$designs, K = fixture2$K, rank = 2,
                   w_method = "none")

  # All key outputs must be identical
  expect_equal(fit1$U, fit2$U, tolerance = 0)
  expect_equal(fit1$Chat, fit2$Chat, tolerance = 0)
  expect_equal(fit1$R, fit2$R, tolerance = 0)
  expect_equal(fit1$weights, fit2$weights, tolerance = 0)
  expect_equal(fit1$eigenvalues, fit2$eigenvalues, tolerance = 0)
})

test_that("dkge_fit is deterministic with w_method='mfa_sigma1'", {
  # MFA weighting uses power iteration which has random init
  # but the fixture provides consistent data, so results should match

  fixture1 <- make_fit_fixture(S = 4, q = 4, P = 6, T = 20, seed = 4002)
  withr::local_seed(100)
  fit1 <- dkge_fit(fixture1$betas, fixture1$designs, K = fixture1$K, rank = 2,
                   w_method = "mfa_sigma1", w_tau = 0.3)

  fixture2 <- make_fit_fixture(S = 4, q = 4, P = 6, T = 20, seed = 4002)
  withr::local_seed(100)
  fit2 <- dkge_fit(fixture2$betas, fixture2$designs, K = fixture2$K, rank = 2,
                   w_method = "mfa_sigma1", w_tau = 0.3)

  expect_equal(fit1$U, fit2$U, tolerance = 0)
  expect_equal(fit1$Chat, fit2$Chat, tolerance = 0)
  expect_equal(fit1$R, fit2$R, tolerance = 0)
  expect_equal(fit1$weights, fit2$weights, tolerance = 0)
  expect_equal(fit1$eigenvalues, fit2$eigenvalues, tolerance = 0)
})

test_that("dkge_fit is deterministic with w_method='energy'", {
  fixture1 <- make_fit_fixture(S = 4, q = 4, P = 6, T = 20, seed = 4003)
  fit1 <- dkge_fit(fixture1$betas, fixture1$designs, K = fixture1$K, rank = 2,
                   w_method = "energy", w_tau = 0.5)

  fixture2 <- make_fit_fixture(S = 4, q = 4, P = 6, T = 20, seed = 4003)
  fit2 <- dkge_fit(fixture2$betas, fixture2$designs, K = fixture2$K, rank = 2,
                   w_method = "energy", w_tau = 0.5)

  expect_equal(fit1$U, fit2$U, tolerance = 0)
  expect_equal(fit1$Chat, fit2$Chat, tolerance = 0)
  expect_equal(fit1$R, fit2$R, tolerance = 0)
  expect_equal(fit1$weights, fit2$weights, tolerance = 0)
  expect_equal(fit1$eigenvalues, fit2$eigenvalues, tolerance = 0)
})

test_that("dkge_fit with solver='pooled' is deterministic", {
  fixture1 <- make_fit_fixture(S = 3, q = 5, P = 6, T = 20, seed = 4004)
  fit1 <- dkge_fit(fixture1$betas, fixture1$designs, K = fixture1$K, rank = 3,
                   solver = "pooled", w_method = "none")

  fixture2 <- make_fit_fixture(S = 3, q = 5, P = 6, T = 20, seed = 4004)
  fit2 <- dkge_fit(fixture2$betas, fixture2$designs, K = fixture2$K, rank = 3,
                   solver = "pooled", w_method = "none")

  expect_equal(fit1$U, fit2$U, tolerance = 0)
  expect_equal(fit1$Chat, fit2$Chat, tolerance = 0)
  expect_equal(fit1$eigenvalues, fit2$eigenvalues, tolerance = 0)
})

test_that("dkge_fit determinism holds with ridge regularization", {
  fixture1 <- make_fit_fixture(S = 3, q = 4, P = 5, T = 20, seed = 4005)
  fit1 <- dkge_fit(fixture1$betas, fixture1$designs, K = fixture1$K, rank = 2,
                   ridge = 0.1, w_method = "none")

  fixture2 <- make_fit_fixture(S = 3, q = 4, P = 5, T = 20, seed = 4005)
  fit2 <- dkge_fit(fixture2$betas, fixture2$designs, K = fixture2$K, rank = 2,
                   ridge = 0.1, w_method = "none")

  expect_equal(fit1$U, fit2$U, tolerance = 0)
  expect_equal(fit1$Chat, fit2$Chat, tolerance = 0)
  expect_equal(fit1$R, fit2$R, tolerance = 0)
})

test_that("dkge_fit determinism holds with Omega_list spatial weights", {
  fixture1 <- make_fit_fixture(S = 3, q = 4, P = 5, T = 20, seed = 4006)
  set.seed(4006)
  Omega1 <- lapply(fixture1$betas, function(B) runif(ncol(B)))

  fixture2 <- make_fit_fixture(S = 3, q = 4, P = 5, T = 20, seed = 4006)
  set.seed(4006)
  Omega2 <- lapply(fixture2$betas, function(B) runif(ncol(B)))

  fit1 <- dkge_fit(fixture1$betas, fixture1$designs, K = fixture1$K, rank = 2,
                   Omega_list = Omega1, w_method = "none")
  fit2 <- dkge_fit(fixture2$betas, fixture2$designs, K = fixture2$K, rank = 2,
                   Omega_list = Omega2, w_method = "none")

  expect_equal(fit1$U, fit2$U, tolerance = 0)
  expect_equal(fit1$Chat, fit2$Chat, tolerance = 0)
  expect_equal(fit1$R, fit2$R, tolerance = 0)
})

# -------------------------------------------------------------------------
# Partial coverage: Chat, training blocks, and fold re-pooling agree
# -------------------------------------------------------------------------

make_partial_coverage_data <- function(seed = 77) {
  set.seed(seed)
  effects <- c("e1", "e2", "e3", "e4")
  design_full <- qr.Q(qr(matrix(rnorm(24 * 4), 24, 4)))
  colnames(design_full) <- effects
  subs <- list(
    list(id = "s1", rows = c(1L, 2L, 3L)),
    list(id = "s2", rows = c(2L, 3L, 4L)),
    list(id = "s3", rows = c(1L, 3L, 4L)),
    list(id = "s4", rows = 1:4)
  )
  subjects <- lapply(subs, function(spec) {
    B <- matrix(rnorm(4 * 6), 4, 6, dimnames = list(effects, NULL))
    dkge_subject(B, design = design_full, id = spec$id,
                 observed_rows = spec$rows)
  })
  K <- matrix(0.3, 4, 4)
  diag(K) <- 1
  dimnames(K) <- list(effects, effects)
  list(data = dkge_data(subjects), K = K)
}

test_that("Chat equals X_concat %*% t(X_concat) under partial coverage", {
  fx <- make_partial_coverage_data()
  fit <- dkge_fit(fx$data, K = fx$K, rank = 3, keep_X = TRUE,
                  w_method = "energy")

  expect_true(any(!unlist(fx$data$provenance$obs_mask)))  # coverage really partial
  expect_equal(fit$Chat, tcrossprod(fit$X_concat),
               tolerance = 1e-10, ignore_attr = TRUE)

  # Contributions still sum to Chat on the linear (missingness="none",
  # no effect weighting) pooling path.
  contrib_sum <- Reduce(`+`, Map(function(w, S) w * S,
                                 fit$weights, fit$contribs))
  expect_equal(fit$Chat, contrib_sum, tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("training blocks match dkge_transform_block under partial coverage", {
  fx <- make_partial_coverage_data()
  fit <- dkge_fit(fx$data, K = fx$K, rank = 3, keep_X = TRUE,
                  w_method = "energy")

  for (s in seq_along(fit$Btil)) {
    idx <- fit$block_indices[[s]]
    expect_equal(
      fit$X_concat[, idx, drop = FALSE],
      dkge_transform_block(fit, fx$data$betas[[s]],
                           Omega_s = fit$Omega[[s]], w_s = fit$weights[s]),
      tolerance = 1e-10, ignore_attr = TRUE
    )
  }
})

test_that(".dkge_fold_weight_context on the full training set reproduces fit$Chat", {
  fx <- make_partial_coverage_data()
  for (miss in c("none", "rescale", "mask", "shrink")) {
    fit <- dkge_fit(fx$data, K = fx$K, rank = 3, w_method = "energy",
                    missingness = miss,
                    miss_args = list(min_pairs = 2L, gamma = 1))
    ctx <- dkge:::.dkge_fold_weight_context(fit, seq_along(fit$Btil))
    expect_equal(ctx$Chat, fit$Chat, tolerance = 1e-12, ignore_attr = TRUE,
                 info = miss)
    expect_equal(ctx$pair_counts, fit$pair_counts, ignore_attr = TRUE,
                 info = miss)
  }
})

test_that("fold re-pooling agrees whether or not stored moments are reused", {
  fx <- make_partial_coverage_data()
  fit <- dkge_fit(fx$data, K = fx$K, rank = 3, w_method = "energy")
  train <- c(1L, 2L, 4L)

  # The fast path must actually be selected here.
  expect_true(dkge:::.dkge_voxel_weights_match(fit, fit$voxel_weights, train))

  fast <- dkge:::.dkge_fold_weight_context(fit, train)
  slow_fit <- fit
  slow_fit$effect_moments <- NULL  # force the recompute-from-betas path
  slow <- dkge:::.dkge_fold_weight_context(slow_fit, train)

  expect_equal(fast$Chat, slow$Chat, tolerance = 1e-12)
  expect_equal(fast$pair_counts, slow$pair_counts)
  expect_equal(fast$pair_ess, slow$pair_ess, tolerance = 1e-12)
})

test_that("design kernel labels are validated against the data effects", {
  fx <- make_partial_coverage_data()
  effects <- fx$data$effects

  K_bad <- fx$K
  dimnames(K_bad) <- list(c("e1", "e2", "e3", "zz"), c("e1", "e2", "e3", "zz"))
  expect_error(dkge_fit(fx$data, K = K_bad, rank = 2, w_method = "none"),
               "do not match the data effects")

  perm <- c(3L, 1L, 4L, 2L)
  K_perm <- fx$K[perm, perm, drop = FALSE]
  expect_message(
    fit_perm <- dkge_fit(fx$data, K = K_perm, rank = 2, w_method = "none"),
    "reordered"
  )
  fit_ref <- dkge_fit(fx$data, K = fx$K, rank = 2, w_method = "none")
  expect_equal(fit_perm$K, fit_ref$K, ignore_attr = TRUE)
  expect_equal(fit_perm$Chat, fit_ref$Chat, tolerance = 1e-12)
  expect_equal(colnames(fit_perm$K), effects)
})

test_that("wrong-length observation masks are rejected instead of ignored", {
  expect_error(dkge:::.dkge_observed_rows(c(TRUE, FALSE), 3L),
               "one entry per effect row")
  expect_error(
    dkge:::.dkge_obs_masks_from_provenance(
      list(obs_mask = list(a = c(TRUE, TRUE))), c("a"), 3L
    ),
    "length 2 but the fit has 3"
  )
})

test_that("kernel with duplicated labels warns instead of silently ignoring names", {
  set.seed(5)
  q <- 4
  betas <- replicate(3, matrix(rnorm(q * 10), q, 10,
                               dimnames = list(paste0("e", 1:q), NULL)),
                     simplify = FALSE)
  designs <- replicate(3, {
    X <- matrix(rnorm(20 * q), 20, q); colnames(X) <- paste0("e", 1:q); X
  }, simplify = FALSE)
  K <- diag(q); dimnames(K) <- list(c("e1", "e1", "e3", "e4"), c("e1", "e1", "e3", "e4"))
  expect_warning(dkge_fit(betas, designs, K = K, rank = 2), "duplicates")
})

test_that("miss_args rejects unknown field names", {
  fx <- make_fit_fixture()
  expect_error(
    dkge_fit(fx$betas, fx$designs, K = fx$K, rank = 2,
             missingness = "shrink", miss_args = list(gama = 3)),
    "gama"
  )
  expect_silent(
    dkge_fit(fx$betas, fx$designs, K = fx$K, rank = 2,
             missingness = "shrink", miss_args = list(gamma = 3))
  )
})

test_that("split-half error names the offending subject", {
  effects <- c("e1", "e2")
  B <- matrix(1:4, 2, 2, dimnames = list(effects, NULL))
  X <- diag(2)
  colnames(X) <- effects
  alice <- suppressWarnings(dkge_subject(B, X, id = "alice"))
  bob <- suppressWarnings(dkge_subject(B, X, id = "bob"))
  data <- dkge_data(list(alice, bob))
  expect_error(
    dkge_fit(data, K = diag(2), rank = 1, debias = "split_half",
             w_method = "none", effect_scaling = "none"),
    "alice"
  )
})

