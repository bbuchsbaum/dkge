# test-plot-helpers.R
# Focused coverage for individual plotting utilities

library(testthat)

set.seed(314)

# Build a fit whose beta rows / design columns carry the kernel's effect labels.
make_kernel_fit <- function(dk, S = 4, P = 8, rank = 3, seed = 11) {
  set.seed(seed)
  labels <- rownames(dk$K)
  q <- nrow(dk$K)
  betas <- replicate(S, {
    B <- matrix(rnorm(q * P), q, P)
    rownames(B) <- labels
    B
  }, simplify = FALSE)
  designs <- replicate(S, {
    X <- diag(q)
    dimnames(X) <- list(labels, labels)
    X
  }, simplify = FALSE)
  dkge(betas, designs, K = dk, rank = rank)
}

test_that("plot helpers return ggplot objects", {
  skip_if_not_installed("ggplot2")

  fx <- make_small_fit(S = 3, q = 4, P = 5, T = 20, rank = 3, seed = 7)
  fit <- fx$fit

  th <- theme_dkge()
  expect_s3_class(th, "theme")

  scree <- dkge_plot_scree(fit, one_se_pick = 1)
  expect_s3_class(scree, "ggplot")

  effects_plot <- dkge_plot_effect_loadings(fit, comps = 1:2)
  expect_s3_class(effects_plot, "ggplot")

  sal <- dkge_component_saliences(fit, comps = 1:2, long = FALSE)
  expect_equal(dim(sal), c(nrow(fit$U), 2L))

  raw_saliences <- dkge_plot_component_saliences(fit, comps = 1:2)
  expect_s3_class(raw_saliences, "ggplot")

  contrast_scores <- dkge_component_contrast_scores(fit, comps = 1:2)
  expect_equal(nrow(contrast_scores), nrow(fit$U) * 2L)

  contrast_plot <- dkge_plot_component_contrast_scores(fit, comps = 1:2)
  expect_s3_class(contrast_plot, "ggplot")

  pooled_proj <- dkge_subject_component_projections(
    fit,
    groups = c("g1", "g1", "g2"),
    mode = "pooled",
    comps = 1:2
  )
  expect_equal(nrow(pooled_proj), 3L * 2L)

  loso_proj <- dkge_subject_component_projections(
    fit,
    groups = c("g1", "g1", "g2"),
    mode = "loso",
    comps = 1:2
  )
  expect_equal(nrow(loso_proj), 3L * 2L)

  subject_plot <- dkge_plot_subject_component_projections(
    fit,
    groups = c("g1", "g1", "g2"),
    mode = "pooled",
    comps = 1:2
  )
  expect_s3_class(subject_plot, "ggplot")

  contrib <- dkge_plot_subject_contrib(fit, comps = 1:2)
  expect_true(is.list(contrib))
  expect_s3_class(contrib$weights, "ggplot")
  expect_s3_class(contrib$energy, "ggplot")

  bases <- list(
    fit$U,
    dkge_k_orthonormalize(fit$U + matrix(1e-4, nrow(fit$U), ncol(fit$U)), fit$K)
  )
  stability <- dkge_plot_subspace_stability(bases, fit$K)
  expect_s3_class(stability, "ggplot")

  info_panels <- dkge_plot_info_anchor(
    info_haufe = list(anchor = c(0.4, -0.2, 0.1)),
    info_loco = list(delta = c(0.05, 0, -0.03)),
    top = 0
  )
  expect_true(all(vapply(info_panels, inherits, logical(1), what = "ggplot")))
})

test_that("component saliences are exactly K %*% U and unit scaling is K-unit", {
  fx <- make_small_fit(S = 3, q = 4, P = 5, T = 20, rank = 3, seed = 7)
  fit <- fx$fit

  sal <- dkge_component_saliences(fit, comps = c(1, 3), long = FALSE)
  expect_equal(unname(sal), unname(fit$K %*% fit$U[, c(1, 3), drop = FALSE]))
  expect_equal(colnames(sal), c("LV1", "LV3"))

  # U is K-orthonormal, so saliences already have unit K-norm: "unit" is a no-op.
  sal_unit <- dkge_component_saliences(fit, comps = c(1, 3), scale = "unit", long = FALSE)
  expect_equal(sal_unit, sal, tolerance = 1e-10)

  long <- dkge_component_saliences(fit, comps = 2, long = TRUE)
  expect_equal(long$salience, as.numeric(fit$K %*% fit$U[, 2]))
  expect_equal(long$component_id, rep(2L, nrow(long)))
})

test_that("component indexer rejects out-of-range components", {
  fx <- make_small_fit(S = 3, q = 4, P = 5, T = 20, rank = 2, seed = 7)
  fit <- fx$fit

  expect_error(dkge_component_saliences(fit, comps = 5), "out of range")
  expect_error(dkge_component_saliences(fit, comps = 0), "out of range")
  expect_error(dkge_component_contrast_scores(fit, comps = c(1, 9)), "out of range")
  expect_error(dkge_subject_component_projections(fit, comps = 3, mode = "pooled"),
               "out of range")
})

test_that("dkge_design_basis columns have unit K-norm by default", {
  factors <- list(task = list(L = 2), measure = list(L = 3))
  dk <- design_kernel(factors, basis = "cell", normalize = "none")
  fit <- make_kernel_fit(dk, S = 4, P = 8, rank = 3, seed = 11)

  C <- dkge_design_basis(fit)
  expect_equal(unname(diag(t(C) %*% fit$K %*% C)), rep(1, ncol(C)), tolerance = 1e-10)

  C_l2 <- dkge_design_basis(fit, normalize = "unit_l2")
  expect_equal(unname(sqrt(colSums(C_l2 * C_l2))), rep(1, ncol(C_l2)), tolerance = 1e-10)

  C_raw <- dkge_design_basis(fit, normalize = "none")
  # Same directions, different scaling: each column is a positive multiple.
  ratios <- unname(C[, 2] / C_raw[, 2])
  expect_equal(ratios, rep(ratios[[1]], length(ratios)), tolerance = 1e-10)
  expect_true(ratios[[1]] > 0)

  # Term labels come from the model-matrix `assign` map, not name prefixes.
  expect_equal(attr(C, "term")[1], "grand_mean")
  expect_true(all(attr(C, "term") %in% c("grand_mean", "task", "measure", "task:measure")))
})

test_that("design basis follows the declared level order stored in the cells", {
  # Level labels deliberately out of alphabetical order. `design_kernel()`
  # records level labels only through its cell table (`info$levels` holds level
  # *counts*), so the coding must be read off the cells' first-appearance
  # order, which is the declared order by construction of the grid.
  factors <- list(task = list(L = 2, levels = c("go", "await")),
                  measure = list(L = 3, levels = c("zebra", "apple", "mango")))
  dk <- design_kernel(factors, basis = "cell", normalize = "none")
  fit <- make_kernel_fit(dk, S = 3, P = 8, rank = 2, seed = 14)

  cells <- as.data.frame(fit$kernel_info$cells, stringsAsFactors = FALSE)
  expect_equal(unique(as.character(cells$measure)), c("zebra", "apple", "mango"))
  # The kernel really does not carry level labels for us to read instead.
  expect_equal(fit$kernel_info$levels, c(task = 2L, measure = 3L))

  C <- dkge_design_basis(fit, normalize = "none")
  # contr.sum on 3 levels: first declared level -> +1, second -> 0, last -> -1.
  expect_true(all(C[cells$measure == "zebra", "measure1"] > 0))
  expect_true(all(C[cells$measure == "apple", "measure1"] == 0))
  expect_true(all(C[cells$measure == "mango", "measure1"] < 0))

  # Reversing the declaration reverses which cells anchor the coding.
  factors_rev <- factors
  factors_rev$measure$levels <- c("mango", "apple", "zebra")
  dk_rev <- design_kernel(factors_rev, basis = "cell", normalize = "none")
  fit_rev <- make_kernel_fit(dk_rev, S = 3, P = 8, rank = 2, seed = 14)
  cells_rev <- as.data.frame(fit_rev$kernel_info$cells, stringsAsFactors = FALSE)
  C_rev <- dkge_design_basis(fit_rev, normalize = "none")
  expect_true(all(C_rev[cells_rev$measure == "mango", "measure1"] > 0))
  expect_true(all(C_rev[cells_rev$measure == "zebra", "measure1"] < 0))
})

test_that("the identity fallback basis honours `normalize`", {
  # An effect-basis kernel has cell metadata whose row count does not match the
  # kernel dimension, so `dkge_design_basis()` falls back to the identity basis
  # over effect rows. That fallback is still a default basis and must obey the
  # documented `unit_K` default.
  factors <- list(A = list(L = 2), B = list(L = 2), time = list(L = 4))
  dk <- design_kernel(factors, basis = "effect")
  fit <- make_kernel_fit(dk, S = 3, P = 9, rank = 3, seed = 31)
  expect_false(nrow(as.data.frame(fit$kernel_info$cells)) == nrow(fit$K))

  C <- dkge_design_basis(fit)
  expect_equal(attr(C, "source"), "identity")
  expect_equal(unname(diag(t(C) %*% fit$K %*% C)), rep(1, ncol(C)),
               tolerance = 1e-10)

  C_raw <- dkge_design_basis(fit, normalize = "none")
  # The unnormalized fallback is genuinely off unit K-norm, so the check above
  # is not vacuous.
  expect_false(isTRUE(all.equal(unname(diag(t(C_raw) %*% fit$K %*% C_raw)),
                                rep(1, ncol(C_raw)))))

  C_l2 <- dkge_design_basis(fit, normalize = "unit_l2")
  expect_equal(unname(sqrt(colSums(C_l2 * C_l2))), rep(1, ncol(C_l2)),
               tolerance = 1e-10)
})

test_that("design basis term labels survive factor names that share a prefix", {
  factors <- list(task = list(L = 2), task2 = list(L = 3))
  dk <- design_kernel(factors, basis = "cell", normalize = "none")
  fit <- make_kernel_fit(dk, S = 3, P = 8, rank = 2, seed = 12)

  C <- dkge_design_basis(fit)
  terms <- attr(C, "term")
  expect_true(all(terms %in% c("grand_mean", "task", "task2", "task:task2")))
  expect_equal(sum(terms == "task"), 1L)
  expect_equal(sum(terms == "task2"), 2L)
  expect_equal(sum(terms == "task:task2"), 2L)
})

test_that("subject projections are signed and match a hand-computed reference", {
  fx <- make_small_fit(S = 4, q = 4, P = 6, T = 20, rank = 3, seed = 21)
  fit <- fx$fit
  subjects <- fit$subject_ids

  proj <- dkge_subject_component_projections(fit, mode = "pooled", comps = 1:2)
  expect_equal(nrow(proj), 4L * 2L)

  # Reference: <K u_j, cluster-mean of the subject's whitened betas>.
  for (s in seq_len(4)) {
    for (j in 1:2) {
      expected <- sum((fit$K %*% fit$U[, j]) * rowMeans(fit$Btil[[s]]))
      got <- proj$projection[proj$subject == subjects[[s]] & proj$component_id == j]
      expect_equal(got, expected)
    }
  }

  # Signed: negating one subject's betas flips that subject's scores only.
  flipped <- fit
  flipped$Btil[[1]] <- -flipped$Btil[[1]]
  proj_flip <- dkge_subject_component_projections(flipped, mode = "pooled", comps = 1:2)
  first <- proj$subject == subjects[[1]]
  expect_equal(proj_flip$projection[first], -proj$projection[first])
  expect_equal(proj_flip$projection[!first], proj$projection[!first])

  # Not the sign-blind energy measure.
  energy <- dkge:::.dkge_energy_by_subject_component(fit)
  expect_false(isTRUE(all.equal(proj$projection[proj$component_id == 1],
                                unname(energy[, 1]))))
  expect_true(any(proj$projection < 0))
})

test_that("LOSO subject projections use the held-out basis and stay signed", {
  fx <- make_small_fit(S = 4, q = 4, P = 6, T = 20, rank = 3, seed = 22)
  fit <- fx$fit
  subjects <- fit$subject_ids

  loso <- dkge_subject_component_projections(fit, mode = "loso", comps = 1:2)
  expect_equal(nrow(loso), 4L * 2L)

  folds <- dkge:::.dkge_build_fold_bases(fit,
                                         assignments = lapply(seq_len(4), function(s) s),
                                         ridge = 0,
                                         align = TRUE,
                                         loader_scope = "heldout",
                                         verbose = FALSE)
  for (fold in folds$folds) {
    s <- fold$subjects[[1]]
    # Pin the rotation convention: the aligned basis is the fold basis times
    # the returned rotation, so a transposed rotation would be caught here
    # rather than silently cancelling out below (rank 3, where it is visible).
    expect_equal(fold$basis_aligned, fold$basis %*% fold$rotation)

    # The documented reference is the *pooled* basis, so the rotation applied
    # to the loader scores is a direct K-Procrustes of the fold basis onto
    # fit$U -- not the fold builder's rotation, which references fold 1.
    R_pooled <- dkge_procrustes_K(fit$U, fold$basis, fit$K)$R
    basis_pooled <- fold$basis %*% R_pooled
    A_aligned <- fold$loaders[[as.character(s)]]$A %*% R_pooled
    for (j in 1:2) {
      got <- loso$projection[loso$subject == subjects[[s]] & loso$component_id == j]
      expect_equal(got, mean(A_aligned[, j]))
      # Same value computed straight from the pooled-aligned held-out basis,
      # with no reference to how the rotation is represented.
      scores <- as.numeric(t(fit$Btil[[s]]) %*% fit$K %*% basis_pooled[, j])
      expect_equal(got, mean(scores))
    }
  }

  # The two references really are different: aligning fold 1 to fold 1 is the
  # identity, so the fold builder's rotation cannot be standing in for a
  # rotation onto the pooled axes.
  fold1 <- folds$folds[[1]]
  expect_equal(fold1$rotation, diag(ncol(fit$U)))
  expect_false(isTRUE(all.equal(dkge_procrustes_K(fit$U, fold1$basis, fit$K)$R,
                                diag(ncol(fit$U)))))

  # align = FALSE leaves the fold's own eigen-basis untouched.
  loso_raw <- dkge_subject_component_projections(fit, mode = "loso", comps = 1:2,
                                                 align = FALSE)
  for (fold in folds$folds) {
    s <- fold$subjects[[1]]
    A_raw <- fold$loaders[[as.character(s)]]$A
    for (j in 1:2) {
      got <- loso_raw$projection[loso_raw$subject == subjects[[s]] &
                                   loso_raw$component_id == j]
      expect_equal(got, mean(A_raw[, j]))
    }
  }

  flipped <- fit
  flipped$Btil[[2]] <- -flipped$Btil[[2]]
  loso_flip <- dkge_subject_component_projections(flipped, mode = "loso", comps = 1:2)
  second <- loso$subject == subjects[[2]]
  expect_equal(loso_flip$projection[second], -loso$projection[second])
})

test_that("subject projection groups are matched by name, not position", {
  fx <- make_small_fit(S = 3, q = 3, P = 5, T = 20, rank = 2, seed = 23)
  fit <- fx$fit
  subjects <- fit$subject_ids

  named <- stats::setNames(c("late", "mid", "early"), rev(subjects))
  proj <- dkge_subject_component_projections(fit, groups = named,
                                             mode = "pooled", comps = 1)
  expect_equal(as.character(proj$group),
               unname(named[proj$subject]))
  expect_equal(as.character(proj$group), c("early", "mid", "late"))

  df_groups <- data.frame(subject = rev(subjects),
                          group = c("late", "mid", "early"),
                          stringsAsFactors = FALSE)
  proj_df <- dkge_subject_component_projections(fit, groups = df_groups,
                                                mode = "pooled", comps = 1)
  expect_equal(as.character(proj_df$group), c("early", "mid", "late"))

  # Unnamed vectors stay positional.
  proj_pos <- dkge_subject_component_projections(fit, groups = c("a", "b", "c"),
                                                 mode = "pooled", comps = 1)
  expect_equal(as.character(proj_pos$group), c("a", "b", "c"))

  # Incomplete coverage is an error, not silent NA.
  expect_error(
    dkge_subject_component_projections(fit,
                                       groups = stats::setNames(c("a", "b"), subjects[1:2]),
                                       mode = "pooled", comps = 1),
    "do not cover subject"
  )
  expect_error(
    dkge_subject_component_projections(fit,
                                       groups = data.frame(subject = subjects[1:2],
                                                           group = c("a", "b")),
                                       mode = "pooled", comps = 1),
    "does not cover subject"
  )
})

test_that("subject-projection voxel weights follow the fit-time policy", {
  fx <- make_small_fit(S = 3, q = 3, P = 6, T = 20, rank = 2, seed = 25)
  fit <- fx$fit

  # Length-3 varying weights against P = 6 clusters no longer recycle: the same
  # `.dkge_subject_voxel_weights()` policy dkge_fit() uses rejects stretching
  # one subject's spatial profile onto another's parcels.
  fit$voxel_weights <- c(1, 2, 3)
  expect_error(
    dkge_subject_component_projections(fit, mode = "pooled", comps = 1),
    "do not apply to subject"
  )

  # A constant vector carries no per-column information and is re-expanded.
  fit$voxel_weights <- rep(2, 3)
  expect_warning(
    dkge_subject_component_projections(fit, mode = "pooled", comps = 1),
    NA
  )

  # A scalar weight is a legitimate broadcast, not a recycling accident.
  fit$voxel_weights <- 2
  expect_warning(
    dkge_subject_component_projections(fit, mode = "pooled", comps = 1),
    NA
  )

  # Per-subject lists are indexed, never recycled across subjects.
  fit$voxel_weights <- NULL
  fit$voxel_weights_subject <- list(rep(1, 6), c(3, rep(1, 5)), rep(1, 6))
  proj <- dkge_subject_component_projections(fit, mode = "pooled", comps = 1)
  w2 <- c(3, rep(1, 5))
  expected2 <- sum((fit$K %*% fit$U[, 1]) *
                     rowMeans(sweep(fit$Btil[[2]], 2L, sqrt(w2 / mean(w2)), "*")))
  expect_equal(proj$projection[[2]], expected2)
})

test_that("subject projections are invariant to the voxel-weight scale", {
  fx <- make_small_fit(S = 3, q = 3, P = 6, T = 20, rank = 2, seed = 25)
  fit <- fx$fit
  base <- dkge_subject_component_projections(fit, mode = "pooled", comps = 1:2)

  # A uniform weight only rescales every cluster equally, so it must not move
  # the score at all -- previously a uniform weight of 4 doubled it.
  for (scale in c(0.25, 1, 4, 100)) {
    weighted <- fit
    weighted$voxel_weights <- rep(scale, ncol(fit$Btil[[1]]))
    got <- dkge_subject_component_projections(weighted, mode = "pooled", comps = 1:2)
    expect_equal(got$projection, base$projection)
  }

  # A non-uniform weight does change the score, but only through its shape.
  w <- c(1, 2, 3, 1, 2, 3)
  shaped <- fit
  shaped$voxel_weights <- w
  scaled <- fit
  scaled$voxel_weights <- 10 * w
  p_shaped <- dkge_subject_component_projections(shaped, mode = "pooled", comps = 1)
  p_scaled <- dkge_subject_component_projections(scaled, mode = "pooled", comps = 1)
  expect_equal(p_shaped$projection, p_scaled$projection)
  expect_false(isTRUE(all.equal(p_shaped$projection,
                                base$projection[base$component_id == 1])))

  expected <- sum((fit$K %*% fit$U[, 1]) *
                    rowMeans(sweep(fit$Btil[[1]], 2L, sqrt(w / mean(w)), "*")))
  expect_equal(p_shaped$projection[[1]], expected)
})

test_that("plotting subject projections accepts a precomputed data frame", {
  skip_if_not_installed("ggplot2")
  fx <- make_small_fit(S = 3, q = 3, P = 5, T = 20, rank = 2, seed = 24)
  fit <- fx$fit

  proj <- dkge_subject_component_projections(fit, mode = "pooled", comps = 1:2)
  p <- dkge_plot_subject_component_projections(fit, projections = proj)
  expect_s3_class(p, "ggplot")
  expect_equal(nrow(p$data), nrow(proj))
  expect_error(
    dkge_plot_subject_component_projections(fit, projections = proj[, c("subject", "group")]),
    "must contain columns"
  )

  # `fit` is unnecessary when projections are supplied, but one of the two must be.
  expect_s3_class(dkge_plot_subject_component_projections(projections = proj), "ggplot")
  expect_error(dkge_plot_subject_component_projections(),
               "Supply either `fit` or a precomputed")
})

test_that("supplied projections are labelled by their own mode column", {
  skip_if_not_installed("ggplot2")
  fx <- make_small_fit(S = 3, q = 3, P = 5, T = 20, rank = 2, seed = 24)
  fit <- fx$fit

  pooled <- dkge_subject_component_projections(fit, mode = "pooled", comps = 1:2)
  p_pooled <- dkge_plot_subject_component_projections(projections = pooled)
  expect_match(p_pooled$labels$title, "POOLED", fixed = TRUE)

  loso <- dkge_subject_component_projections(fit, mode = "loso", comps = 1:2)
  p_loso <- dkge_plot_subject_component_projections(projections = loso)
  expect_match(p_loso$labels$title, "LOSO", fixed = TRUE)

  # A frame without a `mode` column carries no claim about the basis, so it
  # must not inherit the "loso" default of the ignored `mode` argument.
  bare <- pooled
  bare$mode <- NULL
  p_bare <- dkge_plot_subject_component_projections(projections = bare)
  expect_false(grepl("LOSO", p_bare$labels$title, fixed = TRUE))
  expect_false(grepl("POOLED", p_bare$labels$title, fixed = TRUE))
  expect_match(p_bare$labels$subtitle, "not recorded")

  # A frame mixing modes is equally unlabelled.
  mixed <- rbind(pooled, loso)
  p_mixed <- dkge_plot_subject_component_projections(projections = mixed)
  expect_false(grepl("LOSO", p_mixed$labels$title, fixed = TRUE))

  # The explicit `mode` argument cannot override the frame.
  p_forced <- dkge_plot_subject_component_projections(projections = pooled,
                                                      mode = "loso")
  expect_match(p_forced$labels$title, "POOLED", fixed = TRUE)

  # Required-column and type validation.
  expect_error(
    dkge_plot_subject_component_projections(projections = pooled[0, , drop = FALSE]),
    "at least one row"
  )
  bad <- pooled
  bad$projection <- as.character(bad$projection)
  expect_error(dkge_plot_subject_component_projections(projections = bad),
               "must be numeric")
})

test_that("contrast basis defaults to factorial coordinates for cell kernels", {
  factors <- list(task = list(L = 2), measure = list(L = 3))
  dk <- design_kernel(factors, basis = "cell", normalize = "none")
  q <- nrow(dk$K)
  fit <- make_kernel_fit(dk, S = 4, P = 8, rank = 3, seed = 13)

  basis <- dkge_design_basis(fit)

  expect_equal(nrow(basis), q)
  expect_true("grand_mean" %in% colnames(basis))
  expect_true("task" %in% attr(basis, "term"))
  expect_true("measure" %in% attr(basis, "term"))
  expect_true("task:measure" %in% attr(basis, "term"))

  scores <- dkge_component_contrast_scores(fit, comps = 1:2)
  expect_equal(nrow(scores), ncol(basis) * 2L)
})

# Restore kernel_info to the original design_kernel() cell order while leaving
# fit$effects / K / U in the permuted data order. Consumers must rematch by
# name; positional pairing of cells against effects is wrong on every row.
.restore_kernel_order_info <- function(fit, fx) {
  info <- fx$kernel$info
  fit$kernel_info$cell_labels <- info$cell_labels
  fit$kernel_info$cells <- info$cells
  fit$kernel_info$blocks <- info$blocks
  fit
}

.oracle_cells_from_effects <- function(effects) {
  parts <- do.call(rbind, strsplit(effects, ":", fixed = TRUE))
  data.frame(
    task = parts[, 1],
    measure = parts[, 2],
    stringsAsFactors = FALSE
  )
}

.oracle_task_measure_basis <- function(effects, K, normalize = "unit_K") {
  cells <- .oracle_cells_from_effects(effects)
  for (nm in c("task", "measure")) {
    cells[[nm]] <- factor(cells[[nm]], levels = unique(as.character(cells[[nm]])))
  }
  mm <- stats::model.matrix(
    ~ task * measure,
    data = cells,
    contrasts.arg = list(
      task = stats::contr.sum(nlevels(cells$task)),
      measure = stats::contr.sum(nlevels(cells$measure))
    )
  )
  rownames(mm) <- effects
  if (identical(normalize, "none")) return(mm)
  norms <- sqrt(pmax(colSums(mm * (K %*% mm)), 0))
  sweep(mm, 2L, norms, "/")
}

test_that("design basis rematches permuted cell metadata by name", {
  fx <- make_permuted_cell_fit(seed = 7)
  fit <- .restore_kernel_order_info(fx$fit, fx)
  expect_false(identical(fit$kernel_info$cell_labels, fit$effects))

  C <- dkge_design_basis(fit, normalize = "none")
  mm <- .oracle_task_measure_basis(fit$effects, fit$K, normalize = "none")
  expect_equal(unname(as.matrix(C)), unname(as.matrix(mm)),
               tolerance = 1e-10, ignore_attr = TRUE)
  expect_equal(rownames(C), fit$effects)

  C_k <- dkge_design_basis(fit)
  mm_k <- .oracle_task_measure_basis(fit$effects, fit$K, normalize = "unit_K")
  expect_equal(unname(as.matrix(C_k)), unname(as.matrix(mm_k)),
               tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("component contrast scores use the name-matched task column", {
  fx <- make_permuted_cell_fit(seed = 7)
  fit <- .restore_kernel_order_info(fx$fit, fx)
  scores <- dkge_component_contrast_scores(fit, comps = 1:3)
  task <- scores[scores$term == "task", , drop = FALSE]
  expect_equal(nrow(task), 3L)

  C <- .oracle_task_measure_basis(fit$effects, fit$K, normalize = "unit_K")
  assign <- attr(
    stats::model.matrix(
      ~ task * measure,
      data = {
        cells <- .oracle_cells_from_effects(fit$effects)
        for (nm in c("task", "measure")) {
          cells[[nm]] <- factor(cells[[nm]],
                                levels = unique(as.character(cells[[nm]])))
        }
        cells
      },
      contrasts.arg = list(task = stats::contr.sum(2), measure = stats::contr.sum(3))
    ),
    "assign"
  )
  C_task <- C[, assign == 1L, drop = FALSE]
  oracle <- as.numeric(crossprod(C_task, fit$K %*% fit$U))
  expect_equal(task$score, oracle, tolerance = 1e-10)
})

test_that("design basis quotes non-syntactic factor names and drops 1-level factors", {
  fx <- make_permuted_cell_fit(seed = 3)
  fit <- fx$fit

  info <- fit$kernel_info
  names(info$cells)[names(info$cells) == "task"] <- "my factor"
  info$factor_names[info$factor_names == "task"] <- "my factor"
  names(info$factor_scope)[names(info$factor_scope) == "task"] <- "my factor"
  fit$kernel_info <- info
  C_space <- dkge_design_basis(fit, normalize = "none")
  expect_true("my factor" %in% attr(C_space, "term"))
  expect_equal(attr(C_space, "source"), "factorial")

  fit_one <- fx$fit
  fit_one$kernel_info$cells$measure <- "only"
  C_drop <- dkge_design_basis(fit_one, normalize = "none")
  expect_false("measure" %in% attr(C_drop, "term"))
  expect_true("task" %in% attr(C_drop, "term"))
  expect_equal(attr(C_drop, "source"), "factorial")

  fit_none <- fx$fit
  fit_none$kernel_info$cells$task <- "only"
  fit_none$kernel_info$cells$measure <- "only"
  C_id <- dkge_design_basis(fit_none, normalize = "none")
  expect_equal(attr(C_id, "source"), "identity")
  expect_equal(unname(as.matrix(C_id)), diag(length(fit_none$effects)),
               ignore_attr = TRUE)
})

test_that("unmatched cell_labels fall back to the identity basis with a message", {
  fx <- make_permuted_cell_fit(seed = 3)
  fit <- fx$fit
  fit$kernel_info$cell_labels <- paste0("other", seq_along(fit$effects))
  expect_message(
    C <- dkge_design_basis(fit, normalize = "none"),
    "cell_labels"
  )
  expect_equal(attr(C, "source"), "identity")
})

test_that("user term/source attributes survive design-basis validation", {
  fx <- make_small_fit(S = 3, q = 4, P = 5, T = 20, rank = 2, seed = 7)
  fit <- fx$fit
  M <- diag(nrow(fit$U))
  rownames(M) <- fit$effects
  attr(M, "term") <- paste0("user", seq_len(ncol(M)))
  attr(M, "source") <- "custom_user"
  C <- dkge_design_basis(fit, basis = M)
  expect_equal(attr(C, "term"), attr(M, "term"))
  expect_equal(attr(C, "source"), "custom_user")

  expect_error(
    dkge_design_basis(fit, basis = list(basis_matrix = M)),
    "must have one row"
  )
})

test_that("projection and energy both respond to voxel weights", {
  set.seed(31)
  q <- 3; P <- 4; S <- 3
  effects <- paste0("e", seq_len(q))
  betas <- replicate(S, {
    B <- matrix(rnorm(q * P), q, P)
    dimnames(B) <- list(effects, paste0("v", seq_len(P)))
    B
  }, simplify = FALSE)
  designs <- replicate(S, { X <- diag(q); colnames(X) <- effects; X },
                       simplify = FALSE)
  prior <- c(0.2, 1, 3, 0.5)
  fit <- dkge_fit(betas, designs, K = diag(q), rank = 2, w_method = "none",
                  weights = dkge_weights(prior = prior, adapt = "none"))
  proj_w <- dkge_subject_component_projections(fit, mode = "pooled", comps = 1:2)
  energy_w <- dkge:::.dkge_energy_by_subject_component(fit)

  fit$voxel_weights <- rep(1, P)
  fit$voxel_weights_subject <- NULL
  proj_u <- dkge_subject_component_projections(fit, mode = "pooled", comps = 1:2)
  energy_u <- dkge:::.dkge_energy_by_subject_component(fit)

  expect_false(isTRUE(all.equal(proj_w$projection, proj_u$projection)))
  expect_false(isTRUE(all.equal(unname(energy_w), unname(energy_u))))
})

test_that("subject-contrib panels share subject order and block separators follow min(idx)", {
  skip_if_not_installed("ggplot2")
  fx <- make_small_fit(S = 3, q = 4, P = 5, T = 20, rank = 2, seed = 7)
  fit <- fx$fit
  fit$subject_ids <- c("zeta", "alpha", "mu")
  panels <- dkge_plot_subject_contrib(fit, comps = 1:2)
  w_build <- ggplot2::ggplot_build(panels$weights)
  e_build <- ggplot2::ggplot_build(panels$energy)
  w_x <- as.character(w_build$layout$panel_params[[1]]$x$get_labels())
  e_y <- as.character(e_build$layout$panel_params[[1]]$y$get_labels())
  expect_equal(w_x, fit$subject_ids)
  expect_equal(e_y, fit$subject_ids)

  factors <- list(A = list(L = 2), B = list(L = 3))
  dk <- design_kernel(factors, basis = "effect", normalize = "none")
  labels <- rownames(dk$K)
  q <- length(labels)
  set.seed(41)
  betas <- replicate(3, {
    B <- matrix(rnorm(q * 6), q, 6)
    rownames(B) <- labels
    B
  }, simplify = FALSE)
  designs <- replicate(3, {
    X <- diag(q)
    dimnames(X) <- list(labels, labels)
    X
  }, simplify = FALSE)
  fit_ab <- dkge(betas, designs, K = dk, rank = 2)
  # Reverse block list order so positional cumsum is wrong.
  fit_ab$kernel_info$blocks <- rev(fit_ab$kernel_info$blocks)
  p <- dkge_plot_effect_loadings(fit_ab, comps = 1)
  built <- ggplot2::ggplot_build(p)
  hlines <- built$data[[2]]$yintercept
  # Effect-basis A (1) + B (2) + A:B (2), ordered by min(idx): 1.5 and 3.5.
  expect_equal(sort(hlines), c(1.5, 3.5))
})

