# Reference transcription of the pre-`rowsum()` `.dkge_aggregate_from_spec()`
# body: one data-frame row per subject-cell, rbind, then a per-aggregate-row
# string-mask scan. Kept here so the vectorised implementation stays pinned to
# the semantics it replaced. Reads values out of the stacked `spec$V`.
old_aggregate_from_spec <- function(spec,
                                    subject_data = spec$subject_data,
                                    subject_indices = seq_along(spec$subject_ids)) {
  subject_indices <- as.integer(subject_indices)
  subject_data <- as.data.frame(subject_data, stringsAsFactors = FALSE)
  row_vars <- c(spec$group_vars, spec$cell_vars)
  src_value <- function(src, cell_idx) spec$V[(src - 1L) * spec$q + cell_idx, ]

  source_rows <- vector("list", length(subject_indices) * spec$q)
  k <- 0L
  for (draw_pos in seq_along(subject_indices)) {
    src <- subject_indices[[draw_pos]]
    subj_meta <- subject_data[src, spec$group_vars, drop = FALSE]
    if (!ncol(subj_meta)) {
      subj_meta <- data.frame(.dummy = 1L)[, FALSE, drop = FALSE]
    }
    for (cell_idx in seq_len(spec$q)) {
      k <- k + 1L
      cell_meta <- spec$cell_data[cell_idx, spec$cell_vars, drop = FALSE]
      w <- if (is.null(spec$weight_flat)) {
        1
      } else {
        spec$weight_flat[(src - 1L) * spec$q + cell_idx]
      }
      source_rows[[k]] <- cbind(
        data.frame(subject_id = spec$subject_ids[[src]],
                   draw_index = draw_pos,
                   source_index = src,
                   cell_index = cell_idx,
                   weight = w,
                   stringsAsFactors = FALSE),
        cbind(subj_meta, cell_meta, stringsAsFactors = FALSE),
        stringsAsFactors = FALSE
      )
    }
  }
  map <- do.call(rbind, source_rows)
  row_meta <- unique(map[row_vars])
  row_meta <- row_meta[do.call(order, row_meta), , drop = FALSE]
  rownames(row_meta) <- NULL
  row_ids <- do.call(paste, c(row_meta[row_vars], sep = spec$row_sep))

  Y <- matrix(NA_real_, nrow = nrow(row_meta), ncol = spec$p,
              dimnames = list(row_ids, spec$feature_ids))
  map$aggregate_row <- NA_integer_
  for (i in seq_len(nrow(row_meta))) {
    sel <- rep(TRUE, nrow(map))
    for (rv in row_vars) {
      sel <- sel & as.character(map[[rv]]) == as.character(row_meta[[rv]][[i]])
    }
    idx <- which(sel)
    source_mat <- t(vapply(idx,
                           function(j) src_value(map$source_index[[j]],
                                                 map$cell_index[[j]]),
                           numeric(spec$p)))
    w <- if (identical(spec$aggregate, "weighted_mean")) {
      map$weight[idx]
    } else {
      rep(1, length(idx))
    }
    Y[i, ] <- colSums(source_mat * w) / sum(w)
    map$aggregate_row[idx] <- i
  }
  map <- map[, c("subject_id", "draw_index", "source_index", "cell_index",
                 "aggregate_row", "weight", row_vars), drop = FALSE]
  rownames(map) <- NULL
  list(Y = Y, row_ids = row_ids, row_metadata = row_meta, subject_to_row = map)
}

# Two-group x task x measure fixture with a group-modulated cell contrast.
make_aggregate_fixture <- function(n_per_group = 5, features = 3, sd = 0.1,
                                   group_effect = 3, use_cell_id = FALSE) {
  subject_data <- data.frame(
    subject_id = paste0("s", seq_len(2 * n_per_group)),
    group = rep(c("ctl", "pat"), each = n_per_group),
    stringsAsFactors = FALSE
  )
  cell_data <- data.frame(
    task = rep(c("recog", "nback"), each = 2),
    measure = rep(c("low", "high"), 2),
    stringsAsFactors = FALSE
  )
  cell_ids <- paste(cell_data$task, cell_data$measure, sep = "_")
  feature_ids <- paste0("f", seq_len(features))
  cell_contrast <- c(1, -1, -1, 1)
  feature_pattern <- seq_len(features) - (features + 1) / 2 + 0.5
  values <- lapply(seq_len(nrow(subject_data)), function(i) {
    group_sign <- if (subject_data$group[[i]] == "ctl") 1 else -1
    M <- group_sign * group_effect * cell_contrast %*% t(feature_pattern) +
      matrix(stats::rnorm(4 * features, sd = sd), 4, features)
    dimnames(M) <- list(cell_ids, feature_ids)
    M
  })
  names(values) <- subject_data$subject_id
  cell_id_col <- NULL
  if (isTRUE(use_cell_id)) {
    cell_data$cell <- cell_ids
    cell_id_col <- "cell"
  }
  target <- dkge_aggregate_target(
    values,
    subject_data = subject_data,
    cell_data = cell_data,
    group_vars = "group",
    cell_vars = c("task", "measure"),
    cell_id_col = cell_id_col
  )
  K <- diag(nrow(target$Y))
  dimnames(K) <- list(rownames(target$Y), rownames(target$Y))
  within <- with(target$row_metadata,
                 ifelse(task == "recog" & measure == "low", 1,
                        ifelse(task == "recog" & measure == "high", -1,
                               ifelse(task == "nback" & measure == "low", -1, 1))))
  row_contrast <- within * ifelse(target$row_metadata$group == "ctl", 1, -1)
  names(row_contrast) <- rownames(target$Y)
  list(target = target, K = K, row_contrast = row_contrast,
       subject_data = subject_data, cell_data = cell_data, values = values)
}

# Independent unaligned oracle for a single subject-label permutation.
# Rebuilds the permuted target the same way the worker does, then takes the
# statistic on the raw refit (no Procrustes).
unaligned_perm_stat <- function(target, K, idx, statistic, rank = 3L,
                                contrast = NULL, component = 1L,
                                group_vars = NULL) {
  spec <- target$resample_spec
  group_vars <- group_vars %||% target$group_vars
  subject_data_b <- spec$subject_data
  subject_data_b[group_vars] <- subject_data_b[idx, group_vars, drop = FALSE]
  target_b <- dkge:::.dkge_aggregate_from_spec(spec, subject_data = subject_data_b)
  fit_b <- dkge_aggregate_fit(target_b, K = K, rank = rank, center = "none")
  args <- list(fit = fit_b, statistic = statistic, component = component)
  if (!is.null(contrast)) args$contrast <- contrast
  do.call(dkge_aggregate_stat, args)
}


test_that("aggregate target builds reusable group-by-cell means", {
  values <- list(
    s1 = matrix(c(1, 2, 3, 4,
                  2, 3, 4, 5), 4, 2,
                dimnames = list(NULL, c("v1", "v2"))),
    s2 = matrix(c(2, 3, 4, 5,
                  3, 4, 5, 6), 4, 2,
                dimnames = list(NULL, c("v1", "v2"))),
    s3 = matrix(c(10, 20, 30, 40,
                  11, 21, 31, 41), 4, 2,
                dimnames = list(NULL, c("v1", "v2"))),
    s4 = matrix(c(12, 22, 32, 42,
                  13, 23, 33, 43), 4, 2,
                dimnames = list(NULL, c("v1", "v2")))
  )
  subject_data <- data.frame(
    subject_id = paste0("s", 1:4),
    group = rep(c("ctl", "pat"), each = 2),
    stringsAsFactors = FALSE
  )
  cell_data <- data.frame(
    task = rep(c("recog", "nback"), each = 2),
    measure = rep(c("low", "high"), 2),
    stringsAsFactors = FALSE
  )

  target <- dkge_aggregate_target(
    values,
    subject_data = subject_data,
    cell_data = cell_data,
    group_vars = "group",
    cell_vars = c("task", "measure")
  )

  expect_s3_class(target, "dkge_aggregate_target")
  expect_equal(dim(target$Y), c(8L, 2L))
  expect_equal(
    rownames(target$Y),
    c("ctl:nback:high", "ctl:nback:low", "ctl:recog:high", "ctl:recog:low",
      "pat:nback:high", "pat:nback:low", "pat:recog:high", "pat:recog:low")
  )
  expect_equal(
    unname(target$Y["ctl:recog:low", ]),
    unname(colMeans(rbind(values$s1[1, ], values$s2[1, ])))
  )
  expect_equal(nrow(target$subject_to_row), length(values) * nrow(cell_data))
  expect_equal(target$provenance$inference_unit, "subject")

  wtarget <- dkge_aggregate_target(
    values,
    subject_data = subject_data,
    cell_data = cell_data,
    group_vars = "group",
    cell_vars = c("task", "measure"),
    weights = c(3, 1, 1, 1),
    aggregate = "weighted_mean"
  )
  expect_equal(
    unname(wtarget$Y["ctl:recog:low", ]),
    unname((3 * values$s1[1, ] + values$s2[1, ]) / 4)
  )

  bad_subjects <- subject_data
  bad_subjects$subject_id[[1]] <- "missing"
  expect_error(
    dkge_aggregate_target(values, bad_subjects, cell_data,
                          group_vars = "group",
                          cell_vars = c("task", "measure")),
    "matching values"
  )
})

test_that("aggregate target rejects missing values in row-defining variables", {
  values <- replicate(4, matrix(seq_len(8), 4, 2), simplify = FALSE)
  names(values) <- paste0("s", 1:4)
  subject_data <- data.frame(
    subject_id = paste0("s", 1:4),
    group = c("ctl", NA, "pat", "pat"),
    stringsAsFactors = FALSE
  )
  cell_data <- data.frame(
    task = rep(c("recog", "nback"), each = 2),
    measure = rep(c("low", "high"), 2),
    stringsAsFactors = FALSE
  )

  expect_error(
    dkge_aggregate_target(values, subject_data, cell_data,
                          group_vars = "group",
                          cell_vars = c("task", "measure")),
    "Subject grouping variable 'group' contains missing values"
  )

  subject_data$group <- rep(c("ctl", "pat"), each = 2)
  bad_cells <- cell_data
  bad_cells$measure[[3]] <- NA
  expect_error(
    dkge_aggregate_target(values, subject_data, bad_cells,
                          group_vars = "group",
                          cell_vars = c("task", "measure")),
    "Cell variable 'measure' contains missing values"
  )
})

test_that("aggregate target validates and reorders cell rows by cell IDs", {
  cell_data <- data.frame(
    cell = c("recog_low", "recog_high", "nback_low", "nback_high"),
    task = rep(c("recog", "nback"), each = 2),
    measure = rep(c("low", "high"), 2),
    stringsAsFactors = FALSE
  )
  subject_data <- data.frame(
    subject_id = c("s1", "s2"),
    group = c("ctl", "ctl"),
    stringsAsFactors = FALSE
  )
  base <- matrix(c(1, 10, 2, 20,
                   3, 30, 4, 40), 4, 2,
                 dimnames = list(cell_data$cell, c("v1", "v2")))
  perm <- c("nback_high", "recog_low", "nback_low", "recog_high")
  values <- list(
    s1 = base[perm, , drop = FALSE],
    s2 = 2 * base[perm, , drop = FALSE]
  )

  target <- dkge_aggregate_target(
    values,
    subject_data = subject_data,
    cell_data = cell_data,
    group_vars = "group",
    cell_vars = c("task", "measure"),
    cell_id_col = "cell"
  )

  expect_equal(unname(target$Y["ctl:recog:low", ]),
               unname(colMeans(rbind(base["recog_low", ],
                                      2 * base["recog_low", ]))))

  bad_values <- values
  rownames(bad_values$s1)[[1]] <- "unknown"
  expect_error(
    dkge_aggregate_target(bad_values, subject_data, cell_data,
                          group_vars = "group",
                          cell_vars = c("task", "measure"),
                          cell_id_col = "cell"),
    "contain all"
  )

  dup_cell <- cell_data
  dup_cell$cell[[2]] <- dup_cell$cell[[1]]
  expect_error(
    dkge_aggregate_target(values, subject_data, dup_cell,
                          group_vars = "group",
                          cell_vars = c("task", "measure"),
                          cell_id_col = "cell"),
    "unique non-empty"
  )
})

test_that("aggregate target matches cells and features by name", {
  A <- matrix(c(1, 2, 3, 4), 2, 2,
              dimnames = list(c("c1", "c2"), c("f1", "f2")))
  subject_data <- data.frame(
    subject_id = c("s1", "s2"),
    g = c("g", "g"),
    stringsAsFactors = FALSE
  )

  target <- dkge_aggregate_target(
    list(s1 = A, s2 = A),
    subject_data,
    cell_data = NULL,
    group_vars = "g"
  )
  expect_equal(target$row_ids, c("g:c1", "g:c2"))
  expect_equal(unname(target$Y), unname(A), tolerance = 1e-12)

  # Collapse is opt-in: cell_vars = character(0) drops the inferred cell.
  collapsed <- dkge_aggregate_target(
    list(s1 = A, s2 = A),
    subject_data,
    cell_data = NULL,
    group_vars = "g",
    cell_vars = character(0)
  )
  expect_equal(collapsed$row_ids, "g")
  expect_equal(as.numeric(collapsed$Y), unname(colMeans(A)), tolerance = 1e-12)

  row_perm <- dkge_aggregate_target(
    list(s1 = A, s2 = A[2:1, , drop = FALSE]),
    subject_data,
    cell_data = NULL,
    group_vars = "g"
  )
  expect_equal(unname(row_perm$Y), unname(A), tolerance = 1e-12)

  s1 <- matrix(c(1, 10, 3, 20, 5, 30), 2, 3,
               dimnames = list(c("c1", "c2"), c("f1", "f2", "f3")))
  s2 <- matrix(c(2, 12, 4, 24, 6, 36), 2, 3,
               dimnames = list(c("c1", "c2"), c("f1", "f2", "f3")))
  s2 <- s2[, c("f3", "f1", "f2"), drop = FALSE]
  col_perm <- dkge_aggregate_target(
    list(s1 = s1, s2 = s2),
    subject_data,
    cell_data = NULL,
    group_vars = "g"
  )
  expect_equal(unname(col_perm$Y["g:c1", ]), c(1.5, 3.5, 5.5),
               tolerance = 1e-12)

  bad_cols <- s2
  colnames(bad_cols) <- c("f3", "f1", "fx")
  expect_error(
    dkge_aggregate_target(list(s1 = s1, s2 = bad_cols), subject_data,
                          cell_data = NULL, group_vars = "g"),
    "feature columns of subject `s2` do not match subject 1"
  )

  partial <- A
  rownames(partial) <- NULL
  expect_error(
    dkge_aggregate_target(list(s1 = A, s2 = partial), subject_data,
                          cell_data = NULL, group_vars = "g"),
    "every subject or for none"
  )

  fx_id <- make_aggregate_fixture(n_per_group = 2, features = 2, sd = 0,
                                  use_cell_id = TRUE)
  expect_equal(fx_id$target$resample_spec$cell_id_col, "cell")
  expect_equal(rownames(fx_id$values[[1]]),
               c("recog_low", "recog_high", "nback_low", "nback_high"))
})

test_that("vectorised aggregation reproduces the row-scan implementation", {
  set.seed(4)
  n <- 9L
  q <- 6L
  p <- 5L
  subject_data <- data.frame(
    subject_id = paste0("s", seq_len(n)),
    # unequal cell counts: 4 / 3 / 2 subjects, crossed with an unbalanced site
    group = c(rep("a", 4), rep("b", 3), rep("c", 2)),
    site = c("x", "y", "x", "y", "x", "x", "y", "y", "x"),
    stringsAsFactors = FALSE
  )
  cell_data <- data.frame(
    task = rep(c("t1", "t2", "t3"), each = 2),
    measure = rep(c("lo", "hi"), 3),
    stringsAsFactors = FALSE
  )
  values <- lapply(seq_len(n), function(i) matrix(rnorm(q * p), q, p))
  names(values) <- subject_data$subject_id
  W <- matrix(runif(n * q, 0.1, 3), n, q,
              dimnames = list(subject_data$subject_id, NULL))

  for (agg in c("mean", "weighted_mean")) {
    target <- dkge_aggregate_target(values, subject_data, cell_data,
                                    group_vars = c("group", "site"),
                                    cell_vars = c("task", "measure"),
                                    weights = W, aggregate = agg)
    spec <- target$resample_spec
    ref <- old_aggregate_from_spec(spec)
    expect_equal(target$Y, ref$Y, tolerance = 1e-12)
    expect_identical(target$row_ids, ref$row_ids)
    expect_equal(target$row_metadata, ref$row_metadata)
    map <- target$subject_to_row
    rownames(map) <- NULL
    expect_equal(map, ref$subject_to_row)

    # bootstrap-style resample (duplicated subjects)
    idx <- c(1L, 1L, 3L, 4L, 4L, 4L, 7L, 8L, 9L)
    boot_new <- dkge:::.dkge_aggregate_from_spec(spec, subject_indices = idx)
    boot_old <- old_aggregate_from_spec(spec, subject_indices = idx)
    expect_equal(boot_new$Y, boot_old$Y, tolerance = 1e-12)
    map <- boot_new$subject_to_row
    rownames(map) <- NULL
    expect_equal(map, boot_old$subject_to_row)

    # permuted subject metadata
    perm_data <- spec$subject_data
    pidx <- c(5L, 2L, 9L, 1L, 7L, 3L, 8L, 4L, 6L)
    perm_data[c("group", "site")] <- perm_data[pidx, c("group", "site"),
                                               drop = FALSE]
    perm_new <- dkge:::.dkge_aggregate_from_spec(spec, subject_data = perm_data)
    perm_old <- old_aggregate_from_spec(spec, subject_data = perm_data)
    expect_equal(perm_new$Y, perm_old$Y, tolerance = 1e-12)
    expect_identical(perm_new$row_ids, perm_old$row_ids)
  }

  # factor-valued row variables must keep their level ordering
  sd_fac <- subject_data
  sd_fac$group <- factor(sd_fac$group, levels = c("c", "b", "a"))
  cd_fac <- cell_data
  cd_fac$task <- factor(cd_fac$task, levels = c("t3", "t2", "t1"))
  target <- dkge_aggregate_target(values, sd_fac, cd_fac, group_vars = "group",
                                  cell_vars = c("task", "measure"))
  ref <- old_aggregate_from_spec(target$resample_spec)
  expect_equal(target$Y, ref$Y, tolerance = 1e-12)
  expect_identical(target$row_ids, ref$row_ids)
  expect_identical(target$row_ids[[1]], "c:t3:hi")
})

test_that("bootstrap draws record the resampled subjects with duplicates", {
  fx <- make_aggregate_fixture(n_per_group = 2, features = 2, sd = 0)
  spec <- fx$target$resample_spec
  expect_identical(fx$target$subject_ids, spec$subject_ids)

  idx <- c(1L, 1L, 3L, 4L)
  drawn <- dkge:::.dkge_aggregate_from_spec(spec, subject_indices = idx)
  expect_identical(drawn$subject_ids, spec$subject_ids[idx])
  expect_identical(drawn$subject_indices, idx)
  expect_identical(drawn$source_subject_ids, spec$subject_ids)
  expect_equal(table(drawn$subject_to_row$subject_id)[["s1"]], 2L * spec$q)
})

test_that("aggregate fit reproduces the SVD of K^{1/2} Y", {
  set.seed(31)
  row_ids <- paste0("r", 1:6)
  Y <- matrix(rnorm(6 * 4), 6, 4,
              dimnames = list(row_ids, paste0("v", 1:4)))
  A <- matrix(rnorm(36), 6, 6)
  K <- crossprod(A) + diag(6)
  dimnames(K) <- list(row_ids, row_ids)

  fit <- dkge_aggregate_fit(Y, K = K, rank = 4)
  ek <- eigen((K + t(K)) / 2, symmetric = TRUE)
  Khalf <- ek$vectors %*% diag(sqrt(pmax(ek$values, 0))) %*% t(ek$vectors)
  reference <- svd(Khalf %*% Y)

  expect_equal(unname(fit$singular_values), reference$d[1:4], tolerance = 1e-10)
  expect_equal(unname(sqrt(fit$eig_values[1:4])), reference$d[1:4],
               tolerance = 1e-10)
  # K-orthonormal basis, and salience scores carry the same energy
  expect_equal(unname(t(fit$U) %*% K %*% fit$U), diag(4), tolerance = 1e-10)
  expect_equal(unname(sqrt(colSums(fit$scores_feature^2))),
               reference$d[1:4], tolerance = 1e-10)
})

test_that("aggregate fit matches kernels on whichever dimnames are present", {
  set.seed(32)
  row_ids <- paste0("r", 1:5)
  Y <- matrix(rnorm(15), 5, 3, dimnames = list(row_ids, paste0("v", 1:3)))
  A <- matrix(rnorm(25), 5, 5)
  K <- crossprod(A) + diag(5)
  dimnames(K) <- list(row_ids, row_ids)
  reference <- dkge_aggregate_fit(Y, K = K, rank = 2)

  ord <- c(3L, 1L, 5L, 2L, 4L)
  K_rows_only <- K[ord, ord]
  rownames(K_rows_only) <- row_ids[ord]
  colnames(K_rows_only) <- NULL
  expect_equal(dkge_aggregate_fit(Y, K = K_rows_only, rank = 2)$singular_values,
               reference$singular_values, tolerance = 1e-10)

  K_cols_only <- K[ord, ord]
  rownames(K_cols_only) <- NULL
  colnames(K_cols_only) <- row_ids[ord]
  expect_equal(dkge_aggregate_fit(Y, K = K_cols_only, rank = 2)$singular_values,
               reference$singular_values, tolerance = 1e-10)

  bad <- K
  rownames(bad) <- paste0("z", 1:5)
  colnames(bad) <- NULL
  expect_error(dkge_aggregate_fit(Y, K = bad, rank = 2),
               "Kernel dimnames must match aggregate row IDs")

  dup <- K
  rownames(dup) <- c("r1", "r1", "r3", "r4", "r5")
  colnames(dup) <- NULL
  expect_error(dkge_aggregate_fit(Y, K = dup, rank = 2),
               "unique and non-empty")
})

test_that("aggregate statistics are signed and match a direct SVD", {
  row_ids <- c("A:t1:c1", "A:t1:c2", "A:t2:c1", "A:t2:c2",
               "B:t1:c1", "B:t1:c2", "B:t2:c1", "B:t2:c2")
  contrast <- c(1, -1, -1, 1, -1, 1, 1, -1)
  names(contrast) <- row_ids
  Y <- contrast %*% matrix(c(3, -2, 1), 1, 3) +
    matrix(seq_along(row_ids) / 100, length(row_ids), 3)
  rownames(Y) <- row_ids
  colnames(Y) <- paste0("v", 1:3)
  K <- diag(length(row_ids))
  dimnames(K) <- list(row_ids, row_ids)

  fit <- dkge_aggregate_fit(Y, K = K, rank = 2)
  stat <- dkge_aggregate_stat(fit, "contrast_score", contrast = contrast)
  salience_stat <- dkge_aggregate_stat(fit, "salience_contrast",
                                       contrast = contrast, component = 1)
  score_stat <- dkge_aggregate_stat(fit, "component_score_contrast",
                                    contrast = contrast, component = 1)
  bridge_stat <- dkge_aggregate_stat(fit, "between_group_contrast",
                                     contrast = contrast, component = 1)

  # With K = I the saliences are the left singular vectors of Y, so the
  # statistics have closed forms up to the arbitrary component sign.
  reference <- svd(Y)
  expect_equal(abs(salience_stat), abs(sum(contrast * reference$u[, 1])),
               tolerance = 1e-10)
  expect_equal(abs(score_stat),
               abs(sum(contrast * reference$u[, 1]) * reference$d[[1]]),
               tolerance = 1e-10)
  expect_equal(stat, salience_stat)
  expect_equal(bridge_stat, score_stat)
  expect_equal(unname(fit$singular_values), reference$d[1:2], tolerance = 1e-10)

  # Signed, not folded: negating the contrast negates the statistic.
  expect_equal(
    dkge_aggregate_stat(fit, "salience_contrast", contrast = -contrast),
    -salience_stat
  )
  expect_equal(
    dkge_aggregate_stat(fit, "component_score_contrast", contrast = -contrast),
    -score_stat
  )
  expect_true(abs(salience_stat) > 1e-8)

  lv2_stat <- dkge_aggregate_stat(fit, "component_score_contrast",
                                  contrast = contrast, component = 2)
  expect_equal(abs(lv2_stat),
               abs(sum(contrast * reference$u[, 2]) * reference$d[[2]]),
               tolerance = 1e-10)
  expect_false(isTRUE(all.equal(score_stat, lv2_stat)))

  # Row order is not identified up to component sign, so compare magnitudes.
  ord <- c(5:8, 1:4)
  fit_reordered <- dkge_aggregate_fit(Y[ord, ], K = K[ord, ord], rank = 2)
  stat_reordered <- dkge_aggregate_stat(fit_reordered, "contrast_score",
                                        contrast = contrast)
  expect_equal(abs(stat), abs(stat_reordered), tolerance = 1e-8)

  expect_equal(
    dkge_aggregate_stat(fit, function(x) x$singular_values[[1]]),
    unname(fit$singular_values[[1]])
  )
  expect_equal(fit$estimand, "aggregate_cell_mean")
  expect_error(dkge_aggregate_stat(fit, "salience_contrast"),
               "`contrast` is required")
})

test_that("aggregate alignment handles signs, swaps, ties, and truncation", {
  row_ids <- paste0("r", 1:4)
  K <- diag(4)
  dimnames(K) <- list(row_ids, row_ids)
  Y <- diag(4)
  rownames(Y) <- row_ids
  fit <- dkge_aggregate_fit(Y, K = K, rank = 2)

  flipped <- fit
  flipped$U <- -flipped$U
  flipped$saliences <- -flipped$saliences
  flipped$scores_feature <- -flipped$scores_feature
  aligned1 <- dkge_aggregate_align(fit, flipped, rank = 1)
  expect_equal(aligned1$alignment$method, "sign")
  expect_equal(aligned1$U, fit$U[, 1, drop = FALSE], tolerance = 1e-8)

  swapped <- fit
  swapped$U <- swapped$U[, 2:1, drop = FALSE]
  swapped$saliences <- swapped$saliences[, 2:1, drop = FALSE]
  swapped$scores_feature <- swapped$scores_feature[, 2:1, drop = FALSE]
  aligned2 <- dkge_aggregate_align(fit, swapped)
  expect_equal(aligned2$alignment$method, "K-procrustes")
  expect_equal(unname(t(fit$U) %*% K %*% aligned2$U),
               diag(2), tolerance = 1e-8)
  expect_true(aligned2$alignment$near_tie)

  # Singular values must travel with the rotated components.
  expect_equal(unname(aligned2$singular_values),
               unname(sqrt(colSums(aligned2$scores_feature^2))))
  expect_identical(colnames(aligned2$U), colnames(fit$U))
  expect_identical(names(aligned2$singular_values), colnames(fit$U))

  truncated <- dkge_aggregate_align(fit, swapped, rank = 1)
  expect_equal(truncated$rank, 1L)
  expect_equal(ncol(truncated$U), 1L)
  expect_equal(unname(truncated$singular_values),
               unname(sqrt(colSums(truncated$scores_feature^2))))
})

test_that("alignment rotates singular values with the components", {
  set.seed(77)
  row_ids <- paste0("r", 1:5)
  Y <- matrix(rnorm(5 * 6), 5, 6,
              dimnames = list(row_ids, paste0("v", 1:6)))
  K <- diag(5)
  dimnames(K) <- list(row_ids, row_ids)
  fit <- dkge_aggregate_fit(Y, K = K, rank = 3)
  # A genuinely non-degenerate spectrum, so a mis-sliced singular value is
  # numerically obvious.
  expect_gt(min(abs(diff(fit$singular_values))), 0.1)

  swap_component <- function(fit, perm) {
    out <- fit
    out$U <- fit$U[, perm, drop = FALSE]
    colnames(out$U) <- colnames(fit$U)
    out$saliences <- fit$K %*% out$U
    dimnames(out$saliences) <- dimnames(out$U)
    out$scores_feature <- t(fit$Y) %*% out$saliences
    colnames(out$scores_feature) <- colnames(out$U)
    out$singular_values <- sqrt(colSums(out$scores_feature^2))
    out
  }
  scrambled <- swap_component(fit, c(3L, 1L, 2L))
  expect_equal(unname(scrambled$singular_values),
               unname(fit$singular_values[c(3, 1, 2)]), tolerance = 1e-10)

  aligned <- dkge_aggregate_align(fit, scrambled)
  expect_equal(unname(aligned$U), unname(fit$U), tolerance = 1e-8)
  expect_equal(unname(aligned$singular_values), unname(fit$singular_values),
               tolerance = 1e-8)
  expect_equal(unname(aligned$singular_values),
               unname(sqrt(colSums(aligned$scores_feature^2))),
               tolerance = 1e-12)

  contrast <- stats::setNames(c(1, -1, 2, 0, -2), row_ids)
  for (k in 1:3) {
    expect_equal(
      dkge_aggregate_stat(aligned, "component_score_contrast",
                          contrast = contrast, component = k),
      dkge_aggregate_stat(fit, "component_score_contrast",
                          contrast = contrast, component = k),
      tolerance = 1e-8
    )
  }

  # `Chat`/`eig_values` are rotation invariant and are documented as such.
  expect_identical(aligned$Chat, fit$Chat)
  expect_identical(aligned$eig_values, fit$eig_values)
})

test_that("aggregate alignment requires matching row identities", {
  row_ids <- paste0("r", 1:4)
  K <- diag(4)
  dimnames(K) <- list(row_ids, row_ids)
  Y <- diag(4)
  rownames(Y) <- row_ids
  fit <- dkge_aggregate_fit(Y, K = K, rank = 2)

  Y2 <- Y
  rownames(Y2) <- rev(row_ids)
  K2 <- K
  dimnames(K2) <- list(rev(row_ids), rev(row_ids))
  other <- dkge_aggregate_fit(Y2, K = K2, rank = 2)

  expect_error(dkge_aggregate_align(fit, other),
               "same aggregate row IDs")
})

test_that("aggregate permutation and bootstrap resample subjects and refit", {
  set.seed(11)
  fx <- make_aggregate_fixture(n_per_group = 5, features = 3, sd = 0.1)
  target <- fx$target

  perm <- dkge_aggregate_permute(
    target,
    K = fx$K,
    statistic = "contrast_score",
    contrast = fx$row_contrast,
    B = 39,
    rank = 1,
    seed = 101
  )
  expect_s3_class(perm, "dkge_aggregate_permutation")
  expect_equal(length(perm$null), 39L)
  expect_equal(perm$resampling$inference_unit, "subject")
  expect_equal(perm$alternative, "two.sided")
  expect_named(perm$alignment_summary,
               c("min_cosine", "median_cosine", "n_weak", "n_near_tie",
                 "weak_cosine", "weak_alignment"))
  expect_lte(perm$p, 0.2)
  expect_equal(perm$p,
               (1 + sum(abs(perm$null) >= abs(perm$observed) -
                          sqrt(.Machine$double.eps) * max(1, abs(perm$observed)))) / 40)

  # `contrast_score` on the observed fit has a closed form (K = I).
  reference <- svd(target$Y)
  expect_equal(abs(perm$observed),
               abs(sum(fx$row_contrast * reference$u[, 1])), tolerance = 1e-8)

  boot <- dkge_aggregate_bootstrap(
    target,
    K = fx$K,
    statistic = "contrast_score",
    contrast = fx$row_contrast,
    component_contrasts = matrix(fx$row_contrast, ncol = 1,
                                 dimnames = list(names(fx$row_contrast), "gxint")),
    return_features = TRUE,
    B = 12,
    rank = 1,
    seed = 102
  )
  expect_s3_class(boot, "dkge_aggregate_bootstrap")
  expect_equal(length(boot$statistics), 12L)
  expect_true(all(is.finite(boot$statistics)))
  expect_equal(boot$resampling$inference_unit, "subject")
  expect_named(boot$alignment_summary,
               c("min_cosine", "median_cosine", "n_weak", "n_near_tie",
                 "weak_cosine", "weak_alignment"))
  expect_equal(dim(boot$feature_bootstrap$z), c(3L, 1L))
  expect_equal(boot$feature_bootstrap$n, 12L)
  expect_equal(nrow(boot$component_contrasts$summary), 1L)
  expect_true(all(is.finite(boot$component_contrasts$summary$observed)))
  expect_true(boot$interval[[1]] <= boot$interval[[2]])
  expect_equal(boot$excludes_zero,
               boot$interval[[1]] > 0 || boot$interval[[2]] < 0)
  # Bootstrap ratios are observed / bootstrap SD.
  expect_equal(boot$feature_bootstrap$z,
               boot$feature_bootstrap$observed / boot$feature_bootstrap$sd)
})

test_that("rank > 1 permutation and bootstrap track the reference components", {
  set.seed(13)
  fx <- make_aggregate_fixture(n_per_group = 6, features = 4, sd = 0.4)
  target <- fx$target
  reference <- svd(target$Y)

  perm <- dkge_aggregate_permute(
    target,
    K = fx$K,
    statistic = "component_score_contrast",
    contrast = fx$row_contrast,
    component = 2,
    B = 25,
    rank = 3,
    seed = 303
  )
  expect_equal(perm$fit$rank, 3L)
  expect_equal(unname(perm$fit$singular_values), reference$d[1:3],
               tolerance = 1e-10)
  expect_equal(abs(perm$observed),
               abs(sum(fx$row_contrast * reference$u[, 2]) * reference$d[[2]]),
               tolerance = 1e-8)
  expect_true(all(is.finite(perm$null)))
  expect_true(perm$p > 0 && perm$p <= 1)
  expect_equal(length(perm$alignments[[1]]$cosines), 3L)

  boot <- dkge_aggregate_bootstrap(
    target,
    K = fx$K,
    statistic = "component_score_contrast",
    contrast = fx$row_contrast,
    component = 1,
    component_contrasts = matrix(fx$row_contrast, ncol = 1,
                                 dimnames = list(names(fx$row_contrast), "gxint")),
    return_features = TRUE,
    B = 20,
    rank = 3,
    seed = 304
  )
  expect_true(all(is.finite(boot$statistics)))
  expect_equal(dim(boot$component_contrasts$samples), c(1L, 3L, 20L))
  expect_equal(dim(boot$feature_bootstrap$z), c(4L, 3L))
  expect_equal(abs(boot$observed),
               abs(sum(fx$row_contrast * reference$u[, 1]) * reference$d[[1]]),
               tolerance = 1e-8)
  expect_true(!is.null(boot$alignment_summary$min_cosine))
})

test_that("stratified bootstrap resamples within strata", {
  strata_index <- split(1:4, factor(c("b", "b", "a", "b")))
  set.seed(5)
  draws <- replicate(200, dkge:::.dkge_aggregate_bootstrap_indices(strata_index),
                     simplify = FALSE)
  # Stratum "a" holds only subject 3; every draw must return exactly that
  # subject (the historical `sample(ii, 1)` would have sampled from 1:3).
  expect_true(all(vapply(draws, function(x) identical(x[[1]], 3L), logical(1))))
  expect_true(all(vapply(draws, function(x) all(x[-1] %in% c(1L, 2L, 4L)),
                         logical(1))))
  expect_true(all(vapply(draws, length, integer(1)) == 4L))
  # ... and the larger stratum really does vary.
  expect_gt(length(unique(vapply(draws, function(x) paste(x[-1], collapse = "-"),
                                 character(1)))), 1L)

  empty <- split(integer(0), factor(character(0), levels = "z"))
  expect_identical(dkge:::.dkge_aggregate_bootstrap_indices(empty), integer(0))
})

test_that("singleton strata reproduce the observed fit and are warned about", {
  set.seed(17)
  fx <- make_aggregate_fixture(n_per_group = 3, features = 3, sd = 0.2)
  target <- fx$target

  expect_warning(
    boot <- dkge_aggregate_bootstrap(
      target,
      K = fx$K,
      statistic = "component_score_contrast",
      contrast = fx$row_contrast,
      strata = target$resample_spec$subject_ids,
      B = 8,
      rank = 2,
      seed = 501
    ),
    "fewer than 2 subjects"
  )
  # Every stratum is a singleton, so each draw is the observed sample.
  expect_equal(boot$statistics, rep(boot$observed, 8L), tolerance = 1e-10)
  expect_equal(boot$interval, rep(boot$observed, 2L), tolerance = 1e-10)

  expect_error(
    dkge_aggregate_bootstrap(target, K = fx$K, statistic = "singular_value",
                             strata = c("a", NA, "a", "b", "b", "b"),
                             B = 4, rank = 1),
    "must not contain missing values"
  )
  expect_error(
    dkge_aggregate_bootstrap(target, K = fx$K, statistic = "singular_value",
                             strata = c("a", "a", "b"), B = 4, rank = 1),
    "must match the number of source subjects"
  )
})

test_that("bootstrap interval on a null contrast can contain zero", {
  set.seed(23)
  subject_data <- data.frame(
    subject_id = paste0("s", 1:12),
    group = rep(c("ctl", "pat"), each = 6),
    stringsAsFactors = FALSE
  )
  cell_data <- data.frame(
    task = rep(c("recog", "nback"), each = 2),
    measure = rep(c("low", "high"), 2),
    stringsAsFactors = FALSE
  )
  cell_contrast <- c(1, -1, -1, 1)
  feature_pattern <- c(1.5, -0.75, 0.5)
  # Shared cell structure, no group difference: the group-by-cell contrast is
  # null by construction.
  values <- lapply(seq_len(nrow(subject_data)), function(i) {
    cell_contrast %*% t(feature_pattern) +
      matrix(rnorm(4 * 3, sd = 0.25), 4, 3)
  })
  names(values) <- subject_data$subject_id

  target <- dkge_aggregate_target(values, subject_data, cell_data,
                                  group_vars = "group",
                                  cell_vars = c("task", "measure"))
  K <- diag(nrow(target$Y))
  dimnames(K) <- list(rownames(target$Y), rownames(target$Y))
  within <- with(target$row_metadata,
                 ifelse(task == "recog" & measure == "low", 1,
                        ifelse(task == "recog" & measure == "high", -1,
                               ifelse(task == "nback" & measure == "low", -1, 1))))
  row_contrast <- within * ifelse(target$row_metadata$group == "ctl", 1, -1)
  names(row_contrast) <- rownames(target$Y)

  boot <- dkge_aggregate_bootstrap(
    target,
    K = K,
    statistic = "between_group_contrast",
    contrast = row_contrast,
    B = 200,
    rank = 1,
    seed = 606
  )

  expect_true(any(boot$statistics < 0) && any(boot$statistics > 0))
  expect_lt(boot$interval[[1]], 0)
  expect_gt(boot$interval[[2]], 0)
  expect_false(boot$excludes_zero)
  expect_gte(boot$observed, boot$interval[[1]])
  expect_lte(boot$observed, boot$interval[[2]])
})

test_that("aggregate permutation is calibrated at smoke scale under a null", {
  set.seed(21)
  subject_data <- data.frame(
    subject_id = paste0("s", 1:8),
    group = rep(c("ctl", "pat"), each = 4),
    stringsAsFactors = FALSE
  )
  cell_data <- data.frame(
    task = rep(c("recog", "nback"), each = 2),
    measure = rep(c("low", "high"), 2),
    stringsAsFactors = FALSE
  )
  cell_contrast <- c(1, -1, -1, 1)
  feature_pattern <- c(1.5, -0.75)
  values <- lapply(seq_len(nrow(subject_data)), function(i) {
    subject_shift <- rnorm(1, sd = 0.02)
    cell_contrast %*% t(feature_pattern) +
      subject_shift +
      matrix(rnorm(4 * 2, sd = 0.03), 4, 2)
  })
  names(values) <- subject_data$subject_id

  target <- dkge_aggregate_target(
    values,
    subject_data = subject_data,
    cell_data = cell_data,
    group_vars = "group",
    cell_vars = c("task", "measure")
  )
  K <- diag(nrow(target$Y))
  dimnames(K) <- list(rownames(target$Y), rownames(target$Y))
  within <- with(target$row_metadata, ifelse(task == "recog" & measure == "low", 1,
                                             ifelse(task == "recog" & measure == "high", -1,
                                                    ifelse(task == "nback" & measure == "low", -1, 1))))
  row_contrast <- within * ifelse(target$row_metadata$group == "ctl", 1, -1)
  names(row_contrast) <- rownames(target$Y)

  perm <- dkge_aggregate_permute(
    target,
    K = K,
    statistic = "between_group_contrast",
    contrast = row_contrast,
    B = 39,
    rank = 1,
    seed = 202
  )

  expect_s3_class(perm, "dkge_aggregate_permutation")
  expect_gt(perm$p, 0.2)
  expect_true(all(is.finite(perm$null)))
  # A two-sided null must place mass on both sides of zero.
  expect_true(any(perm$null < 0) && any(perm$null > 0))

  one_sided <- dkge_aggregate_permute(
    target, K = K, statistic = "between_group_contrast",
    contrast = row_contrast, B = 39, rank = 1, alternative = "greater",
    seed = 202
  )
  expect_equal(one_sided$null, perm$null)
  expect_equal(one_sided$alternative, "greater")
  expect_equal(one_sided$p,
               (1 + sum(perm$null >= perm$observed -
                          sqrt(.Machine$double.eps) * max(1, abs(perm$observed)))) / 40)

  # Rank 2 is the first case where K-Procrustes alignment of the null would
  # shrink signed statistics; the smoke null must still look like a null.
  perm2 <- dkge_aggregate_permute(
    target, K = K, statistic = "between_group_contrast",
    contrast = row_contrast, B = 39, rank = 2, seed = 202
  )
  expect_gt(perm2$p, 0.2)
  expect_true(all(is.finite(perm2$null)))
  expect_true(any(perm2$null < 0) && any(perm2$null > 0))
})

test_that("bootstrap alignment stays on the reference of a strong-signal fit", {
  set.seed(19)
  fx <- make_aggregate_fixture(n_per_group = 6, features = 4, sd = 0.05,
                               group_effect = 5)
  # Rank 1 is the identified group-by-cell direction; higher components are
  # noise and their Procrustes cosines are not a test of alignment quality.
  boot <- dkge_aggregate_bootstrap(
    fx$target, K = fx$K, statistic = "singular_value",
    B = 20, rank = 1, seed = 304
  )
  expect_gt(boot$alignment_summary$min_cosine, 0.9)
})

test_that("permutation null statistic equals the unaligned refit", {
  fx <- null_aggregate_target(seed = 88, n = 16, q = 4, p = 8)
  target <- fx$target
  n_subj <- nrow(target$resample_spec$subject_data)
  # Replay the seeded index draw: set.seed then the same fit that permute
  # evaluates before sample.int().
  set.seed(505)
  invisible(dkge_aggregate_fit(target, K = fx$K, rank = 3, center = "none"))
  idx <- lapply(seq_len(5), function(b) sample.int(n_subj))

  stats <- c("singular_value", "salience_contrast", "component_score_contrast")
  for (stat in stats) {
    extra <- if (identical(stat, "singular_value")) {
      list()
    } else {
      list(contrast = fx$contrast)
    }
    pm <- do.call(dkge_aggregate_permute, c(
      list(target, K = fx$K, statistic = stat, B = 5, rank = 3, seed = 505),
      extra
    ))
    oracle <- vapply(idx, function(i) {
      unaligned_perm_stat(target, fx$K, i, stat, rank = 3,
                          contrast = fx$contrast)
    }, numeric(1))
    expect_equal(unname(pm$null), unname(oracle), tolerance = 1e-12,
                 info = stat)
    expect_true(!is.null(pm$alignment_summary))
  }

  # Alignment is diagnostic only: flipping a reference component must not
  # change the permutation statistic.
  observed_fit <- dkge_aggregate_fit(target, K = fx$K, rank = 3, center = "none")
  reference <- dkge:::.dkge_aggregate_alignment_reference(observed_fit)
  flipped <- reference
  flipped$U[, 1] <- -flipped$U[, 1]
  worker <- dkge:::.dkge_aggregate_permute_worker(
    spec = target$resample_spec, reference = reference,
    row_ids = target$row_ids, group_vars = "grp", center = "none",
    statistic = "component_score_contrast",
    stat_args = list(component = 1L, contrast = fx$contrast)
  )
  worker_flip <- dkge:::.dkge_aggregate_permute_worker(
    spec = target$resample_spec, reference = flipped,
    row_ids = target$row_ids, group_vars = "grp", center = "none",
    statistic = "component_score_contrast",
    stat_args = list(component = 1L, contrast = fx$contrast)
  )
  expect_equal(worker(idx[[1]])$stat, worker_flip(idx[[1]])$stat,
               tolerance = 1e-12)
})

test_that("permutation type-I at rank 3 is nominal", {
  skip_on_cran()
  stats <- c("singular_value", "salience_contrast", "component_score_contrast")
  seeds <- 4000L + seq_len(100L)
  for (stat in stats) {
    rate <- type1_rate(function(s) {
      fx <- null_aggregate_target(seed = s)
      extra <- if (identical(stat, "singular_value")) {
        list()
      } else {
        list(contrast = fx$contrast)
      }
      pm <- do.call(dkge_aggregate_permute, c(
        list(fx$target, K = fx$K, statistic = stat,
             B = 99, rank = 3, seed = s + 10000L),
        extra
      ))
      pm$p
    }, n_rep = 100L, alpha = 0.05, seeds = seeds)
    expect_true(
      rate >= 0.01 && rate <= 0.11,
      info = sprintf("%s type-I rate = %.3f (want [0.01, 0.11])", stat, rate)
    )
  }
})

test_that("permutation rejects group_vars that would change the row set", {
  set.seed(29)
  n <- 8L
  subject_data <- data.frame(
    subject_id = paste0("s", seq_len(n)),
    group = rep(c("ctl", "pat"), each = 4),
    site = rep(c("x", "y"), 4),
    stringsAsFactors = FALSE
  )
  cell_data <- data.frame(cell = c("c1", "c2"), stringsAsFactors = FALSE)
  values <- lapply(seq_len(n), function(i) matrix(rnorm(4), 2, 2))
  names(values) <- subject_data$subject_id
  target <- dkge_aggregate_target(values, subject_data, cell_data,
                                  group_vars = c("group", "site"),
                                  cell_vars = "cell")

  expect_error(
    dkge_aggregate_permute(target, statistic = "singular_value",
                           group_vars = "group", B = 5, rank = 1),
    "would stay fixed"
  )
  expect_error(
    dkge_aggregate_permute(target, statistic = "singular_value",
                           group_vars = c("group", "site", "nope"),
                           B = 5, rank = 1),
    "Unknown subject grouping variable"
  )
  # Permuting all row-defining variables jointly is fine.
  perm <- dkge_aggregate_permute(target, statistic = "singular_value",
                                 group_vars = c("group", "site"),
                                 B = 5, rank = 1, seed = 9)
  expect_equal(length(perm$null), 5L)
})

test_that("feature bootstrap ratios are NA when the bootstrap SD collapses", {
  observed <- matrix(c(2, -3, 0.5, 1), 2, 2,
                     dimnames = list(c("v1", "v2"), c("LV1", "LV2")))
  state <- list(n = 0L,
                mean = matrix(0, 2, 2, dimnames = dimnames(observed)),
                m2 = matrix(0, 2, 2, dimnames = dimnames(observed)))
  for (b in 1:5) {
    draw <- observed
    draw[1, 1] <- observed[1, 1] + b * 0.1
    state <- dkge:::.dkge_welford_update_matrix(state, draw)
  }
  out <- dkge:::.dkge_feature_bootstrap_finalize(state, observed = observed)

  expect_equal(out$sd[2, 1], 0)
  expect_true(is.na(out$z[2, 1]))
  expect_true(all(is.na(out$z[, 2])))
  expect_equal(out$z[1, 1], observed[1, 1] / out$sd[1, 1])

  # A single draw has no SD at all.
  single <- dkge:::.dkge_welford_update_matrix(
    list(n = 0L,
         mean = matrix(0, 2, 2, dimnames = dimnames(observed)),
         m2 = matrix(0, 2, 2, dimnames = dimnames(observed))),
    observed
  )
  out1 <- dkge:::.dkge_feature_bootstrap_finalize(single, observed = observed)
  expect_true(all(is.na(out1$z)))
})

test_that("parallel resampling reproduces the serial result", {
  skip_if_not_installed("future.apply")
  # future re-signals namespace build-version warnings from its workers; mute
  # only those, so any real warning still surfaces.
  without_build_warnings <- function(expr) {
    withCallingHandlers(expr, warning = function(w) {
      if (grepl("built under R version", conditionMessage(w), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    })
  }
  set.seed(41)
  fx <- make_aggregate_fixture(n_per_group = 3, features = 3, sd = 0.3)

  serial <- dkge_aggregate_permute(fx$target, K = fx$K,
                                   statistic = "contrast_score",
                                   contrast = fx$row_contrast,
                                   B = 6, rank = 1, seed = 77)
  parallel <- without_build_warnings(
    dkge_aggregate_permute(fx$target, K = fx$K,
                           statistic = "contrast_score",
                           contrast = fx$row_contrast,
                           B = 6, rank = 1, seed = 77,
                           parallel = TRUE)
  )
  expect_equal(parallel$null, serial$null)
  expect_equal(parallel$p, serial$p)

  serial_boot <- suppressWarnings(
    dkge_aggregate_bootstrap(fx$target, K = fx$K, statistic = "contrast_score",
                             contrast = fx$row_contrast, return_features = TRUE,
                             B = 6, rank = 1, seed = 78)
  )
  parallel_boot <- without_build_warnings(suppressWarnings(
    dkge_aggregate_bootstrap(fx$target, K = fx$K, statistic = "contrast_score",
                             contrast = fx$row_contrast, return_features = TRUE,
                             B = 6, rank = 1, seed = 78, parallel = TRUE)
  ))
  expect_equal(parallel_boot$statistics, serial_boot$statistics)
  expect_equal(parallel_boot$feature_bootstrap$sd,
               serial_boot$feature_bootstrap$sd)
})

test_that("aggregate row IDs cannot be made ambiguous by a label separator", {
  set.seed(41)
  subject_data <- data.frame(
    subject_id = paste0("s", 1:4),
    # "a:b" + "c" and "a" + "b:c" both paste to "a:b:c".
    group = c("a:b", "a:b", "a", "a"),
    stringsAsFactors = FALSE
  )
  cell_data <- data.frame(task = c("c", "b:c"), stringsAsFactors = FALSE)
  values <- lapply(seq_len(4), function(i) matrix(stats::rnorm(2 * 3), 2, 3))
  names(values) <- subject_data$subject_id

  expect_error(
    dkge_aggregate_target(values, subject_data, cell_data,
                          group_vars = "group", cell_vars = "task"),
    "containing the row separator"
  )

  # A separator that cannot occur in the labels is accepted, and the four
  # distinct rows keep four distinct IDs.
  ok <- dkge_aggregate_target(values, subject_data, cell_data,
                              group_vars = "group", cell_vars = "task",
                              row_sep = "\r")
  expect_equal(nrow(ok$Y), 4L)
  expect_false(anyDuplicated(rownames(ok$Y)) > 0L)

  # Renaming the offending levels is the other way out.
  clean <- subject_data
  clean$group <- c("ab", "ab", "a", "a")
  clean_cells <- data.frame(task = c("c", "bc"), stringsAsFactors = FALSE)
  ok2 <- dkge_aggregate_target(values, clean, clean_cells,
                               group_vars = "group", cell_vars = "task")
  expect_equal(nrow(ok2$Y), 4L)
  expect_false(anyDuplicated(rownames(ok2$Y)) > 0L)

  # A cell label alone is enough to trigger the check.
  expect_error(
    dkge_aggregate_target(values, clean, cell_data,
                          group_vars = "group", cell_vars = "task"),
    "Cell variable 'task' has level"
  )

  # `row_sep` itself is validated.
  expect_error(
    dkge_aggregate_target(values, clean, clean_cells, group_vars = "group",
                          cell_vars = "task", row_sep = c(":", "-")),
    "single non-empty string"
  )
  expect_error(
    dkge_aggregate_target(values, clean, clean_cells, group_vars = "group",
                          cell_vars = "task", row_sep = ""),
    "single non-empty string"
  )

  # The row-ID builder refuses a collision even if one ever reached it.
  meta <- data.frame(group = c("a:b", "a"), task = c("c", "b:c"),
                     stringsAsFactors = FALSE)
  expect_error(
    dkge:::.dkge_aggregate_row_ids(meta, c("group", "task"), ":"),
    "must be unique and non-empty"
  )
})

test_that("bootstrap strata must nest within the target's group variables", {
  set.seed(43)
  fx <- make_aggregate_fixture(n_per_group = 3, features = 3, sd = 0.2)
  target <- fx$target
  subjects <- target$resample_spec$subject_ids

  # A single stratum spanning both groups can drop a whole group level, which
  # would change the aggregate row set mid-run and surface as a kernel error.
  expect_error(
    dkge_aggregate_bootstrap(target, K = fx$K, statistic = "singular_value",
                             strata = rep("all", length(subjects)),
                             B = 5, rank = 1, seed = 11),
    "must nest within"
  )

  # Cross-cutting strata are rejected for the same reason.
  expect_error(
    dkge_aggregate_bootstrap(target, K = fx$K, statistic = "singular_value",
                             strata = rep(c("x", "y"), length.out = length(subjects)),
                             B = 5, rank = 1, seed = 11),
    "must nest within"
  )

  # A stratification finer than `group_vars` is fine. (n_per_group = 3, so the
  # split below leaves strata of size 2 and 1; the singleton warning is a
  # different check, tested above.)
  finer <- paste(target$resample_spec$subject_data$group,
                 rep(c("odd", "odd", "even"), length.out = length(subjects)),
                 sep = "_")
  boot <- suppressWarnings(
    dkge_aggregate_bootstrap(target, K = fx$K, statistic = "singular_value",
                             strata = finer, B = 5, rank = 1, seed = 11)
  )
  expect_equal(length(boot$statistics), 5L)

  # The default (interaction of group_vars) obviously nests.
  boot_default <- dkge_aggregate_bootstrap(target, K = fx$K,
                                           statistic = "singular_value",
                                           B = 5, rank = 1, seed = 11)
  expect_equal(length(boot_default$statistics), 5L)

  # And the per-draw guard reports the row set, not a kernel dimname mismatch.
  dropped <- list(row_ids = setdiff(target$row_ids, target$row_ids[[1]]))
  expect_error(
    dkge:::.dkge_aggregate_check_row_ids(dropped, target$row_ids, "bootstrap"),
    "must nest within"
  )
  expect_error(
    dkge:::.dkge_aggregate_check_row_ids(dropped, target$row_ids, "bootstrap"),
    target$row_ids[[1]],
    fixed = TRUE
  )
  expect_error(
    dkge:::.dkge_aggregate_check_row_ids(dropped, target$row_ids, "permute"),
    "permute all row-defining subject variables jointly"
  )
})

test_that("misspelled resampling arguments are errors, not silent no-ops", {
  set.seed(45)
  fx <- make_aggregate_fixture(n_per_group = 3, features = 3, sd = 0.2)
  target <- fx$target

  # Every setting is a named formal, so a typo lands in `...` and would
  # previously have been swallowed by dkge_aggregate_stat()'s own `...`.
  expect_error(
    dkge_aggregate_permute(target, K = fx$K, statistic = "singular_value",
                           B = 3, rank = 1, sed = 1),
    "Unused argument"
  )
  expect_error(
    dkge_aggregate_bootstrap(target, K = fx$K, statistic = "singular_value",
                             B = 3, rank = 1, sed = 1),
    "Unused argument"
  )
  expect_error(
    dkge_aggregate_stat(dkge_aggregate_fit(target, K = fx$K, rank = 1),
                        statistic = "singular_value", nonsense = TRUE),
    "Unused argument"
  )

  # `contrast` and `component` are the two that legitimately apply.
  fit <- dkge_aggregate_fit(target, K = fx$K, rank = 2)
  expect_silent(dkge_aggregate_stat(fit, statistic = "salience_contrast",
                                    component = 2L,
                                    contrast = fx$row_contrast))

  # User-supplied statistic functions still receive `...`.
  scaled <- function(fit, scale = 1) scale * fit$singular_values[[1]]
  perm1 <- dkge_aggregate_permute(target, K = fx$K, statistic = scaled,
                                  B = 3, rank = 1, seed = 3)
  perm2 <- dkge_aggregate_permute(target, K = fx$K, statistic = scaled,
                                  B = 3, rank = 1, seed = 3, scale = 2)
  expect_equal(perm2$observed, 2 * perm1$observed)
  expect_equal(perm2$null, 2 * perm1$null)
})

test_that("resampling workers capture only what a draw needs", {
  set.seed(47)
  # Enough features that the stacked value matrix dominates the sizes below.
  fx <- make_aggregate_fixture(n_per_group = 3, features = 200, sd = 0.2)
  target <- fx$target
  spec <- target$resample_spec
  observed_fit <- dkge_aggregate_fit(target, K = fx$K, rank = 2)
  reference <- dkge:::.dkge_aggregate_alignment_reference(observed_fit)

  # The alignment reference keeps the three fields dkge_aggregate_align() reads
  # and drops the payloads that made the observed fit a second copy of the
  # stacked value matrix.
  expect_true(all(c("U", "K", "rank") %in% names(reference)))
  expect_null(reference$target)
  expect_null(reference$Y)
  expect_null(reference$Chat)
  expect_null(reference$scores_feature)
  expect_s3_class(reference, "dkge_aggregate_fit")

  # ... and it still aligns exactly like the full fit does.
  boot_idx <- c(1L, 1L, 2L, 4L, 5L, 6L)
  target_b <- dkge:::.dkge_aggregate_from_spec(spec, subject_indices = boot_idx)
  fit_b <- dkge_aggregate_fit(target_b, K = observed_fit$K,
                              rank = observed_fit$rank, center = "none")
  expect_equal(dkge_aggregate_align(reference, fit_b)$U,
               dkge_aggregate_align(observed_fit, fit_b)$U)

  worker <- dkge:::.dkge_aggregate_bootstrap_worker(
    spec = spec, reference = reference, row_ids = target$row_ids,
    center = "none", statistic = "singular_value",
    stat_args = list(component = 1L), component_contrasts = NULL,
    component_scale = "score", want_contrast = FALSE, want_features = FALSE
  )
  env <- environment(worker)
  expect_setequal(ls(env),
                  c("spec", "reference", "row_ids", "center", "statistic",
                    "stat_args", "component_contrasts", "component_scale",
                    "want_contrast", "want_features"))
  # The stacked values are unavoidable; nothing else of consequence rides along.
  exported <- length(serialize(as.list(env), NULL))
  spec_size <- length(serialize(spec, NULL))
  expect_lt(exported, 1.15 * spec_size)

  # What an inline closure over the whole calling frame would have exported to
  # every parallel chunk: the target and the observed fit each carry their own
  # reference to `spec`, plus the pre-drawn index list and the accumulators.
  naive <- length(serialize(list(spec = spec, target = target,
                                 observed_fit = observed_fit,
                                 boot_idx = replicate(64, boot_idx,
                                                      simplify = FALSE),
                                 stats = numeric(64),
                                 alignments = vector("list", 64)), NULL))
  expect_lt(exported, 0.5 * naive)

  perm_worker <- dkge:::.dkge_aggregate_permute_worker(
    spec = spec, reference = reference, row_ids = target$row_ids,
    group_vars = "group", center = "none", statistic = "singular_value",
    stat_args = list(component = 1L)
  )
  expect_setequal(ls(environment(perm_worker)),
                  c("spec", "reference", "row_ids", "group_vars", "center",
                    "statistic", "stat_args"))
  expect_lt(length(serialize(as.list(environment(perm_worker)), NULL)),
            1.15 * spec_size)

  # The worker still produces exactly what the entry point reports.
  boot <- dkge_aggregate_bootstrap(target, K = fx$K,
                                   statistic = "singular_value",
                                   B = 4, rank = 2, seed = 909)
  set.seed(909)
  strata_index <- split(seq_along(spec$subject_ids),
                        droplevels(as.factor(spec$subject_data$group)))
  idx <- lapply(seq_len(4), function(b) {
    dkge:::.dkge_aggregate_bootstrap_indices(strata_index)
  })
  expect_equal(vapply(idx, function(i) worker(i)$stat, numeric(1)),
               boot$statistics)
})

test_that("shared seed and B helpers behave identically in both layers", {
  set.seed(51)
  fx <- make_aggregate_fixture(n_per_group = 3, features = 3, sd = 0.2)
  target <- fx$target

  for (bad in list(0L, -1L, NA_integer_, c(2L, 3L))) {
    expect_error(dkge_aggregate_permute(target, K = fx$K,
                                        statistic = "singular_value",
                                        B = bad, rank = 1),
                 "`B` must be a positive integer")
    expect_error(dkge_aggregate_bootstrap(target, K = fx$K,
                                          statistic = "singular_value",
                                          B = bad, rank = 1),
                 "`B` must be a positive integer")
  }

  # Both layers restore the caller's RNG stream after a seeded run.
  set.seed(1234)
  before <- .Random.seed
  dkge_aggregate_permute(target, K = fx$K, statistic = "singular_value",
                         B = 3, rank = 1, seed = 99)
  expect_identical(.Random.seed, before)
  dkge_aggregate_bootstrap(target, K = fx$K, statistic = "singular_value",
                           B = 3, rank = 1, seed = 99)
  expect_identical(.Random.seed, before)

  # ... including the case where the caller had no RNG stream at all.
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    rm(".Random.seed", envir = globalenv())
  }
  dkge_aggregate_permute(target, K = fx$K, statistic = "singular_value",
                         B = 3, rank = 1, seed = 99)
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
})
