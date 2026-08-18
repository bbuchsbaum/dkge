# Design-Kernel Group Embedding on the SDAM first-level MVPA AUC maps.
#
# Mirrors plsrri/inst/examples/sdam_firstlevel_task_pls.R so the same dataset
# (32 subjects, 4 MVPA AUC maps each: task in {nback, recog} crossed with
# measure in {low_mid, high_sem}) can be reanalyzed with DKGE. The maps are
# AUC - 0.5, i.e. activity-like classification scores already centred at zero,
# not GLM betas. DKGE accepts them as the rows of B_s without modification.
#
# What this script demonstrates that classical task PLS does not:
#   * Explicit factorial design kernel `K` over the 2x2 cells (task, measure,
#     interaction) so components rotate to respect the design, not just raw
#     covariance.
#   * Leave-one-subject-out debiased contrasts via `dkge_contrast(method =
#     "loso")` for the task and measure main effects and the interaction.
#   * A subject-level GLM (`~ group`) on the latent target via
#     `dkge_between_rrr()` + `dkge_between_permute()` instead of mean-centering
#     groups inside SVD.
#   * Equal subject block weighting (`w_method = "none"`) so AUC-minus-chance
#     maps keep their absolute scale.
#   * BSR-equivalent z-maps from a small q-space bootstrap.
#
# Usage (from the dkge repo root):
#   Rscript inst/examples/sdam_firstlevel_task_dkge.R
# Override defaults with environment variables:
#   DKGE_SDAM_TESTDATA   path to the testdata folder (default: testdata/)
#   DKGE_SDAM_NPERM      between-subject permutation reps (default: 999)
#   DKGE_SDAM_NBOOT      bootstrap reps for BSR-like maps (default: 200)
#   DKGE_SDAM_MAX_VOXELS subset of voxels for fast development (default: Inf)

SDAM_TASKS    <- c("recog", "nback")
SDAM_MEASURES <- c("low_mid", "high_sem")
SDAM_GROUPS   <- c("control", "sdam")
SDAM_SEED     <- 20260509L

sdam_time <- function(label, expr) {
  message(sprintf("[%s] %s...", format(Sys.time(), "%H:%M:%S"), label))
  elapsed <- system.time(value <- force(expr))
  message(sprintf("[%s] %s done (%.1fs)",
                  format(Sys.time(), "%H:%M:%S"), label, elapsed[["elapsed"]]))
  value
}

# -----------------------------------------------------------------------------
# Manifest construction (parallel to the PLS example)
# -----------------------------------------------------------------------------

sdam_testdata_root <- function(root = Sys.getenv("DKGE_SDAM_TESTDATA", "testdata")) {
  root <- normalizePath(root, winslash = "/", mustWork = FALSE)
  if (!dir.exists(root)) {
    stop(
      "Could not find the SDAM testdata directory. ",
      "Set DKGE_SDAM_TESTDATA or run from the dkge repository root.",
      call. = FALSE
    )
  }
  root
}

sdam_cell_table <- function() {
  data.frame(
    task    = rep(SDAM_TASKS, each = length(SDAM_MEASURES)),
    measure = rep(SDAM_MEASURES, times = length(SDAM_TASKS)),
    stringsAsFactors = FALSE
  ) |>
    transform(cell = paste(task, measure, sep = "_"))
}

read_sdam_participants <- function(root = sdam_testdata_root()) {
  participants <- utils::read.table(
    file.path(root, "participants.tsv"),
    header = TRUE, sep = "", quote = "",
    comment.char = "", stringsAsFactors = FALSE
  )
  participants <- participants[, c("participant_id", "group")]
  names(participants)[names(participants) == "participant_id"] <- "subject"
  participants$subject <- as.character(participants$subject)
  participants$group   <- factor(participants$group, levels = SDAM_GROUPS)
  participants <- participants[order(participants$group, participants$subject), ]
  rownames(participants) <- NULL
  participants
}

build_sdam_manifest <- function(root = sdam_testdata_root()) {
  root <- sdam_testdata_root(root)
  participants <- read_sdam_participants(root)
  cells <- sdam_cell_table()

  rows <- vector("list", nrow(participants) * nrow(cells))
  k <- 1L
  for (subject_id in participants$subject) {
    group <- as.character(participants$group[participants$subject == subject_id])
    for (i in seq_len(nrow(cells))) {
      rows[[k]] <- data.frame(
        subject = subject_id,
        group   = group,
        task    = cells$task[[i]],
        measure = cells$measure[[i]],
        cell    = cells$cell[[i]],
        file    = file.path(
          root,
          sprintf("%s_vgg_model_rsa_%s", subject_id, cells$task[[i]]),
          "maps",
          sprintf("%s.nii.gz", cells$measure[[i]])
        ),
        stringsAsFactors = FALSE
      )
      k <- k + 1L
    }
  }

  manifest <- do.call(rbind, rows)
  manifest$group <- factor(manifest$group, levels = SDAM_GROUPS)
  manifest$cell  <- factor(manifest$cell,  levels = cells$cell)

  missing_files <- manifest$file[!file.exists(manifest$file)]
  if (length(missing_files) > 0L) {
    stop(
      "Manifest references missing maps: ",
      paste(utils::head(missing_files, 5), collapse = ", "),
      call. = FALSE
    )
  }
  manifest
}

# -----------------------------------------------------------------------------
# Mask + B_s assembly. Each subject contributes a 4 x P matrix whose rows are
# the four cells in canonical order. This is the q x P_s shape that DKGE wants;
# the per-subject "design" is the q x q identity (one observation per cell).
# -----------------------------------------------------------------------------

derive_common_nonzero_mask <- function(manifest, max_voxels = Inf) {
  if (!requireNamespace("neuroim2", quietly = TRUE)) {
    stop("Package 'neuroim2' is required to derive the SDAM mask.", call. = FALSE)
  }
  files <- unique(as.character(manifest$file))
  first_map <- neuroim2::read_vol(files[[1L]])
  first_arr <- as.array(first_map)
  common <- is.finite(first_arr) & first_arr != 0

  for (file in files[-1L]) {
    arr <- as.array(neuroim2::read_vol(file))
    if (!identical(dim(arr), dim(common))) {
      stop("All SDAM maps must share dimensions.", call. = FALSE)
    }
    common <- common & is.finite(arr) & arr != 0
  }

  if (is.finite(max_voxels)) {
    keep <- utils::head(which(common), as.integer(max_voxels))
    common[] <- FALSE
    common[keep] <- TRUE
  }
  if (!sum(common)) stop("Empty SDAM mask.", call. = FALSE)

  mask_arr <- array(as.integer(common), dim = dim(common))
  neuroim2::NeuroVol(mask_arr, neuroim2::space(first_map))
}

load_sdam_betas <- function(manifest, mask) {
  mask_idx <- which(as.array(mask) > 0)
  cells    <- levels(manifest$cell)
  q        <- length(cells)
  subjects <- unique(as.character(manifest$subject))

  B_list <- vector("list", length(subjects))
  names(B_list) <- subjects
  for (sid in subjects) {
    sub_rows <- manifest[manifest$subject == sid, , drop = FALSE]
    B <- matrix(NA_real_, nrow = q, ncol = length(mask_idx),
                dimnames = list(cells, NULL))
    for (i in seq_len(nrow(sub_rows))) {
      arr <- as.array(neuroim2::read_vol(sub_rows$file[[i]]))
      B[as.character(sub_rows$cell[[i]]), ] <- arr[mask_idx]
    }
    B[!is.finite(B)] <- 0
    B_list[[sid]] <- B
  }
  list(B_list = B_list, cells = cells)
}

# Pseudo-design: one observation per cell. The pooled Cholesky reduces to
# sqrt(S) * I, so row standardisation is a uniform scale -- the design kernel,
# not the design matrix, carries the structural information here.
make_sdam_design_list <- function(B_list) {
  q <- nrow(B_list[[1]])
  design_template <- diag(q)
  colnames(design_template) <- rownames(B_list[[1]])
  rep(list(design_template), length(B_list))
}

# -----------------------------------------------------------------------------
# Design kernel: 2x2 factorial (task x measure). All three terms get an
# explicit weight; the unit-trace normalisation keeps the kernel scale-stable
# when rho changes.
# -----------------------------------------------------------------------------

make_sdam_design_kernel <- function(rho_task = 1, rho_measure = 1, rho_int = 1) {
  kern <- dkge::design_kernel(
    factors = list(
      task    = list(L = length(SDAM_TASKS),    type = "nominal"),
      measure = list(L = length(SDAM_MEASURES), type = "nominal")
    ),
    terms = list("task", "measure", c("task", "measure")),
    rho   = c(task = rho_task, measure = rho_measure, "task:measure" = rho_int),
    basis = "cell",
    normalize = "unit_trace"
  )
  cells <- sdam_cell_table()$cell
  rownames(kern$K) <- cells
  colnames(kern$K) <- cells
  kern
}

# -----------------------------------------------------------------------------
# Fit
# -----------------------------------------------------------------------------

fit_sdam_dkge <- function(B_list, kernel, rank = 3L,
                          w_method = c("none", "mfa_sigma1", "energy"),
                          effect_scaling = c("none", "pooled_design")) {
  w_method <- match.arg(w_method)
  effect_scaling <- match.arg(effect_scaling)
  designs  <- make_sdam_design_list(B_list)
  data     <- dkge::dkge_data(B_list, designs = designs,
                              subject_ids = names(B_list))
  fit <- dkge::dkge(
    data,
    K = kernel,
    rank = rank,
    w_method = w_method,
    effect_scaling = effect_scaling,
    keep_inputs = TRUE
  )
  fit
}

# -----------------------------------------------------------------------------
# Cell-space contrasts (task / measure / interaction). The cell ordering matches
# `sdam_cell_table()`: recog_low_mid, recog_high_sem, nback_low_mid, nback_high_sem.
# -----------------------------------------------------------------------------

sdam_contrast_matrix <- function() {
  # Rows: cells in canonical order. Columns: contrasts of interest.
  C <- cbind(
    task            = c(-0.5, -0.5,  0.5,  0.5),  # nback - recog
    measure         = c(-0.5,  0.5, -0.5,  0.5),  # high_sem - low_mid
    `task:measure`  = c( 0.5, -0.5, -0.5,  0.5),  # interaction
    grand_mean      = c( 0.25, 0.25, 0.25, 0.25)
  )
  rownames(C) <- sdam_cell_table()$cell
  C
}

make_sdam_group_effect_grid <- function() {
  dkge::dkge_effect_grid(
    factors = list(
      group   = SDAM_GROUPS,
      task    = SDAM_TASKS,
      measure = SDAM_MEASURES
    ),
    scope = c(group = "between", task = "within", measure = "within"),
    block_factors = "group"
  )
}

make_sdam_group_design_kernel <- function(rho_group = 1,
                                          rho_task = 1,
                                          rho_measure = 1,
                                          rho_group_task = 1,
                                          rho_group_measure = 1,
                                          rho_task_measure = 1,
                                          rho_group_task_measure = 1) {
  grid <- make_sdam_group_effect_grid()
  dkge::design_kernel(
    grid,
    terms = list(
      "group", "task", "measure",
      c("group", "task"),
      c("group", "measure"),
      c("task", "measure"),
      c("group", "task", "measure")
    ),
    rho = c(
      group = rho_group,
      task = rho_task,
      measure = rho_measure,
      "group:task" = rho_group_task,
      "group:measure" = rho_group_measure,
      "task:measure" = rho_task_measure,
      "group:task:measure" = rho_group_task_measure
    ),
    basis = "cell",
    normalize = "unit_trace"
  )
}

make_sdam_cell_mean_bridge_kernel <- function(rho_group = 1,
                                              rho_task = 1,
                                              rho_measure = 1,
                                              rho_group_task = 1,
                                              rho_group_measure = 1,
                                              rho_task_measure = 1,
                                              rho_group_task_measure = 1) {
  grid <- dkge::dkge_effect_grid(
    factors = list(
      group   = SDAM_GROUPS,
      task    = SDAM_TASKS,
      measure = SDAM_MEASURES
    ),
    scope = c(group = "between", task = "within", measure = "within")
  )
  dkge::design_kernel(
    grid,
    terms = list(
      "group", "task", "measure",
      c("group", "task"),
      c("group", "measure"),
      c("task", "measure"),
      c("group", "task", "measure")
    ),
    rho = c(
      group = rho_group,
      task = rho_task,
      measure = rho_measure,
      "group:task" = rho_group_task,
      "group:measure" = rho_group_measure,
      "task:measure" = rho_task_measure,
      "group:task:measure" = rho_group_task_measure
    ),
    basis = "cell",
    normalize = "unit_trace"
  )
}

make_sdam_group_explicit_data <- function(B_list, manifest) {
  grid <- make_sdam_group_effect_grid()
  participants <- unique(manifest[, c("subject", "group")])
  participants$group <- factor(participants$group, levels = SDAM_GROUPS)
  participants <- participants[order(participants$group, participants$subject), ]
  subjects <- intersect(as.character(participants$subject), names(B_list))
  cell_meta <- sdam_cell_table()

  betas <- lapply(subjects, function(subject_id) {
    group <- as.character(participants$group[participants$subject == subject_id][[1]])
    rows <- paste(group, cell_meta$task, cell_meta$measure, sep = ":")
    B <- B_list[[subject_id]][cell_meta$cell, , drop = FALSE]
    rownames(B) <- rows
    B
  })
  names(betas) <- subjects
  designs <- lapply(betas, function(B) {
    X <- diag(nrow(B))
    colnames(X) <- rownames(B)
    X
  })

  list(
    data = dkge::dkge_data(betas, designs = designs, subject_ids = subjects),
    grid = grid,
    participants = participants[match(subjects, participants$subject), ]
  )
}

sdam_group_contrast_matrix <- function() {
  grid <- make_sdam_group_effect_grid()
  C4 <- sdam_contrast_matrix()
  C <- cbind(
    group = c(rep(-0.25, 4), rep(0.25, 4)),
    `group:task` = c(-C4[, "task"], C4[, "task"]),
    `group:measure` = c(-C4[, "measure"], C4[, "measure"]),
    `group:task:measure` = c(-C4[, "task:measure"], C4[, "task:measure"])
  )
  rownames(C) <- grid$cell_labels
  C
}

sdam_cell_mean_bridge_contrast_matrix <- function() {
  grid <- make_sdam_group_effect_grid()
  C4 <- sdam_contrast_matrix()
  C <- cbind(
    task = c(C4[, "task"], C4[, "task"]) / 2,
    measure = c(C4[, "measure"], C4[, "measure"]) / 2,
    `task:measure` = c(C4[, "task:measure"], C4[, "task:measure"]) / 2,
    group = c(rep(-0.25, 4), rep(0.25, 4)),
    `group:task` = c(-C4[, "task"], C4[, "task"]) / 2,
    `group:measure` = c(-C4[, "measure"], C4[, "measure"]) / 2,
    `group:task:measure` = c(-C4[, "task:measure"], C4[, "task:measure"]) / 2
  )
  rownames(C) <- grid$cell_labels
  C
}

run_sdam_group_explicit_validation <- function(B_list, manifest,
                                               rank = 3L,
                                               folds = 4L,
                                               seed = SDAM_SEED,
                                               w_method = c("none", "mfa_sigma1", "energy"),
                                               effect_scaling = c("none", "pooled_design")) {
  w_method <- match.arg(w_method)
  effect_scaling <- match.arg(effect_scaling)
  payload <- make_sdam_group_explicit_data(B_list, manifest)
  kernel <- make_sdam_group_design_kernel()
  fit <- dkge::dkge(
    payload$data,
    K = kernel,
    rank = rank,
    w_method = w_method,
    effect_scaling = effect_scaling,
    keep_inputs = TRUE,
    missingness = "mask",
    miss_args = list(min_pairs = 1L)
  )
  contrasts <- sdam_group_contrast_matrix()
  set.seed(seed)
  validated <- dkge::dkge_contrast_validated(
    fit,
    contrasts = contrasts,
    folds = folds,
    align = FALSE,
    observed_missingness = "mask",
    observed_args = list(min_pairs = 1L),
    completed_missingness = "shrink",
    completed_args = list(gamma = 1)
  )

  list(
    fit = fit,
    kernel = kernel,
    data = payload$data,
    participants = payload$participants,
    w_method = w_method,
    contrasts = contrasts,
    validation = validated,
    component_summary = summarise_sdam_dkge(fit),
    coverage = payload$data$provenance$coverage,
    pair_counts = payload$data$provenance$pair_counts,
    cell_labels = payload$grid$cell_labels,
    term_scope = kernel$info$term_scope,
    block_factors = kernel$info$block_factors
  )
}

sdam_q8_cell_meta <- function(labels) {
  parts <- strsplit(as.character(labels), ":", fixed = TRUE)
  lens <- lengths(parts)
  if (any(lens != 3L)) {
    stop("q=8 cell labels must have form group:task:measure.", call. = FALSE)
  }
  mat <- do.call(rbind, parts)
  cell_lookup <- sdam_cell_table()
  cell <- paste(mat[, 2], mat[, 3], sep = "_")
  cell <- cell_lookup$cell[match(cell, paste(cell_lookup$task, cell_lookup$measure, sep = "_"))]
  data.frame(
    group = mat[, 1],
    task = mat[, 2],
    measure = mat[, 3],
    cell = cell,
    stringsAsFactors = FALSE
  )
}

sdam_q8_design_loading_table <- function(fit) {
  load_mat <- fit$K %*% fit$U
  labels <- rownames(load_mat) %||% rownames(fit$K) %||% paste0("cell", seq_len(nrow(load_mat)))
  meta <- sdam_q8_cell_meta(labels)
  rows <- lapply(seq_len(ncol(load_mat)), function(j) {
    data.frame(
      component = paste0("LV", j),
      group = meta$group,
      task = meta$task,
      measure = meta$measure,
      cell = labels,
      loading = as.numeric(load_mat[, j]),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$component <- factor(out$component, levels = paste0("LV", seq_len(ncol(load_mat))))
  out$group <- factor(out$group, levels = SDAM_GROUPS)
  out$task <- factor(out$task, levels = SDAM_TASKS)
  out$measure <- factor(out$measure, levels = SDAM_MEASURES)
  out
}

sdam_subject_weight_table <- function(fit, participants) {
  subject_ids <- fit$subject_ids %||% names(fit$weights)
  if (is.null(subject_ids)) subject_ids <- paste0("subject", seq_along(fit$weights))
  weights <- fit$weights %||% rep(1, length(subject_ids))
  out <- data.frame(
    subject = as.character(subject_ids),
    weight = as.numeric(weights),
    stringsAsFactors = FALSE
  )
  if (!is.null(participants)) {
    out$group <- participants$group[match(out$subject, participants$subject)]
  }
  out
}

sdam_subject_weight_summary <- function(weight_table) {
  rows <- lapply(split(weight_table, weight_table$group), function(d) {
    data.frame(
      group = as.character(d$group[[1]]),
      n = nrow(d),
      mean = mean(d$weight),
      sd = stats::sd(d$weight),
      min = min(d$weight),
      median = stats::median(d$weight),
      max = max(d$weight),
      total = sum(d$weight),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

sdam_component_signs_to_reference <- function(fit, reference_fit) {
  if (is.null(reference_fit)) return(rep(1, ncol(fit$U)))
  K <- fit$K
  cross <- t(reference_fit$U) %*% K %*% fit$U
  signs <- sign(diag(cross))
  signs[!is.finite(signs) | signs == 0] <- 1
  signs
}

sdam_q8_component_cell_scores <- function(fit, participants, signs = NULL) {
  common_loadings <- sdam_common_component_loadings(fit)
  if (is.null(signs)) signs <- rep(1, ncol(common_loadings))
  rows <- vector("list", length(fit$Btil) * ncol(common_loadings))
  obs_mask <- fit$provenance$obs_mask %||% vector("list", length(fit$Btil))
  k <- 1L
  for (s in seq_along(fit$Btil)) {
    Bts <- fit$Btil[[s]]
    keep <- obs_mask[[s]]
    if (is.null(keep)) keep <- rowSums(abs(Bts)) > 0
    keep <- as.logical(keep)
    scores <- Bts[keep, , drop = FALSE] %*% common_loadings
    scores <- sweep(scores, 2L, signs, "*")
    labels <- rownames(Bts)[keep] %||% rownames(fit$K)[keep]
    meta <- sdam_q8_cell_meta(labels)
    sid <- fit$subject_ids[[s]] %||% paste0("sub", s)
    subject_group <- participants$group[match(sid, participants$subject)]
    for (j in seq_len(ncol(scores))) {
      rows[[k]] <- data.frame(
        subject = sid,
        group = subject_group,
        cell_group = meta$group,
        cell = meta$cell,
        task = meta$task,
        measure = meta$measure,
        component = paste0("LV", j),
        score = scores[, j],
        stringsAsFactors = FALSE
      )
      k <- k + 1L
    }
  }
  df <- do.call(rbind, rows)
  df$group <- factor(df$group, levels = SDAM_GROUPS)
  df$cell_group <- factor(df$cell_group, levels = SDAM_GROUPS)
  df$cell <- factor(df$cell, levels = sdam_cell_table()$cell)
  df$task <- factor(df$task, levels = SDAM_TASKS)
  df$measure <- factor(df$measure, levels = SDAM_MEASURES)
  df$component <- factor(df$component, levels = paste0("LV", seq_len(ncol(common_loadings))))
  df
}

sdam_q8_fit_summary <- function(fit, participants, w_method, signs = NULL) {
  scores <- sdam_q8_component_cell_scores(fit, participants, signs = signs)
  weights <- sdam_subject_weight_table(fit, participants)
  list(
    w_method = w_method,
    component_signs = signs %||% rep(1, ncol(fit$U)),
    component_summary = summarise_sdam_dkge(fit),
    score_values = scores,
    score_summary = sdam_component_score_summary(scores),
    score_tests = sdam_component_score_tests(scores),
    weight_table = weights,
    weight_summary = sdam_subject_weight_summary(weights)
  )
}

run_sdam_group_explicit_sensitivity <- function(B_list, manifest,
                                                rank = 3L,
                                                methods = c("none", "mfa_sigma1", "energy"),
                                                reference = NULL) {
  out <- vector("list", length(methods))
  names(out) <- methods
  reference_fit <- if (is.null(reference)) NULL else reference$fit
  for (method in methods) {
    if (!is.null(reference) && identical(method, reference$w_method)) {
      fit <- reference$fit
      participants <- reference$participants
    } else {
      payload <- make_sdam_group_explicit_data(B_list, manifest)
      kernel <- make_sdam_group_design_kernel()
      fit <- dkge::dkge(
        payload$data,
        K = kernel,
        rank = rank,
        w_method = method,
        effect_scaling = "none",
        keep_inputs = TRUE,
        missingness = "mask",
        miss_args = list(min_pairs = 1L)
      )
      participants <- payload$participants
    }
    signs <- sdam_component_signs_to_reference(fit, reference_fit)
    out[[method]] <- sdam_q8_fit_summary(fit, participants, method, signs = signs)
  }
  out
}

compute_sdam_loso_contrasts <- function(fit, contrasts = NULL,
                                        method = c("loso", "analytic")) {
  method <- match.arg(method)
  if (is.null(contrasts)) contrasts <- sdam_contrast_matrix()
  dkge::dkge_contrast(fit, contrasts, method = method, align = TRUE)
}

# Group-level voxel summaries from a `dkge_contrasts` object: mean and a
# per-voxel z (mean / SE across subjects). These are the DKGE analogues of
# PLS bootstrap-ratio maps.
sdam_contrast_group_maps <- function(contrasts_obj, manifest) {
  out <- list()
  for (nm in names(contrasts_obj$values)) {
    M <- do.call(rbind, contrasts_obj$values[[nm]])  # subjects x voxels
    rownames(M) <- names(contrasts_obj$values[[nm]])
    n  <- nrow(M)
    mu <- colMeans(M)
    sd <- apply(M, 2L, stats::sd)
    se <- sd / sqrt(n)
    z  <- ifelse(se > 0, mu / se, 0)
    out[[nm]] <- list(mean = mu, sd = sd, se = se, z = z, n = n,
                      subject_ids = rownames(M))
  }
  out
}

# Group difference (sdam - control) on a per-voxel basis from LOSO contrasts.
sdam_contrast_group_diff <- function(contrasts_obj, manifest) {
  meta <- unique(manifest[, c("subject", "group")])
  out <- list()
  for (nm in names(contrasts_obj$values)) {
    M <- do.call(rbind, contrasts_obj$values[[nm]])
    rownames(M) <- names(contrasts_obj$values[[nm]])
    grp <- meta$group[match(rownames(M), meta$subject)]
    keep <- !is.na(grp)
    M <- M[keep, , drop = FALSE]; grp <- grp[keep]
    A <- M[grp == "sdam", , drop = FALSE]
    B <- M[grp == "control", , drop = FALSE]
    if (!nrow(A) || !nrow(B)) next
    mu_diff <- colMeans(A) - colMeans(B)
    pooled_se <- sqrt(apply(A, 2L, stats::var) / nrow(A) +
                      apply(B, 2L, stats::var) / nrow(B))
    z <- ifelse(pooled_se > 0, mu_diff / pooled_se, 0)
    out[[nm]] <- list(mean_diff = mu_diff, z = z,
                      n_sdam = nrow(A), n_control = nrow(B))
  }
  out
}

# -----------------------------------------------------------------------------
# Between-subjects DKGE: cross-fitted voxel contrast maps feed a reduced-rank
# regression on `~ group`. The target is subjects x (contrast blocks x voxels),
# so group is tested on spatial expressions of the repeated-measures effects
# rather than on component means.
# -----------------------------------------------------------------------------

run_sdam_between_subjects <- function(fit, manifest,
                                       contrast_obj = NULL,
                                       contrasts = NULL,
                                       formula = ~ group,
                                       nperm = 999L,
                                       seed  = SDAM_SEED) {
  participants <- unique(manifest[, c("subject", "group")])
  rownames(participants) <- NULL

  if (is.null(contrast_obj)) {
    if (is.null(contrasts)) contrasts <- sdam_contrast_matrix()
    contrast_obj <- dkge::dkge_contrast(fit, contrasts, method = "loso", align = TRUE)
  }

  target <- dkge::dkge_make_target(
    fit,
    type = "transported_maps",
    contrast_obj = contrast_obj,
    crossfit = "loso"
  )
  design <- dkge::dkge_subject_model(formula, data = participants,
                                     subject_id_col = "subject")
  fit_b <- dkge::dkge_between_rrr(target, design, rank = 1L)

  perm <- dkge::dkge_between_permute(
    fit_b,
    terms = setdiff(design$term_names, "(Intercept)"),
    B = nperm, seed = seed,
    scope = "both", feature_adjust = "maxT"
  )

  feature_blocks <- sub(":[^:]+$", "", target$feature_ids)
  int_mask <- feature_blocks == "task:measure"
  focused <- NULL
  if (any(int_mask)) {
    fit_int <- dkge::dkge_between_rrr(target, design, rank = 1L,
                                      feature_mask = int_mask)
    perm_int <- dkge::dkge_between_permute(
      fit_int,
      terms = setdiff(design$term_names, "(Intercept)"),
      B = nperm, seed = seed,
      scope = "both", feature_adjust = "maxT"
    )
    focused <- list(
      contrast_block = "task:measure",
      feature_mask = int_mask,
      fit = fit_int,
      perm = perm_int
    )
  }

  list(target = target, design = design, fit = fit_b, perm = perm,
       focused = focused,
       participants = participants)
}

run_sdam_cell_mean_bridge <- function(B_list, manifest,
                                      rank = 3L,
                                      component = 2L,
                                      nperm = 999L,
                                      nboot = 200L,
                                      seed = SDAM_SEED) {
  participants <- unique(manifest[, c("subject", "group")])
  participants <- participants[match(names(B_list), participants$subject), , drop = FALSE]
  participants$group <- factor(participants$group, levels = SDAM_GROUPS)
  cell_meta <- sdam_cell_table()

  target <- dkge::dkge_aggregate_target(
    values = B_list,
    subject_data = participants,
    cell_data = cell_meta,
    group_vars = "group",
    cell_vars = c("task", "measure"),
    subject_id_col = "subject"
  )
  kernel <- make_sdam_cell_mean_bridge_kernel()
  bridge_contrasts <- sdam_cell_mean_bridge_contrast_matrix()
  bridge_contrasts <- bridge_contrasts[rownames(target$Y), , drop = FALSE]
  contrast <- bridge_contrasts[, "group:task:measure"]

  fit <- dkge::dkge_aggregate_fit(target, K = kernel, rank = rank,
                                  center = "column")
  stat <- dkge::dkge_aggregate_stat(
    fit,
    statistic = "between_group_contrast",
    component = component,
    contrast = contrast
  )
  perm <- dkge::dkge_aggregate_permute(
    target,
    K = kernel,
    statistic = "between_group_contrast",
    component = component,
    contrast = contrast,
    B = nperm,
    rank = rank,
    center = "column",
    seed = seed
  )
  boot <- dkge::dkge_aggregate_bootstrap(
    target,
    K = kernel,
    statistic = "between_group_contrast",
    component = component,
    contrast = contrast,
    component_contrasts = bridge_contrasts,
    component_scale = "score",
    return_features = TRUE,
    B = nboot,
    rank = rank,
    center = "column",
    seed = seed + 1L
  )

  list(
    target = target,
    kernel = kernel,
    fit = fit,
    contrast = contrast,
    component_contrasts = bridge_contrasts,
    statistic = stat,
    permutation = perm,
    bootstrap = boot,
    rank = rank,
    component = component,
    center = "column",
    estimand = "cell_mean_pls_bridge"
  )
}

sdam_aggregate_salience_table <- function(fit) {
  sal <- fit$saliences
  row_meta <- fit$target$row_metadata
  rows <- lapply(seq_len(ncol(sal)), function(j) {
    data.frame(
      component = paste0("LV", j),
      row_meta,
      loading = as.numeric(sal[, j]),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$component <- factor(out$component, levels = paste0("LV", seq_len(ncol(sal))))
  out$group <- factor(out$group, levels = SDAM_GROUPS)
  out$task <- factor(out$task, levels = SDAM_TASKS)
  out$measure <- factor(out$measure, levels = SDAM_MEASURES)
  out
}

run_sdam_stratified_basis_check <- function(B_list, kernel, manifest,
                                            rank = 3L) {
  participants <- unique(manifest[, c("subject", "group")])
  ids_by_group <- split(as.character(participants$subject), participants$group)

  # Two rank-r subspaces of a q-dimensional ambient space share at least
  # max(0, 2r - q) dimensions by the dimension formula, so the first
  # `forced_overlap` principal angles are geometrically zero independent
  # of the data. Cap the stratified rank at floor(q/2) so every reported
  # angle carries information.
  K <- if (is.list(kernel) && !is.null(kernel$K)) kernel$K else kernel
  q <- nrow(K)
  rank_strat <- max(1L, min(as.integer(rank), as.integer(q %/% 2L)))

  fits <- lapply(ids_by_group, function(ids) {
    fit_sdam_dkge(B_list[ids], kernel, rank = rank_strat,
                  w_method = "none", effect_scaling = "none")
  })
  if (length(fits) < 2L) {
    stop("Stratified basis check requires at least two groups.", call. = FALSE)
  }
  group_names <- names(fits)
  fit_a <- fits[[group_names[[1L]]]]
  fit_b <- fits[[group_names[[2L]]]]
  common_rank <- min(ncol(fit_a$U), ncol(fit_b$U))
  angles <- dkge::dkge_principal_angles_K(
    fit_a$U[, seq_len(common_rank), drop = FALSE],
    fit_b$U[, seq_len(common_rank), drop = FALSE],
    fit_a$K
  )
  forced <- max(0L, 2L * common_rank - q)
  angle_table <- data.frame(
    component = paste0("LV", seq_len(common_rank)),
    angle_rad = angles,
    angle_deg = angles * 180 / pi,
    group_a = group_names[[1L]],
    group_b = group_names[[2L]],
    stringsAsFactors = FALSE
  )
  list(
    fits = fits,
    angles = angle_table,
    rank_strat = rank_strat,
    rank_requested = as.integer(rank),
    q = q,
    forced_overlap = forced
  )
}

# -----------------------------------------------------------------------------
# BSR-equivalent: a small q-space bootstrap of a contrast. Returns z = mean/sd
# at every cluster/voxel, ready for an overlay map. Uses the projected
# bootstrap (no transport cache needed because we work in a single common mask).
# -----------------------------------------------------------------------------

sdam_bsr_voxels <- function(fit, contrast, B = 200L, seed = SDAM_SEED) {
  vals <- vector("list", length(fit$Btil))
  for (s in seq_along(fit$Btil)) {
    res <- dkge::dkge_loso_contrast(fit, s = s, contrasts = contrast)
    vals[[s]] <- res$v
  }
  names(vals) <- fit$subject_ids %||% paste0("sub", seq_along(vals))
  boot <- dkge::dkge_bootstrap_projected(vals, B = B, seed = seed,
                                          return_samples = FALSE)
  mu <- boot$medoid$mean
  sdv <- boot$medoid$sd
  list(z = mu / pmax(sdv, 1e-12), mean = mu, sd = sdv, B = B)
}

# -----------------------------------------------------------------------------
# Diagnostics + plotting (mirrors the PLS panel set: scree, design saliences,
# subject brain scores, BSR-like brain map). Surface plots are wired through
# neuroatlas::plot_brain when available.
# -----------------------------------------------------------------------------

sdam_values_to_neurovol <- function(values, mask) {
  values <- as.numeric(values)
  mask_arr <- as.array(mask) > 0
  if (length(values) != sum(mask_arr)) {
    stop("Value vector length does not match mask size.", call. = FALSE)
  }
  vol <- array(NA_real_, dim = dim(mask))
  vol[mask_arr] <- values
  neuroim2::NeuroVol(vol, neuroim2::space(mask))
}

plot_sdam_design_heatmap <- function(fit) {
  # Mirrors PLS's design heatmap. K %*% U gives cell-wise loadings per component.
  load_mat <- fit$K %*% fit$U
  cells <- rownames(load_mat) %||% rownames(fit$K) %||% paste0("cell", seq_len(nrow(load_mat)))
  cell_meta <- sdam_cell_table()
  ord <- match(cells, cell_meta$cell)
  df <- data.frame(
    cell      = factor(rep(cells, ncol(load_mat)), levels = cells),
    component = factor(rep(paste0("LV", seq_len(ncol(load_mat))), each = nrow(load_mat)),
                       levels = paste0("LV", seq_len(ncol(load_mat)))),
    task      = factor(rep(cell_meta$task[ord],    ncol(load_mat))),
    measure   = factor(rep(cell_meta$measure[ord], ncol(load_mat))),
    loading   = as.numeric(load_mat)
  )
  ggplot2::ggplot(df, ggplot2::aes(measure, task, fill = loading)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.4) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%+0.2f", loading)),
                       size = 3.2, colour = "grey20") +
    ggplot2::facet_wrap(~ component, nrow = 1L) +
    ggplot2::scale_fill_gradient2(low = "#3B4CC0", mid = "white", high = "#B40426",
                                   midpoint = 0, name = "K %*% U") +
    ggplot2::labs(
      title = "Design saliences in cell space",
      subtitle = "Each panel is one DKGE component; cells laid out as task x measure",
      x = "measure", y = "task"
    ) +
    dkge::theme_dkge()
}

plot_sdam_design_interaction <- function(fit) {
  load_mat <- fit$K %*% fit$U
  cells <- rownames(load_mat) %||% rownames(fit$K)
  cell_meta <- sdam_cell_table()
  ord <- match(cells, cell_meta$cell)
  df <- data.frame(
    task      = factor(rep(cell_meta$task[ord],    ncol(load_mat))),
    measure   = factor(rep(cell_meta$measure[ord], ncol(load_mat))),
    component = factor(rep(paste0("LV", seq_len(ncol(load_mat))), each = nrow(load_mat)),
                       levels = paste0("LV", seq_len(ncol(load_mat)))),
    loading   = as.numeric(load_mat)
  )
  ggplot2::ggplot(df, ggplot2::aes(task, loading,
                                   colour = measure, group = measure)) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_point(size = 2.2) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey60", linetype = 2) +
    ggplot2::facet_wrap(~ component, nrow = 1L) +
    ggplot2::scale_colour_brewer(palette = "Set1") +
    ggplot2::labs(
      title = "Design interaction (task x measure) per component",
      x = "task", y = "Cell loading"
    ) +
    dkge::theme_dkge()
}

sdam_common_component_loadings <- function(fit) {
  Btil <- fit$Btil
  K <- fit$K
  U <- fit$U
  weights <- fit$weights %||% rep(1, length(Btil))
  weights <- pmax(as.numeric(weights), 0)
  if (!sum(weights)) weights[] <- 1

  loadings <- Map(function(Bts, w) {
    w * (t(Bts) %*% K %*% U)
  }, Btil, weights)
  Reduce(`+`, loadings) / sum(weights)
}

sdam_component_cell_scores <- function(fit, manifest) {
  # PLS-style brain scores: every subject-cell map is projected onto the same
  # pooled component loading map. This differs from subject-specific loading
  # maps A_s = t(B_s) K U, which are spatial-expression objects rather than
  # cell-score axes.
  meta <- unique(manifest[, c("subject", "group")])
  Btil <- fit$Btil
  common_loadings <- sdam_common_component_loadings(fit)
  rows <- vector("list", length(Btil) * nrow(Btil[[1]]) * ncol(common_loadings))
  k <- 1L
  for (s in seq_along(Btil)) {
    Bts <- Btil[[s]]
    scores <- Bts %*% common_loadings
    cells <- rownames(Bts) %||% rownames(fit$K)
    cell_meta <- sdam_cell_table()
    ord <- match(cells, cell_meta$cell)
    sid <- fit$subject_ids[[s]] %||% paste0("sub", s)
    for (j in seq_len(ncol(scores))) {
      rows[[k]] <- data.frame(
        subject  = sid,
        cell     = cells,
        task     = cell_meta$task[ord],
        measure  = cell_meta$measure[ord],
        component = paste0("LV", j),
        score    = scores[, j],
        stringsAsFactors = FALSE
      )
      k <- k + 1L
    }
  }
  df <- do.call(rbind, rows)
  df$group <- meta$group[match(df$subject, meta$subject)]
  df$cell <- factor(df$cell, levels = sdam_cell_table()$cell)
  df$task <- factor(df$task, levels = SDAM_TASKS)
  df$measure <- factor(df$measure, levels = SDAM_MEASURES)
  df$component <- factor(df$component, levels = paste0("LV", seq_len(ncol(common_loadings))))
  df
}

sdam_component_score_summary <- function(score_df) {
  agg <- aggregate(score ~ component + group + cell + task + measure,
                   score_df,
                   function(x) c(mean = mean(x), sd = stats::sd(x), n = length(x)))
  data.frame(
    component = agg$component,
    group = agg$group,
    cell = agg$cell,
    task = agg$task,
    measure = agg$measure,
    mean = agg$score[, "mean"],
    sd = agg$score[, "sd"],
    n = agg$score[, "n"],
    row.names = NULL
  )
}

sdam_component_score_tests <- function(score_df,
                                       contrasts = sdam_contrast_matrix()) {
  cells <- rownames(contrasts)
  split_by_subject <- split(score_df, interaction(score_df$component,
                                                  score_df$subject,
                                                  drop = TRUE))
  rows <- lapply(split_by_subject, function(d) {
    d <- d[match(cells, as.character(d$cell)), , drop = FALSE]
    vals <- vapply(seq_len(ncol(contrasts)), function(j) {
      sum(d$score * contrasts[, j])
    }, numeric(1))
    data.frame(
      component = as.character(d$component[[1]]),
      subject = as.character(d$subject[[1]]),
      group = as.character(d$group[[1]]),
      contrast = colnames(contrasts),
      value = vals,
      stringsAsFactors = FALSE
    )
  })
  long <- do.call(rbind, rows)

  out <- lapply(split(long, interaction(long$component, long$contrast, drop = TRUE)), function(d) {
    tt <- stats::t.test(value ~ group, data = d)
    means <- tapply(d$value, d$group, mean)
    data.frame(
      component = d$component[[1]],
      contrast = d$contrast[[1]],
      mean_control = unname(means[["control"]]),
      mean_sdam = unname(means[["sdam"]]),
      diff_sdam_minus_control = unname(means[["sdam"]] - means[["control"]]),
      statistic = unname(tt$statistic),
      p = tt$p.value,
      stringsAsFactors = FALSE
    )
  })
  res <- do.call(rbind, out)
  res$p_adjusted <- stats::p.adjust(res$p, method = "fdr")
  rownames(res) <- NULL
  res[order(res$component, match(res$contrast, colnames(contrasts))), , drop = FALSE]
}

plot_sdam_subject_scores <- function(fit, manifest, comp = 1L) {
  df <- sdam_component_cell_scores(fit, manifest)
  df <- df[df$component == paste0("LV", comp), , drop = FALSE]
  ggplot2::ggplot(df, ggplot2::aes(measure, score, fill = group)) +
    ggplot2::geom_violin(position = ggplot2::position_dodge(0.6),
                         alpha = 0.55, scale = "width") +
    ggplot2::geom_boxplot(position = ggplot2::position_dodge(0.6),
                          width = 0.18, outlier.shape = NA, alpha = 0.85) +
    ggplot2::facet_wrap(~ task) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey60", linetype = 2) +
    ggplot2::labs(
      title    = sprintf("Subject brain scores on LV%d", comp),
      subtitle = "Cell maps projected onto the pooled DKGE component loading",
      x = "measure", y = "Score"
    ) +
    dkge::theme_dkge()
}

sdam_subject_component_projections <- function(fit, manifest,
                                               mode = c("loso", "pooled")) {
  mode <- match.arg(mode)
  meta <- unique(manifest[, c("subject", "group")])
  subject_ids <- fit$subject_ids %||% names(fit$Btil) %||%
    paste0("subject", seq_along(fit$Btil))
  comps <- paste0("LV", seq_len(ncol(fit$U)))

  if (identical(mode, "pooled")) {
    score_df <- sdam_component_cell_scores(fit, manifest)
    saliences <- fit$K %*% fit$U
    rownames(saliences) <- rownames(fit$K) %||% rownames(fit$U)
    rows <- lapply(split(score_df, interaction(score_df$subject,
                                               score_df$component,
                                               drop = TRUE)), function(d) {
      comp <- as.character(d$component[[1L]])
      j <- match(comp, comps)
      d <- d[match(rownames(saliences), as.character(d$cell)), , drop = FALSE]
      data.frame(
        subject = as.character(d$subject[[1L]]),
        group = as.character(d$group[[1L]]),
        component = comp,
        projection = sum(d$score * saliences[, j]),
        mode = mode,
        stringsAsFactors = FALSE
      )
    })
    out <- do.call(rbind, rows)
  } else {
    assignments <- lapply(seq_along(fit$Btil), function(s) s)
    folds <- dkge:::.dkge_build_fold_bases(
      fit,
      assignments = assignments,
      ridge = 0,
      align = TRUE,
      loader_scope = "heldout",
      verbose = FALSE
    )$folds
    rows <- vector("list", length(folds) * ncol(fit$U))
    k <- 1L
    for (fold in folds) {
      s <- fold$subjects[[1L]]
      sid <- subject_ids[[s]]
      R_align <- fold$rotation %||% diag(ncol(fit$U))
      U_aligned <- fold$basis %*% R_align
      saliences <- fit$K %*% U_aligned
      loader <- fold$loaders[[as.character(s)]]
      Y_aligned <- loader$Y %*% R_align
      for (j in seq_len(ncol(fit$U))) {
        rows[[k]] <- data.frame(
          subject = sid,
          group = as.character(meta$group[match(sid, meta$subject)]),
          component = comps[[j]],
          projection = sum(Y_aligned[, j] * saliences[, j]),
          mode = mode,
          stringsAsFactors = FALSE
        )
        k <- k + 1L
      }
    }
    out <- do.call(rbind, rows)
  }

  out$component <- factor(out$component, levels = comps)
  out$group <- factor(out$group, levels = SDAM_GROUPS)
  out
}

plot_sdam_subject_component_projections <- function(fit, manifest,
                                                    mode = c("loso", "pooled")) {
  mode <- match.arg(mode)
  df <- sdam_subject_component_projections(fit, manifest, mode = mode)
  ggplot2::ggplot(df, ggplot2::aes(component, projection, colour = group)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey65", linetype = 2) +
    ggplot2::geom_point(
      position = ggplot2::position_jitter(width = 0.12, height = 0, seed = 1),
      size = 2.3,
      alpha = 0.9
    ) +
    ggplot2::scale_colour_brewer(palette = "Set1") +
    ggplot2::labs(
      title = sprintf("Subject projections on DKGE components (%s)", toupper(mode)),
      subtitle = if (identical(mode, "loso")) {
        "Each point uses a basis learned without that subject, aligned back to the pooled component axes"
      } else {
        "Each point uses the pooled component axes and is descriptive"
      },
      x = "Component",
      y = "Subject projection"
    ) +
    dkge::theme_dkge()
}

plot_sdam_bsr_volume <- function(z, mask, threshold = 3, n_slices = 9L,
                                  title = "BSR-like z-map",
                                  subtitle = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(NULL)
  arr <- as.array(mask) > 0
  vol <- array(NA_real_, dim = dim(mask))
  vol[arr] <- ifelse(abs(z) >= threshold, z, NA_real_)

  # Pick slices that actually contain mask voxels, evenly spaced.
  zs_with_mask <- which(apply(arr, 3, any))
  if (length(zs_with_mask) < n_slices) {
    slices <- zs_with_mask
  } else {
    slices <- zs_with_mask[round(seq(1, length(zs_with_mask), length.out = n_slices))]
  }

  df_list <- lapply(slices, function(zk) {
    sl <- vol[, , zk]
    base_df <- expand.grid(x = seq_len(nrow(sl)), y = seq_len(ncol(sl)))
    base_df$z   <- as.numeric(sl)
    base_df$mask_present <- as.logical(arr[, , zk])
    base_df$slice <- factor(sprintf("z = %d", zk),
                            levels = sprintf("z = %d", slices))
    base_df
  })
  df <- do.call(rbind, df_list)
  finite_z <- df$z[is.finite(df$z)]
  lim <- if (length(finite_z)) max(abs(finite_z), threshold + 0.5) else threshold + 1

  ggplot2::ggplot(df, ggplot2::aes(x, y)) +
    ggplot2::geom_raster(data = subset(df, mask_present),
                         fill = "grey88") +
    ggplot2::geom_raster(data = subset(df, is.finite(z)),
                         ggplot2::aes(fill = z)) +
    ggplot2::scale_fill_gradient2(low = "#3B4CC0", mid = "white", high = "#B40426",
                                  midpoint = 0, limits = c(-lim, lim)) +
    ggplot2::facet_wrap(~ slice, ncol = 3) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = title,
      subtitle = if (is.null(subtitle)) {
        sprintf("|z| >= %g shown; q-space bootstrap of contrast", threshold)
      } else {
        subtitle
      },
      x = NULL, y = NULL, fill = "z"
    ) +
    dkge::theme_dkge() +
    ggplot2::theme(axis.text = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank(),
                   panel.spacing = ggplot2::unit(2, "pt"),
                   strip.text = ggplot2::element_text(size = 9))
}

plot_sdam_bsr_surface <- function(z, mask,
                                   threshold = 3,
                                   surface_atlas = NULL,
                                   views = c("lateral", "medial"),
                                   palette = "vik",
                                   base_color = "#e8e8e8") {
  if (!requireNamespace("neuroatlas", quietly = TRUE)) return(NULL)
  if (is.null(surface_atlas)) {
    surface_atlas <- neuroatlas::schaefer_surf(
      parcels = 200, networks = 7, space = "fsaverage6", surf = "inflated"
    )
  }
  overlay <- sdam_values_to_neurovol(z, mask)
  base_colors <- stats::setNames(
    rep(base_color, length(surface_atlas$ids)),
    as.character(surface_atlas$ids)
  )
  lim <- max(abs(z[is.finite(z)]), threshold + 1, na.rm = TRUE)
  neuroatlas::plot_brain(
    surface_atlas,
    colors = base_colors,
    overlay = overlay,
    overlay_threshold = threshold,
    overlay_palette = palette,
    overlay_lim = c(-lim, lim),
    overlay_alpha = 0.9,
    views = views,
    interactive = FALSE,
    style = "ggseg_like",
    colorbar = TRUE,
    title = "DKGE bootstrap z-map",
    subtitle = sprintf("|z| >= %g; q-space bootstrap, voxel projection", threshold)
  )
}

# -----------------------------------------------------------------------------
# Assembly: run + save
# -----------------------------------------------------------------------------

summarise_sdam_dkge <- function(fit) {
  ve <- dkge::dkge_variance_explained(fit)
  data.frame(
    component        = ve$component,
    sdev             = round(ve$sdev, 4),
    variance         = round(ve$variance, 4),
    variance_percent = round(100 * ve$prop_var, 2),
    cumulative       = round(100 * ve$cum_prop_var, 2)
  )
}

save_sdam_plots <- function(fit, mask, contrasts_obj, group_diff, between, bsr,
                            output_dir, lv = 1L, score_comps = 1:2, threshold = 3) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  ggsave <- function(name, plot, w = 7, h = 4) {
    if (is.null(plot)) return(invisible(NULL))
    ggplot2::ggsave(file.path(output_dir, name), plot, width = w, height = h, dpi = 150)
  }

  ggsave("singular-values.png", dkge::dkge_plot_scree(fit))
  subj_panels <- dkge::dkge_plot_subject_contrib(fit)
  ggsave("subject-weights.png", subj_panels$weights, w = 8, h = 4)
  ggsave("subject-energy.png", subj_panels$energy, w = 6, h = 5)
  ggsave("design-heatmap.png", plot_sdam_design_heatmap(fit), w = 8, h = 3.4)
  ggsave("design-interaction.png", plot_sdam_design_interaction(fit), w = 8, h = 3.4)
  for (comp in score_comps) {
    ggsave(sprintf("lv%d-brain-scores.png", comp),
           plot_sdam_subject_scores(fit, between$participants[, c("subject", "group")],
                                    comp = comp),
           w = 7, h = 4)
  }
  projection_df <- sdam_subject_component_projections(
    fit,
    between$participants[, c("subject", "group")],
    mode = "loso"
  )
  ggsave(
    "subject-component-projections-loso.png",
    plot_sdam_subject_component_projections(
      fit,
      between$participants[, c("subject", "group")],
      mode = "loso"
    ),
    w = 7,
    h = 4
  )
  utils::write.csv(
    projection_df,
    file.path(output_dir, "subject-component-projections-loso.csv"),
    row.names = FALSE
  )
  ggsave(sprintf("lv%d-bsr-volume.png", lv), plot_sdam_bsr_volume(bsr$z, mask, threshold), w = 6, h = 5)

  if (requireNamespace("neuroatlas", quietly = TRUE)) {
    surf <- try(plot_sdam_bsr_surface(bsr$z, mask, threshold = threshold), silent = TRUE)
    if (!inherits(surf, "try-error")) ggsave("bsr-surface.png", surf, w = 9, h = 6)
  }

  utils::write.csv(
    summarise_sdam_dkge(fit),
    file.path(output_dir, "component-summary.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    between$perm$summary[, c("term", "statistic", "p", "p_adjusted")],
    file.path(output_dir, "between-permutation-summary.csv"),
    row.names = FALSE
  )
  invisible(output_dir)
}

sdam_loso_summary <- function(group_maps) {
  rows <- lapply(names(group_maps), function(nm) {
    gm <- group_maps[[nm]]
    data.frame(
      contrast = nm,
      mean_min = round(min(gm$mean), 3),
      mean_max = round(max(gm$mean), 3),
      z_min    = round(min(gm$z), 2),
      z_max    = round(max(gm$z), 2),
      n_z_ge3  = sum(abs(gm$z) >= 3),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

sdam_bsr_summary <- function(bsr_all) {
  rows <- lapply(names(bsr_all), function(nm) {
    z <- bsr_all[[nm]]$z
    data.frame(contrast = nm,
               min_z = round(min(z), 2),
               max_z = round(max(z), 2),
               pos_ge3 = sum(z >= 3),
               neg_le_neg3 = sum(z <= -3),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

sdam_group_diff_summary <- function(group_diff) {
  rows <- lapply(names(group_diff), function(nm) {
    gd <- group_diff[[nm]]
    data.frame(
      contrast = nm,
      n_sdam = gd$n_sdam, n_control = gd$n_control,
      mean_diff_min = round(min(gd$mean_diff), 3),
      mean_diff_max = round(max(gd$mean_diff), 3),
      z_min = round(min(gd$z), 2),
      z_max = round(max(gd$z), 2),
      n_z_ge3 = sum(abs(gd$z) >= 3),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

run_sdam_dkge <- function(root        = sdam_testdata_root(),
                          rank        = 3L,
                          rho_int     = 1.0,
                          nperm       = 999L,
                          nboot       = 200L,
                          max_voxels  = Inf,
                          contrast_name = "task",
                          w_method    = "none",
                          effect_scaling = "none",
                          seed        = SDAM_SEED) {
  manifest <- sdam_time("Build SDAM manifest", build_sdam_manifest(root))
  mask     <- sdam_time(
    "Derive common nonzero mask",
    derive_common_nonzero_mask(manifest, max_voxels = max_voxels)
  )

  loaded   <- sdam_time("Load masked SDAM maps", load_sdam_betas(manifest, mask))
  kernel   <- make_sdam_design_kernel(rho_int = rho_int)

  set.seed(seed)
  fit <- sdam_time(
    "Fit pooled DKGE",
    fit_sdam_dkge(loaded$B_list, kernel, rank = rank,
                  w_method = w_method, effect_scaling = effect_scaling)
  )

  loso_contrasts <- sdam_time(
    "Compute LOSO contrasts",
    compute_sdam_loso_contrasts(fit, contrasts = sdam_contrast_matrix(),
                                method = "loso")
  )
  group_maps  <- sdam_time("Summarise LOSO group maps",
                           sdam_contrast_group_maps(loso_contrasts, manifest))
  group_diff  <- sdam_time("Compute voxel group differences",
                           sdam_contrast_group_diff(loso_contrasts, manifest))
  between     <- sdam_time(
    "Run voxel-level between-subject DKGE",
    run_sdam_between_subjects(fit, manifest,
                              contrast_obj = loso_contrasts,
                              nperm = nperm, seed = seed)
  )
  basis_check <- sdam_time(
    "Run stratified basis diagnostic",
    run_sdam_stratified_basis_check(loaded$B_list, kernel,
                                    manifest, rank = rank)
  )
  group_explicit <- sdam_time(
    "Run q=8 group-explicit DKGE validation",
    run_sdam_group_explicit_validation(loaded$B_list, manifest,
                                       rank = rank, folds = 4L, seed = seed)
  )
  group_explicit_sensitivity <- sdam_time(
    "Run q=8 weighting sensitivity",
    run_sdam_group_explicit_sensitivity(
      loaded$B_list, manifest,
      rank = rank,
      methods = c("none", "mfa_sigma1", "energy"),
      reference = group_explicit
    )
  )
  cell_mean_bridge <- sdam_time(
    "Run cell-mean DKGE/PLS bridge",
    run_sdam_cell_mean_bridge(loaded$B_list, manifest,
                              rank = rank, nperm = nperm,
                              nboot = nboot, seed = seed)
  )
  component_scores <- sdam_time("Compute DKGE component cell scores",
                                sdam_component_cell_scores(fit, manifest))
  component_score_summary <- sdam_component_score_summary(component_scores)
  component_score_tests <- sdam_component_score_tests(component_scores)

  contrasts <- sdam_contrast_matrix()
  bsr_names <- setdiff(colnames(contrasts), "grand_mean")
  bsr_all <- sdam_time("Compute bootstrap-ratio maps", {
    lapply(bsr_names, function(nm) {
      sdam_bsr_voxels(fit, contrasts[, nm], B = nboot, seed = seed)
    })
  })
  names(bsr_all) <- bsr_names
  bsr <- bsr_all[[contrast_name]]

  list(
    manifest = manifest, mask = mask,
    fit = fit, kernel = kernel,
    loso = loso_contrasts, group_maps = group_maps,
    group_diff = group_diff, between = between,
    basis_check = basis_check,
    group_explicit = group_explicit,
    group_explicit_sensitivity = group_explicit_sensitivity,
    cell_mean_bridge = cell_mean_bridge,
    component_scores = component_scores,
    component_score_summary = component_score_summary,
    component_score_tests = component_score_tests,
    bsr = bsr, bsr_all = bsr_all, contrast_used = contrast_name
  )
}

sdam_slim_fit_for_report <- function(fit) {
  keep <- c(
    "sdev", "U", "evals", "R", "K", "Khalf", "Kihalf", "Chat",
    "weights", "Btil", "subject_ids", "effects", "provenance",
    "kernel_info", "eig_values_full", "rank", "cpca", "solver",
    "weight_spec", "w_method", "w_tau", "ridge_input", "rank_requested",
    "effect_scaling", "effective_rank", "rank_reduced", "Chat_sym", "KU",
    "scores_matrix"
  )
  out <- fit[intersect(keep, names(fit))]
  class(out) <- class(fit)
  out
}

sdam_slim_between_for_report <- function(between) {
  feature_blocks <- sub(":[^:]+$", "", between$target$feature_ids)
  target <- list(
    subject_ids = between$target$subject_ids,
    dim = dim(between$target$Y),
    feature_blocks = factor(feature_blocks, levels = unique(feature_blocks)),
    type = between$target$type,
    feature_space = between$target$feature_space,
    crossfit = between$target$crossfit
  )

  perm <- between$perm
  perm$tests <- NULL
  perm$term_maps <- NULL
  perm$call <- NULL
  if (!is.null(perm$feature_tests)) {
    perm$feature_tests <- lapply(perm$feature_tests, function(ft) {
      blocks <- sub(":[^:]+$", "", ft$feature_ids)
      block_levels <- levels(target$feature_blocks)
      if (is.null(block_levels)) block_levels <- unique(blocks)
      ft$feature_blocks <- factor(blocks, levels = block_levels)
      ft$feature_index <- as.integer(ave(seq_along(blocks), blocks, FUN = seq_along))
      ft$feature_ids <- NULL
      ft$term_map <- NULL
      ft$null_max <- NULL
      names(ft$statistic) <- NULL
      names(ft$p) <- NULL
      names(ft$p_adjusted) <- NULL
      ft
    })
  }

  focused <- NULL
  if (!is.null(between$focused)) {
    fperm <- between$focused$perm
    fperm$tests <- NULL
    fperm$term_maps <- NULL
    fperm$call <- NULL
    if (!is.null(fperm$feature_tests)) {
      fperm$feature_tests <- lapply(fperm$feature_tests, function(ft) {
        ft$term_map <- NULL
        ft$null_max <- NULL
        ft$feature_ids <- NULL
        names(ft$statistic) <- NULL
        names(ft$p) <- NULL
        names(ft$p_adjusted) <- NULL
        ft
      })
    }
    focused <- list(
      contrast_block = between$focused$contrast_block,
      target_dim = c(
        subjects = length(between$focused$fit$target$subject_ids),
        features = sum(between$focused$feature_mask)
      ),
      perm = fperm
    )
  }

  list(
    target = target,
    design = list(
      term_names = between$design$term_names,
      subject_ids = between$design$subject_ids,
      formula = deparse(between$design$formula)
    ),
    perm = perm,
    focused = focused,
    participants = between$participants
  )
}

sdam_slim_cell_mean_bridge_for_report <- function(bridge) {
  perm_align <- bridge$permutation$alignment_summary %||% list()
  boot_align <- bridge$bootstrap$alignment_summary %||% list()

  list(
    estimand = bridge$estimand,
    target_dim = dim(bridge$target$Y),
    row_metadata = bridge$target$row_metadata,
    aggregate = bridge$target$aggregate,
    inference_unit = bridge$target$provenance$inference_unit,
    center = bridge$center,
    rank = bridge$rank,
    component = bridge$component,
    singular_values = bridge$fit$singular_values,
    saliences = sdam_aggregate_salience_table(bridge$fit),
    contrast = bridge$contrast,
    component_contrasts = bridge$component_contrasts,
    statistic = bridge$statistic,
    permutation = list(
      observed = bridge$permutation$observed,
      p = bridge$permutation$p,
      B = bridge$permutation$B,
      null = bridge$permutation$null,
      resampling = bridge$permutation$resampling,
      alignment = perm_align,
      n_near_tie = perm_align$n_near_tie %||% NA_integer_,
      n_weak_alignment = perm_align$n_weak %||% NA_integer_,
      min_alignment_cosine = perm_align$min_cosine %||% NA_real_
    ),
    bootstrap = list(
      observed = bridge$bootstrap$observed,
      interval = bridge$bootstrap$interval,
      conf = bridge$bootstrap$conf,
      B = bridge$bootstrap$B,
      statistics = bridge$bootstrap$statistics,
      component_contrasts = bridge$bootstrap$component_contrasts,
      feature_bootstrap = bridge$bootstrap$feature_bootstrap,
      alignment = boot_align,
      resampling = bridge$bootstrap$resampling
    )
  )
}

sdam_slim_group_explicit_for_report <- function(group_explicit) {
  validation_summary <- group_explicit$validation$summary
  numeric_cols <- vapply(validation_summary, is.numeric, logical(1))
  validation_summary[numeric_cols] <- lapply(validation_summary[numeric_cols], unname)

  list(
    component_summary = group_explicit$component_summary,
    design_loadings = sdam_q8_design_loading_table(group_explicit$fit),
    validation_summary = validation_summary,
    contrast_estimability = group_explicit$validation$contrast_estimability,
    target_dim = c(
      subjects = length(group_explicit$data$subject_ids),
      effects = group_explicit$data$q
    ),
    cell_labels = group_explicit$cell_labels,
    term_scope = group_explicit$term_scope,
    block_factors = group_explicit$block_factors,
    coverage = group_explicit$coverage,
    pair_counts = group_explicit$pair_counts
  )
}

sdam_slim_q8_sensitivity_for_report <- function(sensitivity) {
  lapply(sensitivity, function(x) {
    list(
      w_method = x$w_method,
      component_summary = x$component_summary,
      score_values = x$score_values,
      score_summary = x$score_summary,
      score_tests = x$score_tests,
      weight_table = x$weight_table,
      weight_summary = x$weight_summary
    )
  })
}

sdam_slim_analysis_for_report <- function(analysis) {
  bsr_z <- lapply(analysis$bsr_all, function(x) {
    list(z = unname(x$z), B = x$B)
  })
  list(
    manifest = analysis$manifest,
    mask = analysis$mask,
    kernel = analysis$kernel,
    component_summary = summarise_sdam_dkge(analysis$fit),
    component_score_summary = analysis$component_score_summary,
    component_score_tests = analysis$component_score_tests,
    loso_summary = sdam_loso_summary(analysis$group_maps),
    bsr_summary = sdam_bsr_summary(analysis$bsr_all),
    group_diff_summary = sdam_group_diff_summary(analysis$group_diff),
    between = sdam_slim_between_for_report(analysis$between),
    basis_check = list(
      angles = analysis$basis_check$angles,
      rank_strat = analysis$basis_check$rank_strat,
      rank_requested = analysis$basis_check$rank_requested,
      q = analysis$basis_check$q,
      forced_overlap = analysis$basis_check$forced_overlap
    ),
    group_explicit = sdam_slim_group_explicit_for_report(analysis$group_explicit),
    group_explicit_sensitivity = sdam_slim_q8_sensitivity_for_report(
      analysis$group_explicit_sensitivity
    ),
    cell_mean_bridge = sdam_slim_cell_mean_bridge_for_report(analysis$cell_mean_bridge),
    bsr_all = bsr_z,
    contrast_used = analysis$contrast_used
  )
}

sdam_save_rds <- function(object, file, compress = "gzip") {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(file, ".tmp")
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(object, tmp, compress = compress)
  if (file.exists(file) && !file.remove(file)) {
    stop("Could not replace existing file: ", file, call. = FALSE)
  }
  if (!file.rename(tmp, file)) {
    stop("Could not move temporary RDS into place: ", file, call. = FALSE)
  }
  invisible(file)
}

sdam_write_q8_sensitivity_tables <- function(sensitivity, output_dir) {
  bind_with_method <- function(field) {
    do.call(rbind, lapply(names(sensitivity), function(method) {
      tab <- sensitivity[[method]][[field]]
      if (is.null(tab)) return(NULL)
      tab$w_method <- method
      tab
    }))
  }
  utils::write.csv(
    bind_with_method("score_tests"),
    file.path(output_dir, "q8-component-score-tests-by-weighting.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    bind_with_method("weight_table"),
    file.path(output_dir, "q8-subject-weights-by-weighting.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    bind_with_method("weight_summary"),
    file.path(output_dir, "q8-subject-weight-summary-by-group.csv"),
    row.names = FALSE
  )
  invisible(output_dir)
}

main <- function() {
  output_dir <- file.path("artifacts", "sdam-firstlevel-task-dkge")
  nperm      <- as.integer(Sys.getenv("DKGE_SDAM_NPERM",  "999"))
  nboot      <- as.integer(Sys.getenv("DKGE_SDAM_NBOOT",  "200"))
  max_voxels <- {
    raw <- Sys.getenv("DKGE_SDAM_MAX_VOXELS", "")
    if (!nzchar(raw)) Inf else as.numeric(raw)
  }

  analysis <- run_sdam_dkge(
    nperm = nperm, nboot = nboot, max_voxels = max_voxels
  )

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  analysis_report <- sdam_slim_analysis_for_report(analysis)
  sdam_time("Save report analysis bundle",
            sdam_save_rds(analysis_report, file.path(output_dir, "analysis.rds")))
  sdam_write_q8_sensitivity_tables(analysis_report$group_explicit_sensitivity,
                                   output_dir)

  save_components <- tolower(Sys.getenv("DKGE_SDAM_SAVE_COMPONENTS", "false")) %in%
    c("1", "true", "yes")
  if (save_components) {
    sdam_time("Save compact between-subject result",
              sdam_save_rds(analysis_report$between,
                            file.path(output_dir, "dkge-between-subjects.rds")))
  }

  save_sdam_plots(analysis$fit, analysis$mask, analysis$loso,
                  analysis$group_diff, analysis$between, analysis$bsr,
                  output_dir = output_dir)

  message("\n--- Component summary ---")
  print(summarise_sdam_dkge(analysis$fit))
  message("\n--- Between-subject permutation (~ group) ---")
  print(analysis$between$perm$summary[, c("term", "statistic", "p", "p_adjusted")])
  message("\n--- q=8 group-explicit validation ---")
  print(analysis$group_explicit$validation$summary[, c("contrast", "estimability",
                                                       "estimate_observed",
                                                       "estimate_completed",
                                                       "sensitivity")])
  message("\n--- q=8 LV2 task:measure score tests by weighting ---")
  lv2_tests <- do.call(rbind, lapply(analysis_report$group_explicit_sensitivity, function(x) {
    tab <- x$score_tests
    tab <- tab[tab$component == "LV2" & tab$contrast == "task:measure", , drop = FALSE]
    tab$w_method <- x$w_method
    tab
  }))
  print(lv2_tests[, c("w_method", "mean_control", "mean_sdam",
                      "diff_sdam_minus_control", "statistic", "p", "p_adjusted")])
  message("\n--- Cell-mean DKGE/PLS bridge permutation ---")
  print(data.frame(
    statistic = analysis_report$cell_mean_bridge$permutation$observed,
    p = analysis_report$cell_mean_bridge$permutation$p,
    B = analysis_report$cell_mean_bridge$permutation$B,
    bootstrap_low = analysis_report$cell_mean_bridge$bootstrap$interval[[1]],
    bootstrap_high = analysis_report$cell_mean_bridge$bootstrap$interval[[2]]
  ))
  invisible(analysis)
}

if (sys.nframe() == 0L) {
  suppressPackageStartupMessages({
    library(dkge)
    library(ggplot2)
  })
  main()
}
