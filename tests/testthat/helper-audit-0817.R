# Shared fixtures for the 2026-08-17 new-paths remediation (phase 07).
# Used by A1 (aggregate permutation), B1/D1 (permuted cell kernel), and
# B5 (partial-coverage cohort mass). Numeric oracles live in the calling
# tests; these helpers only build the data.

# Fixed derangement of the 2 x 3 cell grid (task x measure). No label stays
# in kernel order, so positional kernel_info matching is wrong on every row.
.audit0817_cell_perm <- c(3L, 6L, 1L, 5L, 2L, 4L)

#' 2 x 3 cell-kernel fit whose data effects are a fixed non-identity permutation
#'
#' Builds `design_kernel(basis = "cell")` for `task` (2) x `measure` (3) and
#' subjects whose beta rownames / design colnames are
#' `grid$cell_labels[.audit0817_cell_perm]`. The kernel itself stays in grid
#' order, so `.dkge_align_kernel_effects()` must permute `K` (and, after B1,
#' `kernel_info`) to match `fit$effects`.
#'
#' @return List with `fit`, `perm`, `labels` (kernel order), `permuted_labels`,
#'   `factors`, and `kernel`.
make_permuted_cell_fit <- function(seed = 1L, S = 4L, P = 8L, rank = 3L) {
  set.seed(seed)
  factors <- list(task = list(L = 2L), measure = list(L = 3L))
  kernel <- design_kernel(factors, basis = "cell", normalize = "none")
  labels <- kernel$info$cell_labels
  q <- length(labels)
  perm <- .audit0817_cell_perm
  if (length(perm) != q || identical(perm, seq_len(q))) {
    stop("Internal cell permutation must be a derangement of the 2 x 3 grid.",
         call. = FALSE)
  }
  permuted_labels <- labels[perm]

  betas <- replicate(S, {
    B <- matrix(stats::rnorm(q * P), q, P)
    rownames(B) <- permuted_labels
    B
  }, simplify = FALSE)
  designs <- replicate(S, {
    X <- diag(q)
    dimnames(X) <- list(permuted_labels, permuted_labels)
    X
  }, simplify = FALSE)

  fit <- suppressMessages(dkge(betas, designs, K = kernel, rank = rank))
  list(
    fit = fit,
    perm = perm,
    labels = labels,
    permuted_labels = permuted_labels,
    factors = factors,
    kernel = kernel
  )
}

#' Partial-coverage fit: subject 1 sees all q rows, subjects 2..S see 1:(q-1)
#'
#' Default geometry is S = 6, q = 4 so the unobserved e4 pair is seen by 1 of
#' 6 subjects. Unit `effect_n` is attached so `dkge_effect_weights("count")`
#' is valid. Extra arguments are forwarded to [dkge_fit()].
#'
#' @return List with `fit`, raw `B_list`, `data`, and `effects`.
make_partial_fit <- function(seed = 1L,
                             w_method = "none",
                             S = 6L,
                             q = 4L,
                             P = 8L,
                             rank = 2L,
                             ...) {
  set.seed(seed)
  effects <- paste0("e", seq_len(q))
  B_list <- lapply(seq_len(S), function(s) {
    B <- matrix(stats::rnorm(q * P), q, P,
                dimnames = list(effects, paste0("v", seq_len(P))))
    if (s > 1L) {
      B[q, ] <- 0
    }
    B
  })
  subjects <- lapply(seq_len(S), function(s) {
    rows <- if (s == 1L) seq_len(q) else seq_len(q - 1L)
    X <- diag(q)
    colnames(X) <- effects
    suppressWarnings(dkge_subject(
      B_list[[s]], X,
      id = paste0("s", s),
      observed_rows = rows,
      effect_n = stats::setNames(rep(1, q), effects)
    ))
  })
  data <- dkge_data(subjects)
  K <- diag(q)
  dimnames(K) <- list(effects, effects)
  fit <- dkge_fit(data, K = K, rank = rank, w_method = w_method, ...)
  list(fit = fit, B_list = B_list, data = data, effects = effects)
}

#' Pure-null aggregate target with a `grp` factor for permutation
#'
#' n subjects, q cells, p features; values are i.i.d. Gaussian with no
#' group-modulated signal. `group_vars = "grp"` so [dkge_aggregate_permute()]
#' can shuffle labels. Also returns a row-aligned identity `K` and a
#' group-difference `contrast` for the signed built-in statistics.
#'
#' @return List with `target`, `K`, `contrast`, `subject_data`, `cell_data`,
#'   and `values`.
null_aggregate_target <- function(seed = 1L, n = 16L, q = 4L, p = 25L) {
  set.seed(seed)
  subject_data <- data.frame(
    subject_id = paste0("s", seq_len(n)),
    grp = factor(rep(c("A", "B"), length.out = n)),
    stringsAsFactors = FALSE
  )
  cell_data <- data.frame(
    cell = paste0("c", seq_len(q)),
    stringsAsFactors = FALSE
  )
  values <- lapply(seq_len(n), function(i) {
    matrix(stats::rnorm(q * p), q, p,
           dimnames = list(cell_data$cell, paste0("f", seq_len(p))))
  })
  names(values) <- subject_data$subject_id

  target <- dkge_aggregate_target(
    values,
    subject_data = subject_data,
    cell_data = cell_data,
    group_vars = "grp",
    cell_vars = "cell"
  )
  K <- diag(nrow(target$Y))
  dimnames(K) <- list(rownames(target$Y), rownames(target$Y))
  contrast <- ifelse(as.character(target$row_metadata$grp) == "A", 1, -1)
  names(contrast) <- rownames(target$Y)
  list(
    target = target,
    K = K,
    contrast = contrast,
    subject_data = subject_data,
    cell_data = cell_data,
    values = values
  )
}

#' Rejection rate of a seeded p-value function at level `alpha`
#'
#' `fun(seed)` must return a single finite p-value. Used only in
#' `skip_on_cran()` calibration tests with a fixed `seeds` vector.
type1_rate <- function(fun, n_rep, alpha = 0.05, seeds = NULL) {
  n_rep <- as.integer(n_rep)
  if (length(n_rep) != 1L || is.na(n_rep) || n_rep < 1L) {
    stop("`n_rep` must be a positive integer.", call. = FALSE)
  }
  if (is.null(seeds)) {
    seeds <- seq_len(n_rep)
  } else {
    seeds <- as.integer(seeds)
    if (length(seeds) != n_rep) {
      stop("`seeds` must have length `n_rep`.", call. = FALSE)
    }
  }
  pvals <- vapply(seeds, function(s) {
    p <- fun(s)
    if (!is.numeric(p) || length(p) != 1L || !is.finite(p)) {
      stop("`fun` must return a single finite p-value.", call. = FALSE)
    }
    as.numeric(p)
  }, numeric(1))
  mean(pvals < alpha)
}
