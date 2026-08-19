# dkge-aggregate.R
# General aggregate/cell-mean targets and subject-resampling bridge inference.

#' Build an aggregate target from subject-level repeated-measures data
#'
#' Constructs row-by-feature aggregate targets such as
#' `group:task:measure` cell means from subject-level repeated-measures
#' matrices. The returned object stores the source metadata required to
#' recompute aggregates under subject-label permutations or subject bootstrap
#' resampling. The unit of inference remains subjects; aggregate rows are never
#' resampled directly.
#'
#' @param values Named list of subject matrices. Each matrix is cells x
#'   features. When every subject matrix has rownames, rows are matched by
#'   those names (an error if the name sets differ or only some subjects
#'   have rownames). When subject 1 has colnames, every other subject must
#'   have the same column-name set; columns are reordered to subject 1.
#' @param subject_data Data frame with one row per subject.
#' @param cell_data Data frame with one row per repeated-measures cell. If
#'   `NULL`, a cell factor named `cell` is inferred from matrix row names
#'   (or `cell1`, `cell2`, ...) and is kept in the aggregate rows. Pass
#'   `cell_vars = character(0)` to collapse those cells.
#' @param group_vars Character vector of subject-level grouping variables.
#' @param cell_vars Character vector of cell-level variables. Defaults to
#'   the inferred `cell` column when `cell_data` is `NULL`, otherwise to all
#'   columns in `cell_data` except `cell_id_col`.
#' @param subject_id_col Subject identifier column in `subject_data`.
#' @param cell_id_col Optional cell identifier column in `cell_data`.
#' @param weights Optional subject or subject-by-cell weights.
#' @param aggregate Aggregation method. `"mean"` ignores `weights`;
#'   `"weighted_mean"` uses them.
#' @param row_sep Separator used to build aggregate row labels. No level of a
#'   `group_vars`/`cell_vars` variable may contain it, otherwise two distinct
#'   aggregate rows could be given the same row ID; violations are rejected.
#'
#' @return Object of class `dkge_aggregate_target`.
#' @export
dkge_aggregate_target <- function(values,
                                  subject_data,
                                  cell_data = NULL,
                                  group_vars = NULL,
                                  cell_vars = NULL,
                                  subject_id_col = "subject_id",
                                  cell_id_col = NULL,
                                  weights = NULL,
                                  aggregate = c("mean", "weighted_mean"),
                                  row_sep = ":") {
  aggregate <- match.arg(aggregate)
  spec <- .dkge_aggregate_prepare_spec(
    values = values,
    subject_data = subject_data,
    cell_data = cell_data,
    group_vars = group_vars,
    cell_vars = cell_vars,
    subject_id_col = subject_id_col,
    cell_id_col = cell_id_col,
    weights = weights,
    aggregate = aggregate,
    row_sep = row_sep
  )
  .dkge_aggregate_from_spec(spec)
}

.dkge_aggregate_prepare_spec <- function(values,
                                         subject_data,
                                         cell_data,
                                         group_vars,
                                         cell_vars,
                                         subject_id_col,
                                         cell_id_col,
                                         weights,
                                         aggregate,
                                         row_sep) {
  if (!is.list(values) || !length(values)) {
    stop("`values` must be a non-empty named list of subject matrices.",
         call. = FALSE)
  }
  values <- lapply(values, as.matrix)
  if (any(!vapply(values, is.numeric, logical(1)))) {
    stop("Each entry of `values` must be numeric.", call. = FALSE)
  }
  if (any(vapply(values, function(x) any(!is.finite(x)), logical(1)))) {
    stop("`values` contains non-finite entries.", call. = FALSE)
  }

  dims <- vapply(values, dim, integer(2))
  if (length(unique(dims[1, ])) != 1L || length(unique(dims[2, ])) != 1L) {
    stop("All subject matrices must have the same cell and feature dimensions.",
         call. = FALSE)
  }
  q <- unname(dims[1, 1])
  p <- unname(dims[2, 1])

  subject_data <- as.data.frame(subject_data, stringsAsFactors = FALSE)
  if (!subject_id_col %in% names(subject_data)) {
    stop("`subject_id_col` is not present in `subject_data`.", call. = FALSE)
  }
  subject_ids <- as.character(subject_data[[subject_id_col]])
  if (length(subject_ids) != nrow(subject_data) ||
      any(!nzchar(subject_ids)) ||
      any(duplicated(subject_ids))) {
    stop("Subject identifiers must be unique and non-empty.", call. = FALSE)
  }

  value_ids <- names(values)
  if (is.null(value_ids) || any(!nzchar(value_ids))) {
    if (length(values) != length(subject_ids)) {
      stop("Unnamed `values` must be in one-to-one subject_data order.",
           call. = FALSE)
    }
    names(values) <- subject_ids
  } else {
    ord <- match(subject_ids, value_ids)
    if (anyNA(ord)) {
      stop("Every subject in `subject_data` must have a matching values entry.",
           call. = FALSE)
    }
    values <- values[ord]
  }

  cell_data_inferred <- is.null(cell_data)
  if (cell_data_inferred) {
    cell_ids <- rownames(values[[1]]) %||% paste0("cell", seq_len(q))
    cell_data <- data.frame(cell = cell_ids, stringsAsFactors = FALSE)
    cell_id_col <- "cell"
  } else {
    cell_data <- as.data.frame(cell_data, stringsAsFactors = FALSE)
    if (nrow(cell_data) != q) {
      stop("`cell_data` must have one row per matrix row.", call. = FALSE)
    }
  }
  if (!is.null(cell_id_col) && !cell_id_col %in% names(cell_data)) {
    stop("`cell_id_col` is not present in `cell_data`.", call. = FALSE)
  }

  has_rn <- vapply(values, function(x) !is.null(rownames(x)), logical(1))
  if (any(has_rn) && !all(has_rn)) {
    stop("Subject matrix row names must be present for every subject or for none.",
         call. = FALSE)
  }
  if (!is.null(cell_id_col)) {
    aligned <- .dkge_aggregate_align_cell_rows(values, cell_data, cell_id_col)
    values <- aligned$values
    cell_data <- aligned$cell_data
  } else if (all(has_rn)) {
    values <- .dkge_aggregate_align_rows_by_rownames(values, subject_ids)
  }
  values <- .dkge_aggregate_align_feature_cols(values, subject_ids)

  if (is.null(cell_vars)) {
    cell_vars <- if (cell_data_inferred) {
      "cell"
    } else {
      setdiff(names(cell_data), cell_id_col %||% character(0))
    }
  }
  group_vars <- as.character(group_vars %||% character(0))
  cell_vars <- as.character(cell_vars %||% character(0))
  missing_group <- setdiff(group_vars, names(subject_data))
  missing_cell <- setdiff(cell_vars, names(cell_data))
  if (length(missing_group)) {
    stop("Unknown subject grouping variable(s): ",
         paste(missing_group, collapse = ", "), call. = FALSE)
  }
  if (length(missing_cell)) {
    stop("Unknown cell variable(s): ", paste(missing_cell, collapse = ", "),
         call. = FALSE)
  }
  if (!length(c(group_vars, cell_vars))) {
    stop("At least one `group_vars` or `cell_vars` entry is required.",
         call. = FALSE)
  }
  for (v in group_vars) {
    if (anyNA(subject_data[[v]])) {
      stop(sprintf("Subject grouping variable '%s' contains missing values.", v),
           call. = FALSE)
    }
  }
  for (v in cell_vars) {
    if (anyNA(cell_data[[v]])) {
      stop(sprintf("Cell variable '%s' contains missing values.", v),
           call. = FALSE)
    }
  }

  row_sep <- .dkge_aggregate_validate_row_sep(row_sep)
  # Row IDs are the row-defining labels pasted with `row_sep`, so a label that
  # itself contains the separator makes two distinct aggregate rows share an ID
  # (`a:b` + `c` and `a` + `b:c`). Reject the ambiguity at construction rather
  # than shipping a target whose kernel row matching is silently wrong.
  .dkge_aggregate_check_label_separator(subject_data, group_vars, row_sep,
                                        "Subject grouping variable")
  .dkge_aggregate_check_label_separator(cell_data, cell_vars, row_sep,
                                        "Cell variable")

  feature_ids <- colnames(values[[1]]) %||% paste0("feature", seq_len(p))
  feature_ids <- as.character(feature_ids)
  if (length(feature_ids) != p || any(!nzchar(feature_ids)) ||
      any(duplicated(feature_ids))) {
    stop("Feature identifiers must be unique and non-empty.", call. = FALSE)
  }
  values <- lapply(values, function(x) {
    colnames(x) <- feature_ids
    x
  })

  weights <- .dkge_aggregate_validate_weights(weights, subject_ids, q)

  # Stack the subject matrices exactly once. Row (s - 1) * q + c of `V` holds
  # cell `c` of source subject `s`; every resample is then a row-index gather
  # plus a single `rowsum()` rather than a per-cell data-frame scan.
  V <- do.call(rbind, values)
  dimnames(V) <- NULL
  weight_flat <- if (is.null(weights)) {
    NULL
  } else if (is.matrix(weights)) {
    as.numeric(t(weights))
  } else {
    rep(as.numeric(weights), each = q)
  }

  list(
    V = V,
    weight_flat = weight_flat,
    subject_data = subject_data,
    cell_data = cell_data,
    group_vars = group_vars,
    cell_vars = cell_vars,
    subject_id_col = subject_id_col,
    cell_id_col = cell_id_col,
    weights = weights,
    aggregate = aggregate,
    row_sep = row_sep,
    subject_ids = subject_ids,
    feature_ids = feature_ids,
    q = q,
    p = p
  )
}

.dkge_aggregate_validate_row_sep <- function(row_sep) {
  if (!is.character(row_sep) || length(row_sep) != 1L || is.na(row_sep) ||
      !nzchar(row_sep)) {
    stop("`row_sep` must be a single non-empty string.", call. = FALSE)
  }
  row_sep
}

.dkge_aggregate_check_label_separator <- function(data, vars, row_sep, what) {
  for (v in vars) {
    vals <- as.character(data[[v]])
    hit <- grepl(row_sep, vals, fixed = TRUE)
    if (any(hit)) {
      stop(sprintf(paste0("%s '%s' has level(s) containing the row separator ",
                          "'%s': %s. Aggregate row IDs would be ambiguous; ",
                          "rename the level(s) or pick a different `row_sep`."),
                   what, v, row_sep,
                   paste(sprintf("'%s'", unique(vals[hit])), collapse = ", ")),
           call. = FALSE)
    }
  }
  invisible(NULL)
}

.dkge_aggregate_align_cell_rows <- function(values, cell_data, cell_id_col) {
  cell_ids <- as.character(cell_data[[cell_id_col]])
  if (length(cell_ids) != nrow(cell_data) ||
      any(!nzchar(cell_ids)) ||
      any(duplicated(cell_ids))) {
    stop("`cell_id_col` must identify unique non-empty cells.", call. = FALSE)
  }

  values <- lapply(values, function(x) {
    rn <- rownames(x)
    if (is.null(rn)) {
      rownames(x) <- cell_ids
      return(x)
    }
    rn <- as.character(rn)
    if (length(rn) != length(cell_ids) ||
        any(!nzchar(rn)) ||
        any(duplicated(rn))) {
      stop("Subject matrix row names must be unique non-empty cell IDs when `cell_id_col` is provided.",
           call. = FALSE)
    }
    ord <- match(cell_ids, rn)
    if (anyNA(ord)) {
      stop("Subject matrix row names must contain all `cell_id_col` values.",
           call. = FALSE)
    }
    x[ord, , drop = FALSE]
  })
  list(values = values, cell_data = cell_data)
}

.dkge_aggregate_align_rows_by_rownames <- function(values, subject_ids) {
  ref <- as.character(rownames(values[[1]]))
  if (any(!nzchar(ref)) || anyDuplicated(ref)) {
    stop("Subject matrix row names must be unique and non-empty.",
         call. = FALSE)
  }
  lapply(seq_along(values), function(s) {
    x <- values[[s]]
    rn <- as.character(rownames(x))
    if (any(!nzchar(rn)) || anyDuplicated(rn)) {
      stop(sprintf("Row names of subject `%s` must be unique and non-empty.",
                   subject_ids[[s]]), call. = FALSE)
    }
    missing <- setdiff(ref, rn)
    extra <- setdiff(rn, ref)
    if (length(missing) || length(extra)) {
      stop(sprintf(
        "Row names of subject `%s` do not match subject 1:%s%s",
        subject_ids[[s]],
        if (length(missing)) sprintf(" missing %s.", paste(missing, collapse = ", ")) else "",
        if (length(extra)) sprintf(" extra %s.", paste(extra, collapse = ", ")) else ""
      ), call. = FALSE)
    }
    x[match(ref, rn), , drop = FALSE]
  })
}

.dkge_aggregate_align_feature_cols <- function(values, subject_ids) {
  cn1 <- colnames(values[[1]])
  if (is.null(cn1)) {
    return(values)
  }
  cn1 <- as.character(cn1)
  lapply(seq_along(values), function(s) {
    x <- values[[s]]
    if (s == 1L) {
      return(x)
    }
    cn <- colnames(x)
    if (is.null(cn)) {
      stop(sprintf(
        "feature columns of subject `%s` do not match subject 1: missing %s.",
        subject_ids[[s]], paste(cn1, collapse = ", ")
      ), call. = FALSE)
    }
    cn <- as.character(cn)
    missing <- setdiff(cn1, cn)
    extra <- setdiff(cn, cn1)
    if (length(missing) || length(extra)) {
      stop(sprintf(
        "feature columns of subject `%s` do not match subject 1:%s%s",
        subject_ids[[s]],
        if (length(missing)) sprintf(" missing %s.", paste(missing, collapse = ", ")) else "",
        if (length(extra)) sprintf(" extra %s.", paste(extra, collapse = ", ")) else ""
      ), call. = FALSE)
    }
    x[, match(cn1, cn), drop = FALSE]
  })
}

.dkge_aggregate_validate_weights <- function(weights, subject_ids, q) {
  if (is.null(weights)) {
    return(NULL)
  }
  if (is.matrix(weights) || is.data.frame(weights)) {
    weights <- as.matrix(weights)
    if (!is.numeric(weights) || any(!is.finite(weights)) ||
        any(weights < 0) || ncol(weights) != q) {
      stop("Matrix `weights` must be non-negative, finite, and have one column per cell.",
           call. = FALSE)
    }
    if (is.null(rownames(weights))) {
      if (nrow(weights) != length(subject_ids)) {
        stop("Unnamed matrix `weights` must have one row per subject.",
             call. = FALSE)
      }
      rownames(weights) <- subject_ids
    }
    ord <- match(subject_ids, rownames(weights))
    if (anyNA(ord)) {
      stop("Matrix `weights` row names must contain all subject IDs.",
           call. = FALSE)
    }
    return(weights[ord, , drop = FALSE])
  }

  weights <- as.numeric(weights)
  if (length(weights) != length(subject_ids) ||
      any(!is.finite(weights)) ||
      any(weights < 0)) {
    stop("Vector `weights` must be non-negative, finite, and match subjects.",
         call. = FALSE)
  }
  weights
}

.dkge_aggregate_from_spec <- function(spec,
                                      subject_data = spec$subject_data,
                                      subject_indices = seq_along(spec$subject_ids)) {
  n_src <- length(spec$subject_ids)
  subject_indices <- as.integer(subject_indices)
  if (!length(subject_indices) || anyNA(subject_indices) ||
      any(subject_indices < 1L) ||
      any(subject_indices > n_src)) {
    stop("`subject_indices` must select source subjects.", call. = FALSE)
  }
  if (nrow(subject_data) != n_src) {
    stop("`subject_data` must have one row per source subject.", call. = FALSE)
  }
  q <- spec$q
  n_draw <- length(subject_indices)
  row_vars <- c(spec$group_vars, spec$cell_vars)

  source_index <- rep(subject_indices, each = q)
  cell_index <- rep.int(seq_len(q), n_draw)
  value_row <- (source_index - 1L) * q + cell_index
  weight <- if (is.null(spec$weight_flat)) {
    rep(1, length(value_row))
  } else {
    spec$weight_flat[value_row]
  }

  map <- data.frame(
    subject_id = spec$subject_ids[source_index],
    draw_index = rep(seq_len(n_draw), each = q),
    source_index = source_index,
    cell_index = cell_index,
    aggregate_row = NA_integer_,
    weight = weight,
    stringsAsFactors = FALSE
  )
  if (length(spec$group_vars)) {
    map[spec$group_vars] <- subject_data[source_index, spec$group_vars,
                                         drop = FALSE]
  }
  if (length(spec$cell_vars)) {
    map[spec$cell_vars] <- spec$cell_data[cell_index, spec$cell_vars,
                                          drop = FALSE]
  }

  # Row keys are the character representation of the row variables (matching
  # the historical `as.character(...) == as.character(...)` comparison). Codes
  # are pasted as integers separated by "." so the joined key is injective.
  key <- .dkge_aggregate_row_keys(map, row_vars)
  first <- !duplicated(key)
  row_meta <- map[first, row_vars, drop = FALSE]
  ord <- if (length(row_vars)) {
    do.call(order, unname(as.list(row_meta)))
  } else {
    seq_len(nrow(row_meta))
  }
  row_meta <- row_meta[ord, , drop = FALSE]
  rownames(row_meta) <- NULL
  n_rows <- nrow(row_meta)
  group <- match(key, key[first][ord])
  map$aggregate_row <- group
  row_ids <- .dkge_aggregate_row_ids(row_meta, row_vars, spec$row_sep)

  V <- if (identical(value_row, seq_len(nrow(spec$V)))) {
    spec$V
  } else {
    spec$V[value_row, , drop = FALSE]
  }
  if (identical(spec$aggregate, "weighted_mean")) {
    num <- rowsum(V * weight, group = group, reorder = TRUE)
    den <- as.numeric(rowsum(as.matrix(weight), group = group, reorder = TRUE))
  } else {
    num <- rowsum(V, group = group, reorder = TRUE)
    den <- as.numeric(tabulate(group, nbins = n_rows))
  }
  reorder <- match(as.character(seq_len(n_rows)), rownames(num))
  num <- num[reorder, , drop = FALSE]
  if (identical(spec$aggregate, "weighted_mean")) {
    den <- den[reorder]
  }
  if (any(!is.finite(den)) || any(den <= 0)) {
    stop("Aggregate row has zero total weight.", call. = FALSE)
  }
  Y <- num / den
  dimnames(Y) <- list(row_ids, spec$feature_ids)

  out <- list(
    Y = Y,
    row_metadata = row_meta,
    row_ids = row_ids,
    subject_ids = spec$subject_ids[subject_indices],
    subject_indices = subject_indices,
    source_subject_ids = spec$subject_ids,
    feature_ids = spec$feature_ids,
    subject_to_row = map[, c("subject_id", "draw_index", "source_index",
                             "cell_index", "aggregate_row", "weight", row_vars),
                         drop = FALSE],
    weights = spec$weights,
    aggregate = spec$aggregate,
    group_vars = spec$group_vars,
    cell_vars = spec$cell_vars,
    resample_spec = spec,
    provenance = list(
      estimand = "aggregate_cell_mean",
      inference_unit = "subject",
      row_vars = row_vars
    )
  )
  class(out) <- c("dkge_aggregate_target", "list")
  out
}

# Injective character key per source row, built from integer level codes so
# that pasting can never collide across columns.
.dkge_aggregate_row_keys <- function(map, row_vars) {
  if (!length(row_vars)) {
    return(rep("1", nrow(map)))
  }
  codes <- lapply(row_vars, function(rv) {
    x <- as.character(map[[rv]])
    match(x, unique(x))
  })
  do.call(paste, c(codes, list(sep = ".")))
}

.dkge_aggregate_row_ids <- function(row_meta, row_vars, row_sep) {
  if (!length(row_vars)) {
    return("aggregate")
  }
  ids <- as.character(do.call(paste, c(row_meta[row_vars], sep = row_sep)))
  # Belt and braces: `.dkge_aggregate_prepare_spec()` rejects labels containing
  # `row_sep`, but row IDs are also the handle used to match a user-supplied
  # kernel, so a collision must never escape this function silently.
  if (anyDuplicated(ids) || any(!nzchar(ids))) {
    dup <- unique(ids[duplicated(ids)])
    stop(sprintf(paste0("Aggregate row IDs must be unique and non-empty; ",
                        "'%s' is produced by more than one row. Choose a ",
                        "`row_sep` that cannot occur inside a label."),
                 paste(dup, collapse = "', '")), call. = FALSE)
  }
  ids
}

#' Fit a DKGE-regularized aggregate decomposition
#'
#' Fits a row-space decomposition of an aggregate target using only q x q
#' linear algebra, where q is the number of aggregate rows. This is intended
#' for cell-mean/PLS-bridge analyses, not subject-level DKGE inference.
#'
#' @param target A `dkge_aggregate_target` or row-by-feature numeric matrix.
#' @param K Optional design kernel, or an object returned by
#'   [design_kernel()], for the *aggregate row* space. Note that `q` here is the
#'   number of aggregate rows (for example `group:task:measure` cells), not the
#'   number of subject-level GLM effects used elsewhere in the package: `K` must
#'   therefore be `nrow(target$Y)` by `nrow(target$Y)`. When `K` carries
#'   dimnames they are matched against the aggregate row IDs; a kernel with only
#'   row names (or only column names) is validated and reordered on whichever
#'   labels are present. Rank-deficient PSD kernels keep a true square
#'   root: null directions stay at zero rather than receiving the jitter
#'   that [`.dkge_kernel_roots()`] uses for invertibility elsewhere.
#' @param rank Number of components to retain.
#' @param center Centering applied to the aggregate matrix before fitting.
#'   The default `"none"` keeps the grand mean in the decomposition, which makes
#'   the leading component largely a mean-level effect; under
#'   [dkge_aggregate_permute()] with `statistic = "singular_value"` that mean
#'   survives subject-label permutation, so the resulting test has little power.
#'   Use `"grand"`/`"column"` centering, or a contrast statistic, when the
#'   question is about differences between aggregate rows.
#'
#' @return Object of class `dkge_aggregate_fit`. Notable fields:
#'   \describe{
#'     \item{U}{q x rank K-orthonormal row basis.}
#'     \item{saliences}{\code{K \%*\% U}.}
#'     \item{scores_feature}{p x rank feature-space scores, \code{t(Yc) \%*\% saliences}.}
#'     \item{singular_values}{Per-component energies,
#'       \code{sqrt(colSums(scores_feature^2))}. For an unrotated fit these equal
#'       \code{sqrt(eig_values[seq_len(rank)])} exactly; after
#'       [dkge_aggregate_align()] they are the energies of the *rotated*
#'       components and are no longer sorted.}
#'     \item{eig_values}{Full length-q eigenvalue spectrum of \code{Chat}.}
#'     \item{Chat}{\eqn{K^{1/2} Y_c Y_c^\top K^{1/2}}. Both \code{Chat} and
#'       \code{eig_values} describe the data in the kernel metric and are
#'       invariant to the component rotation applied by [dkge_aggregate_align()];
#'       only \code{U}, \code{saliences}, \code{scores_feature}, and
#'       \code{singular_values} are basis-dependent.}
#'   }
#' @export
dkge_aggregate_fit <- function(target,
                               K = NULL,
                               rank = NULL,
                               center = c("none", "grand", "row", "column")) {
  center <- match.arg(center)
  target <- .dkge_as_aggregate_target(target)
  Y <- target$Y
  Yc <- .dkge_aggregate_center(Y, center)
  q <- nrow(Yc)

  kernel_info <- NULL
  if (is.null(K)) {
    K <- diag(q)
  } else if (is.list(K) && !is.null(K$K)) {
    kernel_info <- K$info %||% NULL
    K <- K$K
  }
  K <- .dkge_match_aggregate_kernel(K, rownames(Yc))
  K <- .dkge_validate_kernel(K)

  rank <- rank %||% min(q, ncol(Yc))
  rank <- as.integer(rank)
  if (length(rank) != 1L || is.na(rank) || rank < 1L) {
    stop("`rank` must be a positive integer.", call. = FALSE)
  }
  rank <- min(rank, q, ncol(Yc))

  roots <- .dkge_aggregate_kernel_roots(K)
  Chat <- roots$Khalf %*% (Yc %*% t(Yc)) %*% roots$Khalf
  Chat <- (Chat + t(Chat)) / 2
  eg <- eigen(Chat, symmetric = TRUE)
  vals <- pmax(eg$values, 0)
  chat_tol <- 1e-10 * max(1, max(vals))
  rank <- min(rank, max(1L, sum(vals > chat_tol)))
  keep <- seq_len(rank)
  U <- roots$Kihalf %*% eg$vectors[, keep, drop = FALSE]
  colnames(U) <- paste0("LV", keep)
  rownames(U) <- rownames(Yc)
  saliences <- K %*% U
  dimnames(saliences) <- dimnames(U)
  scores_feature <- t(Yc) %*% saliences
  rownames(scores_feature) <- colnames(Yc)
  colnames(scores_feature) <- colnames(U)
  # Energy is the feature-score norm so it stays identical after alignment
  # and cannot pick up jitter that K itself annihilates.
  singular_values <- sqrt(colSums(scores_feature^2))
  names(singular_values) <- colnames(U)

  out <- list(
    U = U,
    saliences = saliences,
    scores_feature = scores_feature,
    singular_values = singular_values,
    eig_values = vals,
    Chat = Chat,
    K = K,
    target = target,
    Y = Yc,
    center = center,
    rank = rank,
    kernel_info = kernel_info,
    estimand = "aggregate_cell_mean",
    call = match.call()
  )
  class(out) <- c("dkge_aggregate_fit", "list")
  out
}

#' True PSD square-root / pseudoinverse for aggregate kernels
#'
#' Unlike [`.dkge_kernel_roots()`], null eigenvalues stay zero. Jittering them
#' would put energy into `Chat` in directions that `K` still annihilates, so
#' `singular_values` and `scores_feature` would disagree.
#'
#' @keywords internal
#' @noRd
.dkge_aggregate_kernel_roots <- function(K, tol = NULL) {
  Ksym <- (K + t(K)) / 2
  ee <- eigen(Ksym, symmetric = TRUE)
  vals <- pmax(ee$values, 0)
  scale <- max(1, max(vals))
  tol <- tol %||% (1e-10 * scale)
  pos <- vals > tol
  sqrt_vals <- sqrt(vals)
  inv_sqrt <- numeric(length(vals))
  inv_sqrt[pos] <- 1 / sqrt_vals[pos]
  V <- ee$vectors
  n <- length(vals)
  list(
    Khalf = V %*% diag(sqrt_vals, n) %*% t(V),
    Kihalf = V %*% diag(inv_sqrt, n) %*% t(V),
    evals = vals,
    rank = sum(pos)
  )
}

.dkge_as_aggregate_target <- function(target) {
  if (inherits(target, "dkge_aggregate_target")) {
    return(target)
  }
  if (is.matrix(target) || is.data.frame(target)) {
    Y <- as.matrix(target)
    row_ids <- rownames(Y) %||% paste0("row", seq_len(nrow(Y)))
    feature_ids <- colnames(Y) %||% paste0("feature", seq_len(ncol(Y)))
    rownames(Y) <- row_ids
    colnames(Y) <- feature_ids
    out <- list(
      Y = Y,
      row_metadata = data.frame(row = row_ids, stringsAsFactors = FALSE),
      row_ids = row_ids,
      subject_ids = character(0),
      source_subject_ids = character(0),
      feature_ids = feature_ids,
      subject_to_row = data.frame(),
      weights = NULL,
      aggregate = "matrix",
      group_vars = character(0),
      cell_vars = "row",
      resample_spec = NULL,
      provenance = list(estimand = "aggregate_cell_mean",
                        inference_unit = "unknown")
    )
    class(out) <- c("dkge_aggregate_target", "list")
    return(out)
  }
  stop("`target` must be a dkge_aggregate_target or numeric matrix.",
       call. = FALSE)
}

.dkge_aggregate_center <- function(Y, center) {
  if (identical(center, "none")) {
    return(Y)
  }
  if (identical(center, "grand")) {
    return(Y - mean(Y))
  }
  if (identical(center, "row")) {
    return(Y - rowMeans(Y))
  }
  sweep(Y, 2L, colMeans(Y), "-")
}

.dkge_match_aggregate_kernel <- function(K, row_ids) {
  K <- as.matrix(K)
  rn <- rownames(K)
  cn <- colnames(K)
  if (!is.null(rn) || !is.null(cn)) {
    check_labels <- function(lab, what) {
      if (anyDuplicated(lab) || any(!nzchar(lab))) {
        stop(sprintf("Kernel %s must be unique and non-empty.", what),
             call. = FALSE)
      }
      if (length(lab) != length(row_ids) || !setequal(lab, row_ids)) {
        stop("Kernel dimnames must match aggregate row IDs.", call. = FALSE)
      }
    }
    ridx <- NULL
    cidx <- NULL
    if (!is.null(rn)) {
      check_labels(rn, "row names")
      ridx <- match(row_ids, rn)
    }
    if (!is.null(cn)) {
      check_labels(cn, "column names")
      cidx <- match(row_ids, cn)
    }
    ridx <- ridx %||% cidx
    cidx <- cidx %||% ridx
    K <- K[ridx, cidx, drop = FALSE]
    dimnames(K) <- list(row_ids, row_ids)
  }
  if (nrow(K) != length(row_ids) || ncol(K) != length(row_ids)) {
    stop("Kernel dimensions must match aggregate rows.", call. = FALSE)
  }
  K
}

#' Evaluate aggregate-fit statistics
#'
#' Aggregate statistics operate on a cell-mean decomposition, not on a
#' subject-level second-level model. In particular, `"between_group_contrast"`
#' is a convenience bridge statistic: it applies a supplied group-by-cell row
#' contrast to the selected component's singular-value-scaled aggregate
#' saliences. It should not be interpreted as a replacement for
#' [dkge_between_rrr()] or [dkge_between_permute()] on subject-level targets.
#'
#' @param fit A `dkge_aggregate_fit`.
#' @param statistic Built-in statistic name or a function accepting `fit`.
#'   `"singular_value"` returns the selected component singular value.
#'   `"salience_contrast"` projects the selected component salience vector onto
#'   `contrast`. `"component_score_contrast"` applies the same contrast to the
#'   component's singular-value-scaled salience vector. `"contrast_score"` is a
#'   legacy alias for `"salience_contrast"`. `"between_group_contrast"` is a
#'   bridge-analysis alias for `"component_score_contrast"` and should be
#'   supplied a group-by-within row contrast.
#' @param component Component index.
#' @param contrast Numeric row contrast. Required for every statistic except
#'   `"singular_value"`, i.e. for `"salience_contrast"`,
#'   `"component_score_contrast"`, `"contrast_score"`, and
#'   `"between_group_contrast"`. May be named, in which case it is reordered to
#'   the aggregate row IDs.
#' @param ... Additional arguments passed to a user-supplied statistic
#'   function. Built-in statistics take only `component` and `contrast`, so any
#'   other argument reaching `...` is reported as an error rather than being
#'   silently ignored.
#'
#' @return Numeric scalar statistic. Contrast statistics are returned
#'   **signed**; the sign is only interpretable relative to the component
#'   orientation of the fit the statistic is evaluated on.
#'   [dkge_aggregate_bootstrap()] orients resampled fits with
#'   [dkge_aggregate_align()] before taking the statistic;
#'   [dkge_aggregate_permute()] does not (null draws stay unaligned; two-sided
#'   `abs()` absorbs the arbitrary component sign).
#'   [dkge_aggregate_bootstrap()] reports a percentile interval on the signed
#'   statistic so that it can legitimately contain zero.
#' @export
dkge_aggregate_stat <- function(fit,
                                statistic = c("singular_value",
                                              "salience_contrast",
                                              "component_score_contrast",
                                              "contrast_score",
                                              "between_group_contrast"),
                                component = 1L,
                                contrast = NULL,
                                ...) {
  stopifnot(inherits(fit, "dkge_aggregate_fit"))
  if (is.function(statistic)) {
    out <- statistic(fit, ...)
    if (!is.numeric(out) || length(out) != 1L || !is.finite(out)) {
      stop("User-supplied aggregate statistic must return one finite number.",
           call. = FALSE)
    }
    return(unname(out))
  }
  statistic <- match.arg(statistic)
  # Built-in statistics take exactly `component` and `contrast`. Anything else
  # reaching `...` is a typo (`sed = 1` for `seed = 1`, `contrasts = ` for
  # `contrast = `) that would otherwise be swallowed in silence, including when
  # it was passed to `dkge_aggregate_permute()`/`dkge_aggregate_bootstrap()`.
  dots <- list(...)
  if (length(dots)) {
    nms <- names(dots)
    nms <- if (is.null(nms)) rep("", length(dots)) else nms
    labels <- ifelse(nzchar(nms), nms, "<unnamed>")
    stop(sprintf(paste0("Unused argument(s) for built-in aggregate statistic ",
                        "'%s': %s. Only `component` and `contrast` apply; ",
                        "`...` is forwarded only to user-supplied statistic ",
                        "functions."),
                 statistic, paste(labels, collapse = ", ")), call. = FALSE)
  }
  component <- as.integer(component)
  if (length(component) != 1L || is.na(component) ||
      component < 1L || component > fit$rank) {
    stop("`component` must select a fitted aggregate component.",
         call. = FALSE)
  }
  if (identical(statistic, "singular_value")) {
    return(unname(fit$singular_values[[component]]))
  }

  contrast <- .dkge_aggregate_match_contrast(contrast, rownames(fit$saliences))
  salience_score <- sum(contrast * fit$saliences[, component])
  if (identical(statistic, "component_score_contrast") ||
      identical(statistic, "between_group_contrast")) {
    salience_score <- salience_score * fit$singular_values[[component]]
  }
  unname(salience_score)
}

.dkge_aggregate_match_contrast <- function(contrast, row_ids) {
  if (is.null(contrast)) {
    stop("`contrast` is required for this aggregate statistic.", call. = FALSE)
  }
  names_in <- names(contrast)
  contrast <- as.numeric(contrast)
  if (length(contrast) != length(row_ids)) {
    stop("`contrast` must have one entry per aggregate row.", call. = FALSE)
  }
  if (!is.null(names_in) && all(nzchar(names_in))) {
    ord <- match(row_ids, names_in)
    if (anyNA(ord)) {
      stop("Named `contrast` must contain all aggregate row IDs.",
           call. = FALSE)
    }
    contrast <- contrast[ord]
  }
  contrast
}

#' Align aggregate components to a reference fit
#'
#' Aligns a resampled aggregate fit to a reference using sign alignment for
#' rank one and K-Procrustes alignment for rank greater than one.
#'
#' @param reference Reference `dkge_aggregate_fit`.
#' @param fit Fit to align.
#' @param rank Optional rank truncation before alignment.
#' @param degeneracy_tol Relative singular-value gap below which components
#'   are flagged as near-tied. The gap is computed on the *unrotated* (sorted)
#'   spectrum of `fit`, which is the quantity that governs how identifiable the
#'   rotation is.
#'
#' @return Aligned `dkge_aggregate_fit` with `alignment` metadata. `U`,
#'   `saliences`, `scores_feature`, and `singular_values` are all carried
#'   through the rotation together: `singular_values` is recomputed as
#'   `sqrt(colSums(scores_feature^2))`, which reproduces the unaligned singular
#'   values exactly when the rotation is the identity. `Chat` and `eig_values`
#'   are rotation-invariant properties of the data and are left untouched, so
#'   after alignment `singular_values` no longer equals
#'   `sqrt(eig_values[seq_len(rank)])`.
#' @export
dkge_aggregate_align <- function(reference,
                                 fit,
                                 rank = NULL,
                                 degeneracy_tol = 1e-6) {
  stopifnot(inherits(reference, "dkge_aggregate_fit"),
            inherits(fit, "dkge_aggregate_fit"))
  if (!identical(dim(reference$K), dim(fit$K)) ||
      max(abs(reference$K - fit$K)) > 1e-8 * max(1, max(abs(reference$K)))) {
    stop("Reference and fit kernels must match for aggregate alignment.",
         call. = FALSE)
  }
  if (!identical(rownames(reference$U), rownames(fit$U))) {
    stop("Reference and fit must share the same aggregate row IDs (in the same order) for alignment.",
         call. = FALSE)
  }
  rank <- rank %||% min(reference$rank, fit$rank)
  rank <- as.integer(rank)
  if (length(rank) != 1L || is.na(rank) || rank < 1L) {
    stop("`rank` must be a positive integer.", call. = FALSE)
  }
  rank <- min(rank, reference$rank, fit$rank)

  Uref <- reference$U[, seq_len(rank), drop = FALSE]
  U <- fit$U[, seq_len(rank), drop = FALSE]
  K <- reference$K
  if (rank == 1L) {
    score <- as.numeric(t(Uref) %*% K %*% U)
    sgn <- if (score < 0) -1 else 1
    R <- matrix(sgn, 1, 1)
    U_aligned <- U * sgn
    cosines <- abs(score)
  } else {
    pr <- dkge_procrustes_K(Uref, U, K)
    R <- pr$R
    U_aligned <- pr$U_aligned
    cosines <- pr$cosines
  }
  # `U %*% R` drops the component names, and after alignment the columns
  # correspond to the *reference* components.
  dimnames(U_aligned) <- list(rownames(U), colnames(Uref))

  # Degeneracy is a property of the source spectrum, before rotation.
  sv_source <- fit$singular_values[seq_len(rank)]

  out <- fit
  out$U <- U_aligned
  out$rank <- rank
  out$saliences <- K %*% U_aligned
  dimnames(out$saliences) <- dimnames(U_aligned)
  out$scores_feature <- t(out$Y) %*% out$saliences
  rownames(out$scores_feature) <- colnames(out$Y)
  colnames(out$scores_feature) <- colnames(U_aligned)
  # Rotate the component energies with the components. For the identity
  # rotation this reproduces `sqrt(eig_values[seq_len(rank)])` exactly.
  out$singular_values <- sqrt(colSums(out$scores_feature^2))
  names(out$singular_values) <- colnames(U_aligned)
  gaps <- abs(diff(sv_source))
  near_tie <- length(gaps) > 0L &&
    any(gaps <= degeneracy_tol * max(1, sv_source[[1L]]))
  out$alignment <- list(
    method = if (rank == 1L) "sign" else "K-procrustes",
    R = R,
    cosines = cosines,
    near_tie = near_tie,
    singular_gaps = gaps,
    rank = rank
  )
  out
}

#' Permutation inference for aggregate targets
#'
#' Performs subject-label permutation for aggregate cell-mean targets. Each
#' draw permutes subject metadata, recomputes aggregate rows, refits the
#' aggregate decomposition, and evaluates the selected statistic on that
#' unaligned refit.
#'
#' Because a permuted subject label changes which aggregate row each subject
#' contributes to, `group_vars` must cover every subject-level variable that
#' defines the aggregate rows (`target$group_vars`). Permuting a strict subset
#' would let the joint set of row keys change from draw to draw, so it is
#' rejected up front rather than failing partway through the run.
#'
#' @param target A `dkge_aggregate_target`.
#' @param K Optional aggregate-row design kernel. See [dkge_aggregate_fit()] for
#'   the row-space convention (`q` = number of aggregate rows).
#' @param statistic Built-in statistic name or function.
#' @param B Number of permutation draws.
#' @param group_vars Subject-level labels to permute. Defaults to
#'   `target$group_vars`, and must include all of `target$group_vars`.
#' @param rank,center Passed to [dkge_aggregate_fit()].
#' @param component Component index passed to [dkge_aggregate_stat()] for the
#'   built-in statistics. It is a formal argument rather than part of `...`
#'   because R's partial matching would otherwise bind a `component =` argument
#'   to `component_scale`/`component_contrasts`. User-supplied statistic
#'   functions do not receive it.
#' @param alternative Direction of the permutation p-value. `"two.sided"`
#'   (the default) compares `abs(null)` to `abs(observed)`, which is the
#'   appropriate choice for the signed contrast statistics returned by
#'   [dkge_aggregate_stat()] because the component sign is arbitrary.
#'   `"greater"`/`"less"` compare the signed values, and `"greater"` is the
#'   natural choice for a non-negative statistic such as `"singular_value"`.
#' @param parallel Logical; if `TRUE`, evaluate draws with
#'   `future.apply::future_lapply()`. Permutation indices are drawn up front in
#'   the calling RNG stream, so results are identical for either setting of
#'   `parallel` given the same `seed`.
#' @param seed Optional random seed.
#' @param ... Additional arguments passed to [dkge_aggregate_stat()]. Every
#'   other setting is a named formal, so a misspelled argument (`sed = 1` for
#'   `seed = 1`) lands here and is reported as an error whenever `statistic` is
#'   one of the built-ins; `...` reaches only user-supplied statistic functions.
#'
#' @details
#' Null draws are evaluated **unaligned**. Aligning every permutation to the
#' observed fit shrinks component statistics (aligned
#' \eqn{sv_1 = \sqrt{\sum_k d_k^2 R_{k1}^2} \le d_1}) and invalidates the
#' null at rank greater than one. [dkge_aggregate_align()] is still run so
#' that `alignment_summary` can flag weak or near-tied components, but that
#' rotation is diagnostic only and is not passed to
#' [dkge_aggregate_stat()]. `"two.sided"` is the right default for
#' sign-ambiguous component statistics because the unaligned component sign
#' is arbitrary; `"greater"` remains the natural choice for the non-negative
#' `"singular_value"` statistic. Bootstrap resampling keeps alignment
#' because the observed fit is a legitimate reference for a CI.
#'
#' @return Object of class `dkge_aggregate_permutation`. `observed` and `null`
#'   hold signed statistics; `p` is computed according to `alternative`.
#' @export
dkge_aggregate_permute <- function(target,
                                   K = NULL,
                                   statistic = c("singular_value",
                                                 "salience_contrast",
                                                 "component_score_contrast",
                                                 "contrast_score",
                                                 "between_group_contrast"),
                                   B = 999L,
                                   group_vars = NULL,
                                   rank = NULL,
                                   center = c("none", "grand", "row", "column"),
                                   component = 1L,
                                   alternative = c("two.sided", "greater", "less"),
                                   parallel = FALSE,
                                   seed = NULL,
                                   ...) {
  stopifnot(inherits(target, "dkge_aggregate_target"))
  if (is.null(target$resample_spec)) {
    stop("Aggregate target does not carry resampling provenance.", call. = FALSE)
  }
  center <- match.arg(center)
  alternative <- match.arg(alternative)
  B <- .dkge_validate_resample_B(B)
  spec <- target$resample_spec
  group_vars <- as.character(group_vars %||% target$group_vars)
  if (!length(group_vars)) {
    stop("`group_vars` must identify subject labels to permute.", call. = FALSE)
  }
  unknown <- setdiff(group_vars, names(spec$subject_data))
  if (length(unknown)) {
    stop("Unknown subject grouping variable(s) in `group_vars`: ",
         paste(unknown, collapse = ", "), call. = FALSE)
  }
  fixed <- setdiff(target$group_vars, group_vars)
  if (length(fixed)) {
    stop(sprintf(paste0("`group_vars` must permute every subject variable that ",
                        "defines the aggregate rows; %s would stay fixed, which ",
                        "can change the set of aggregate rows from draw to draw. ",
                        "Permute c(%s), or rebuild the target with only the ",
                        "variable(s) you want to permute."),
                 paste(sprintf("'%s'", fixed), collapse = ", "),
                 paste(sprintf('"%s"', target$group_vars), collapse = ", ")),
         call. = FALSE)
  }

  stat_args <- .dkge_aggregate_stat_args(statistic, component, list(...))
  seed_state <- .dkge_seed_enter(seed)
  on.exit(.dkge_seed_exit(seed_state), add = TRUE)
  observed_fit <- dkge_aggregate_fit(target, K = K, rank = rank, center = center)
  observed <- do.call(dkge_aggregate_stat,
                      c(list(fit = observed_fit, statistic = statistic), stat_args))

  n_subj <- nrow(spec$subject_data)
  perm_idx <- lapply(seq_len(B), function(b) sample.int(n_subj))
  draw <- .dkge_aggregate_permute_worker(
    spec = spec,
    reference = .dkge_aggregate_alignment_reference(observed_fit),
    row_ids = target$row_ids,
    group_vars = group_vars,
    center = center,
    statistic = statistic,
    stat_args = stat_args
  )
  draws <- .dkge_aggregate_apply(perm_idx, draw, parallel = parallel)
  null <- vapply(draws, `[[`, numeric(1), "stat")
  alignments <- lapply(draws, `[[`, "alignment")

  alignment_summary <- .dkge_aggregate_alignment_summary(alignments)
  p <- .dkge_aggregate_permutation_p(observed, null, alternative)
  out <- list(
    observed = observed,
    null = null,
    p = p,
    B = B,
    alternative = alternative,
    statistic = if (is.function(statistic)) "function" else match.arg(statistic),
    fit = observed_fit,
    alignments = alignments,
    alignment_summary = alignment_summary,
    resampling = list(method = "subject_label_permutation",
                      group_vars = group_vars,
                      inference_unit = "subject",
                      seed = seed),
    call = match.call()
  )
  class(out) <- c("dkge_aggregate_permutation", "list")
  out
}

# A resampled fit is refit against the observed kernel/rank. Alignment
# only reads `K`, `U`, and `rank`. Carrying the whole observed fit into the draw
# closure would drag `fit$target$resample_spec` along -- a second reference to
# the stacked value matrix -- plus `Chat`, `Y`, and the feature scores, all of
# which `future.apply` re-serialises for every parallel chunk.
.dkge_aggregate_alignment_reference <- function(fit) {
  ref <- list(
    U = fit$U,
    K = fit$K,
    rank = fit$rank,
    singular_values = fit$singular_values,
    center = fit$center,
    estimand = fit$estimand
  )
  class(ref) <- class(fit)
  ref
}

# Every resample must reproduce the observed aggregate row set: the kernel and
# the reference basis are both indexed by it. Checking here turns an opaque
# "Kernel dimnames must match aggregate row IDs" failure -- possibly hundreds of
# draws in -- into a statement of what actually went wrong.
.dkge_aggregate_check_row_ids <- function(target_b, row_ids, context) {
  if (identical(target_b$row_ids, row_ids)) {
    return(invisible(NULL))
  }
  missing <- setdiff(row_ids, target_b$row_ids)
  detail <- if (length(missing)) {
    sprintf(" Row(s) absent from this draw: %s.",
            paste(sprintf("'%s'", missing), collapse = ", "))
  } else {
    ""
  }
  msg <- switch(
    context,
    permute = paste0("Permuting `group_vars` changed the aggregate row set; ",
                     "permute all row-defining subject variables jointly."),
    bootstrap = paste0("A bootstrap resample dropped an aggregate row. ",
                       "`strata` must nest within the target's `group_vars` so ",
                       "that every draw keeps every group level."),
    "A resample changed the aggregate row set."
  )
  stop(paste0(msg, detail), call. = FALSE)
}

# Workers are built by a factory so that their enclosing environment holds only
# what a draw needs. A closure defined inline in the resampling entry point
# instead captures that whole frame -- the observed fit, the pre-drawn index
# list, the bootstrap accumulators -- and `future.apply` exports all of it to
# every parallel chunk.
.dkge_aggregate_permute_worker <- function(spec,
                                           reference,
                                           row_ids,
                                           group_vars,
                                           center,
                                           statistic,
                                           stat_args) {
  force(spec); force(reference); force(row_ids); force(group_vars)
  force(center); force(statistic); force(stat_args)
  function(idx) {
    subject_data_b <- spec$subject_data
    subject_data_b[group_vars] <- subject_data_b[idx, group_vars, drop = FALSE]
    target_b <- .dkge_aggregate_from_spec(spec, subject_data = subject_data_b)
    .dkge_aggregate_check_row_ids(target_b, row_ids, "permute")
    fit_b <- dkge_aggregate_fit(target_b, K = reference$K,
                                rank = reference$rank, center = center)
    # Alignment is diagnostic only: the statistic is taken on the unaligned
    # refit so the null is not attracted to the observed subspace.
    alignment <- dkge_aggregate_align(reference, fit_b)$alignment
    list(stat = do.call(dkge_aggregate_stat,
                        c(list(fit = fit_b, statistic = statistic), stat_args)),
         alignment = alignment)
  }
}

.dkge_aggregate_bootstrap_worker <- function(spec,
                                             reference,
                                             row_ids,
                                             center,
                                             statistic,
                                             stat_args,
                                             component_contrasts,
                                             component_scale,
                                             want_contrast,
                                             want_features) {
  force(spec); force(reference); force(row_ids); force(center)
  force(statistic); force(stat_args); force(component_contrasts)
  force(component_scale); force(want_contrast); force(want_features)
  function(idx) {
    target_b <- .dkge_aggregate_from_spec(spec, subject_indices = idx)
    .dkge_aggregate_check_row_ids(target_b, row_ids, "bootstrap")
    fit_b <- dkge_aggregate_fit(target_b, K = reference$K,
                                rank = reference$rank, center = center)
    fit_b <- dkge_aggregate_align(reference, fit_b)
    list(
      stat = do.call(dkge_aggregate_stat,
                     c(list(fit = fit_b, statistic = statistic), stat_args)),
      alignment = fit_b$alignment,
      contrast = if (want_contrast) {
        .dkge_aggregate_component_contrasts(fit_b, component_contrasts,
                                            scale = component_scale)
      } else {
        NULL
      },
      scores_feature = if (want_features) fit_b$scores_feature else NULL
    )
  }
}

# `component` is a formal of the resampling entry points (rather than being
# forwarded through `...`) because partial matching would otherwise bind it to
# `component_scale`/`component_contrasts`. User-supplied statistic functions
# never receive it.
.dkge_aggregate_stat_args <- function(statistic, component, dots) {
  if (is.function(statistic)) {
    return(dots)
  }
  dots$component <- component
  dots
}

.dkge_aggregate_permutation_p <- function(observed, null, alternative) {
  B <- length(null)
  tol <- sqrt(.Machine$double.eps) * max(1, abs(observed))
  hits <- switch(
    alternative,
    two.sided = sum(abs(null) >= abs(observed) - tol),
    greater = sum(null >= observed - tol),
    less = sum(null <= observed + tol)
  )
  (1 + hits) / (B + 1)
}

# Serial by default; `future.apply` only when explicitly requested. All RNG
# draws happen in the caller, so the backend never affects reproducibility.
.dkge_aggregate_apply <- function(X, FUN, parallel = FALSE) {
  if (isTRUE(parallel)) {
    .dkge_apply(X, FUN, parallel = TRUE, future.seed = TRUE)
  } else {
    lapply(X, FUN)
  }
}

#' Bootstrap aggregate targets over subjects
#'
#' Resamples subjects with replacement, recomputes aggregate rows, refits the
#' aggregate decomposition, aligns components, and evaluates the statistic.
#'
#' Resampling is **stratified by default**: when `strata` is `NULL` and the
#' target has subject-level `group_vars`, subjects are resampled with
#' replacement *within* `interaction(subject_data[target$group_vars])`, so every
#' draw preserves the observed group sizes and the aggregate row set. Pass
#' `strata` explicitly (a vector, or the name of a column in the subject data)
#' to stratify differently; a stratum with a single subject contributes that
#' subject to every draw, which is reported with a warning.
#'
#' @inheritParams dkge_aggregate_permute
#' @param B Number of bootstrap resamples.
#' @param strata Optional subject-level strata for within-stratum bootstrap:
#'   either a vector with one entry per source subject, or the name of a column
#'   in the subject data. Defaults to the interaction of
#'   `target$group_vars`. Explicit strata **must** nest within `group_vars`
#'   (each stratum lies inside a single group level) and are rejected up front
#'   otherwise: a coarser or cross-cutting stratification can draw a resample in
#'   which an entire group level is absent, which changes the aggregate row set
#'   and fails against the fixed observed kernel. Every draw is additionally
#'   checked against the observed row set.
#' @param conf Confidence level for percentile intervals. Intervals use
#'   `stats::quantile(..., type = 6)` (Weibull plotting positions), which is
#'   slightly wider than the default type 7 in small samples.
#' @param component_contrasts Optional aggregate-row contrast matrix used to
#'   summarize each aligned component in every bootstrap draw.
#' @param component_scale Whether component contrasts are applied to
#'   singular-value-scaled cell scores (`"score"`) or raw saliences
#'   (`"salience"`).
#' @param return_features Logical; if `TRUE`, accumulate aligned
#'   feature-space component maps across bootstrap draws and return streaming
#'   mean, SD, and bootstrap-ratio z maps (`observed / bootstrap SD`, the
#'   standard PLS bootstrap ratio).
#'
#' @return Object of class `dkge_aggregate_bootstrap`. `observed` and
#'   `statistics` are signed, `interval` is the percentile interval of the
#'   signed statistic, and `excludes_zero` reports whether that interval
#'   excludes zero.
#' @export
dkge_aggregate_bootstrap <- function(target,
                                     K = NULL,
                                     statistic = c("singular_value",
                                                   "salience_contrast",
                                                   "component_score_contrast",
                                                   "contrast_score",
                                                   "between_group_contrast"),
                                     B = 999L,
                                     strata = NULL,
                                     rank = NULL,
                                     center = c("none", "grand", "row", "column"),
                                     component = 1L,
                                     conf = 0.95,
                                     component_contrasts = NULL,
                                     component_scale = c("score", "salience"),
                                     return_features = FALSE,
                                     parallel = FALSE,
                                     seed = NULL,
                                     ...) {
  stopifnot(inherits(target, "dkge_aggregate_target"))
  if (is.null(target$resample_spec)) {
    stop("Aggregate target does not carry resampling provenance.", call. = FALSE)
  }
  center <- match.arg(center)
  component_scale <- match.arg(component_scale)
  stat_args <- .dkge_aggregate_stat_args(statistic, component, list(...))
  B <- .dkge_validate_resample_B(B)
  if (!is.numeric(conf) || length(conf) != 1L || conf <= 0 || conf >= 1) {
    stop("`conf` must be a scalar in (0, 1).", call. = FALSE)
  }
  seed_state <- .dkge_seed_enter(seed)
  on.exit(.dkge_seed_exit(seed_state), add = TRUE)
  observed_fit <- dkge_aggregate_fit(target, K = K, rank = rank, center = center)
  observed <- do.call(dkge_aggregate_stat,
                      c(list(fit = observed_fit, statistic = statistic), stat_args))
  contrast_obs <- .dkge_aggregate_component_contrasts(
    observed_fit,
    component_contrasts,
    scale = component_scale
  )
  contrast_samples <- if (is.null(contrast_obs)) {
    NULL
  } else {
    array(NA_real_,
          dim = c(nrow(contrast_obs), ncol(contrast_obs), B),
          dimnames = c(dimnames(contrast_obs), list(sample = seq_len(B))))
  }
  feature_boot <- if (isTRUE(return_features)) {
    p <- nrow(observed_fit$scores_feature)
    r <- observed_fit$rank
    list(
      n = 0L,
      mean = matrix(0, p, r, dimnames = dimnames(observed_fit$scores_feature)),
      m2 = matrix(0, p, r, dimnames = dimnames(observed_fit$scores_feature))
    )
  } else {
    NULL
  }
  spec <- target$resample_spec
  if (is.null(strata)) {
    if (length(target$group_vars)) {
      strata <- interaction(spec$subject_data[target$group_vars],
                            drop = TRUE,
                            sep = target$resample_spec$row_sep)
    } else {
      strata <- rep("all", length(spec$subject_ids))
    }
  } else if (is.character(strata) && length(strata) == 1L &&
             strata %in% names(spec$subject_data)) {
    strata <- spec$subject_data[[strata]]
  }
  if (length(strata) != length(spec$subject_ids)) {
    stop("`strata` must match the number of source subjects.", call. = FALSE)
  }
  if (anyNA(strata)) {
    stop("`strata` must not contain missing values.", call. = FALSE)
  }
  strata <- droplevels(as.factor(strata))
  sizes <- table(strata)
  if (any(sizes < 2L)) {
    singleton <- names(sizes)[sizes < 2L]
    warning(sprintf(paste0("Bootstrap strata with fewer than 2 subjects (%s) ",
                           "contribute the same subject to every draw; their ",
                           "sampling variability is not represented."),
                    paste(singleton, collapse = ", ")), call. = FALSE)
  }
  # A stratification coarser than `group_vars` can draw a resample with a whole
  # group level missing, which changes the aggregate row set and then fails
  # against the fixed observed kernel -- possibly hundreds of draws in, with a
  # message about kernel dimnames. Reject it up front instead.
  .dkge_aggregate_validate_strata_nesting(strata, spec$subject_data,
                                          target$group_vars)

  strata_index <- split(seq_along(spec$subject_ids), strata)
  boot_idx <- lapply(seq_len(B), function(b) {
    .dkge_aggregate_bootstrap_indices(strata_index)
  })
  draw <- .dkge_aggregate_bootstrap_worker(
    spec = spec,
    reference = .dkge_aggregate_alignment_reference(observed_fit),
    row_ids = target$row_ids,
    center = center,
    statistic = statistic,
    stat_args = stat_args,
    component_contrasts = component_contrasts,
    component_scale = component_scale,
    want_contrast = !is.null(contrast_samples),
    want_features = !is.null(feature_boot)
  )

  stats <- numeric(B)
  alignments <- vector("list", B)
  # Chunked so that the feature accumulator never holds more than `chunk`
  # p x rank matrices at once, even when draws are evaluated in parallel.
  chunk <- if (is.null(feature_boot)) B else min(B, 64L)
  for (start in seq(1L, B, by = chunk)) {
    ids <- seq.int(start, min(start + chunk - 1L, B))
    draws <- .dkge_aggregate_apply(boot_idx[ids], draw, parallel = parallel)
    for (k in seq_along(ids)) {
      b <- ids[[k]]
      stats[[b]] <- draws[[k]]$stat
      alignments[[b]] <- draws[[k]]$alignment
      if (!is.null(contrast_samples)) {
        contrast_samples[, , b] <- draws[[k]]$contrast
      }
      if (!is.null(feature_boot)) {
        feature_boot <- .dkge_welford_update_matrix(feature_boot,
                                                    draws[[k]]$scores_feature)
      }
    }
    rm(draws)
  }
  alignment_summary <- .dkge_aggregate_alignment_summary(alignments)
  alpha <- (1 - conf) / 2
  interval <- stats::quantile(stats, probs = c(alpha, 1 - alpha),
                              names = FALSE, type = 6)
  component_contrast_out <- .dkge_aggregate_component_contrast_summary(
    observed = contrast_obs,
    samples = contrast_samples,
    conf = conf
  )
  feature_boot_out <- .dkge_feature_bootstrap_finalize(
    feature_boot,
    observed = observed_fit$scores_feature
  )
  out <- list(
    observed = observed,
    statistics = stats,
    interval = interval,
    excludes_zero = interval[[1]] > 0 || interval[[2]] < 0,
    conf = conf,
    B = B,
    statistic = if (is.function(statistic)) "function" else match.arg(statistic),
    fit = observed_fit,
    alignments = alignments,
    alignment_summary = alignment_summary,
    component_contrasts = component_contrast_out,
    feature_bootstrap = feature_boot_out,
    resampling = list(method = "subject_bootstrap",
                      strata = levels(strata),
                      inference_unit = "subject",
                      seed = seed),
    call = match.call()
  )
  class(out) <- c("dkge_aggregate_bootstrap", "list")
  out
}

# Strata must nest within the subject-level variables that define aggregate
# rows: every stratum has to sit inside a single `group_vars` combination, so
# that resampling within strata preserves each group's size and therefore the
# row set. Finer strata are fine; coarser or cross-cutting ones are not.
.dkge_aggregate_validate_strata_nesting <- function(strata, subject_data,
                                                    group_vars) {
  group_vars <- as.character(group_vars %||% character(0))
  if (!length(group_vars)) {
    return(invisible(NULL))
  }
  group_key <- do.call(paste, c(lapply(group_vars, function(v) {
    as.character(subject_data[[v]])
  }), list(sep = "\r")))
  by_stratum <- split(group_key, strata)
  mixed <- vapply(by_stratum, function(g) length(unique(g)) > 1L, logical(1))
  if (!any(mixed)) {
    return(invisible(NULL))
  }
  bad <- names(by_stratum)[mixed]
  stop(sprintf(paste0("`strata` must nest within the aggregate target's ",
                      "group variable(s) (%s); stratum %s spans more than one ",
                      "group level. A coarser stratification can draw a ",
                      "resample in which a whole group level is absent, which ",
                      "changes the aggregate row set."),
               paste(sprintf("'%s'", group_vars), collapse = ", "),
               paste(sprintf("'%s'", bad), collapse = ", ")),
       call. = FALSE)
}

# One stratified bootstrap resample. `strata_index` is a list of integer
# subject positions per stratum. Indexing `ii` explicitly is required: with a
# length-one stratum, `sample(ii, 1)` would sample from `seq_len(ii)`.
.dkge_aggregate_bootstrap_indices <- function(strata_index) {
  idx <- lapply(strata_index, function(ii) {
    n <- length(ii)
    if (!n) integer(0) else ii[sample.int(n, n, replace = TRUE)]
  })
  as.integer(unlist(idx, use.names = FALSE))
}

.dkge_aggregate_alignment_summary <- function(alignments, weak_cosine = 0.5) {
  if (!length(alignments)) {
    return(list(
      min_cosine = NA_real_,
      median_cosine = NA_real_,
      n_weak = 0L,
      n_near_tie = 0L,
      weak_cosine = weak_cosine,
      weak_alignment = FALSE
    ))
  }
  cosines <- unlist(lapply(alignments, `[[`, "cosines"), use.names = FALSE)
  finite_cos <- cosines[is.finite(cosines)]
  near_tie <- vapply(alignments, function(x) isTRUE(x$near_tie), logical(1))
  n_weak <- sum(finite_cos < weak_cosine)
  list(
    min_cosine = if (length(finite_cos)) min(finite_cos) else NA_real_,
    median_cosine = if (length(finite_cos)) stats::median(finite_cos) else NA_real_,
    n_weak = as.integer(n_weak),
    n_near_tie = as.integer(sum(near_tie)),
    weak_cosine = weak_cosine,
    weak_alignment = isTRUE(n_weak > 0L) || any(near_tie)
  )
}

.dkge_aggregate_component_contrasts <- function(fit,
                                                contrasts = NULL,
                                                scale = c("score", "salience")) {
  if (is.null(contrasts)) {
    return(NULL)
  }
  scale <- match.arg(scale)
  C <- as.matrix(contrasts)
  if (!is.numeric(C) || any(!is.finite(C))) {
    stop("`component_contrasts` must be a finite numeric matrix.", call. = FALSE)
  }
  row_ids <- rownames(fit$saliences)
  if (!is.null(rownames(C))) {
    ord <- match(row_ids, rownames(C))
    if (anyNA(ord)) {
      stop("`component_contrasts` row names must contain all aggregate row IDs.",
           call. = FALSE)
    }
    C <- C[ord, , drop = FALSE]
  }
  if (nrow(C) != length(row_ids)) {
    stop("`component_contrasts` must have one row per aggregate row.",
         call. = FALSE)
  }
  if (is.null(colnames(C))) {
    colnames(C) <- paste0("contrast", seq_len(ncol(C)))
  }

  cell_scores <- fit$saliences
  if (identical(scale, "score")) {
    cell_scores <- sweep(cell_scores, 2L, fit$singular_values, "*")
  }
  out <- t(C) %*% cell_scores
  rownames(out) <- colnames(C)
  colnames(out) <- colnames(fit$saliences)
  out
}

.dkge_aggregate_component_contrast_summary <- function(observed,
                                                       samples,
                                                       conf = 0.95) {
  if (is.null(observed)) {
    return(NULL)
  }
  alpha <- (1 - conf) / 2
  rows <- vector("list", length(observed))
  k <- 0L
  for (i in seq_len(nrow(observed))) {
    for (j in seq_len(ncol(observed))) {
      k <- k + 1L
      vals <- samples[i, j, ]
      qs <- stats::quantile(vals, probs = c(alpha, 1 - alpha),
                            names = FALSE, type = 6)
      rows[[k]] <- data.frame(
        contrast = rownames(observed)[[i]],
        component = colnames(observed)[[j]],
        observed = observed[i, j],
        bootstrap_low = qs[[1]],
        bootstrap_high = qs[[2]],
        excludes_zero = qs[[1]] > 0 || qs[[2]] < 0,
        stringsAsFactors = FALSE
      )
    }
  }
  list(
    observed = observed,
    samples = samples,
    summary = do.call(rbind, rows),
    conf = conf
  )
}

.dkge_welford_update_matrix <- function(state, x) {
  x <- as.matrix(x)
  if (!identical(dim(x), dim(state$mean))) {
    stop("Feature bootstrap matrix dimensions changed across resamples.",
         call. = FALSE)
  }
  state$n <- state$n + 1L
  delta <- x - state$mean
  state$mean <- state$mean + delta / state$n
  state$m2 <- state$m2 + delta * (x - state$mean)
  state
}

.dkge_feature_bootstrap_finalize <- function(state, observed) {
  if (is.null(state)) {
    return(NULL)
  }
  sd <- if (state$n > 1L) {
    sqrt(state$m2 / (state$n - 1L))
  } else {
    matrix(NA_real_, nrow = nrow(state$mean), ncol = ncol(state$mean),
           dimnames = dimnames(state$mean))
  }
  # Standard PLS bootstrap ratio: observed loading over its bootstrap SD.
  # A (near-)zero SD carries no information about stability, so report NA
  # rather than an arbitrarily large ratio.
  scale_ref <- suppressWarnings(max(abs(observed)))
  if (!is.finite(scale_ref) || scale_ref <= 0) {
    scale_ref <- 1
  }
  z <- observed / sd
  z[!is.finite(sd) | sd <= 1e-12 * scale_ref] <- NA_real_
  list(
    observed = observed,
    mean = state$mean,
    sd = sd,
    z = z,
    n = state$n
  )
}

# `.dkge_validate_resample_B()` and `.dkge_seed_enter()`/`.dkge_seed_exit()`
# live in R/dkge-utils.R; the aggregate and between-subject resampling entry
# points share the single definition there.

#' @export
#' @rdname dkge_aggregate_target
print.dkge_aggregate_target <- function(x, ...) {
  cat("<dkge_aggregate_target>\n")
  cat("  rows:     ", nrow(x$Y), "\n", sep = "")
  cat("  features: ", ncol(x$Y), "\n", sep = "")
  cat("  aggregate:", x$aggregate, "\n")
  invisible(x)
}

#' @export
#' @rdname dkge_aggregate_fit
print.dkge_aggregate_fit <- function(x, ...) {
  cat("<dkge_aggregate_fit>\n")
  cat("  estimand: ", x$estimand, "\n", sep = "")
  cat("  rows:     ", nrow(x$Y), "\n", sep = "")
  cat("  rank:     ", x$rank, "\n", sep = "")
  cat("  singular: ", paste(signif(x$singular_values, 4), collapse = ", "), "\n", sep = "")
  invisible(x)
}

#' @export
#' @rdname dkge_aggregate_permute
print.dkge_aggregate_permutation <- function(x, ...) {
  cat("<dkge_aggregate_permutation>\n")
  cat("  statistic: ", signif(x$observed, 4), "\n", sep = "")
  cat("  p:         ", signif(x$p, 4), " (", x$alternative %||% "two.sided", ")\n",
      sep = "")
  cat("  B:         ", x$B, "\n", sep = "")
  invisible(x)
}

#' @export
#' @rdname dkge_aggregate_bootstrap
print.dkge_aggregate_bootstrap <- function(x, ...) {
  cat("<dkge_aggregate_bootstrap>\n")
  cat("  statistic: ", signif(x$observed, 4), "\n", sep = "")
  cat("  interval:  ", paste(signif(x$interval, 4), collapse = ", "),
      if (isTRUE(x$excludes_zero)) "  (excludes 0)" else "  (includes 0)",
      "\n", sep = "")
  cat("  B:         ", x$B, "\n", sep = "")
  invisible(x)
}
