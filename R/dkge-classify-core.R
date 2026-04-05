# dkge-classify-core.R
# Classification add-on for DKGE fits — public API, target prep, S3 methods.

#' Cross-validated classification on DKGE effect patterns
#'
#' @param fit dkge object.
#' @param targets Target specification consumed by [dkge_targets()] or a list of
#'   `dkge_target` objects.
#' @details
#' Anchor-based fits produced by [dkge_anchor_fit()] do not retain the
#' design-factor metadata that `dkge_targets()` expects. In that setting you must
#' supply explicit weight matrices (rows = classes, columns = effects) or
#' ready-made [`dkge_target`] objects—helpers such as
#' [dkge_anchor_targets_from_prototypes()] and
#' [dkge_anchor_targets_from_directions()] can be used to construct them.
#' @param y Optional subject-level labels for delta-mode targets. Can be a vector
#'   (recycled across all delta targets) or a list matching `targets`; values are
#'   coerced to factors using each target's `class_labels`.
#' @param method Classifier backend: "lda" (default) or "logit" (one-vs-rest
#'   ridge logistic regression).
#' @param folds Cross-fitting specification. `NULL` (default) performs LOSO.
#'   Integer values request subject-level K-fold. A `dkge_folds` object is also
#'   accepted. Other coercible inputs are routed through [as_dkge_folds()].
#' @param lambda Optional ridge parameter passed to the classifier backend.
#' @param metric Evaluation metric(s); defaults to c("accuracy", "logloss").
#' @param mode Decoding mode: "auto" (default), "cell", "cell_cross", or
#'   "delta". Cell-cross trains on held-in subjects and tests generalisation to
#'   the held-out subject.
#' @param residualize Forwarded to [dkge_targets()].
#' @param collapse Forwarded to [dkge_targets()] for factor collapsing.
#' @param restrict_factors Optional factor subset for `spec = "fullcell"`.
#' @param n_perm Number of permutations for empirical p-values (default 0).
#' @param scope Override permutation scope (otherwise target scope used).
#' @param class_weights Class weighting scheme ("none", "balanced", "inverse").
#' @param ridge Optional ridge added when recomputing held-out bases (default 0).
#' @param control Optional list of advanced controls (power users). Recognised
#'   entries: `lambda_grid` (numeric vector of candidate penalties) and
#'   `lambda_fun` (function returning a lambda per target/fold with signature
#'   `function(target, fold, method, default)`). Defaults to `NULL`, leaving the
#'   standard `lambda` behaviour unchanged.
#' @param blocks Optional vector identifying within-subject blocks (e.g., runs or
#'   sessions) used to constrain permutations. Length must match the number of
#'   subjects in `fit` when supplied.
#' @param standardize_within_fold Logical indicating whether to z-score features
#'   using training data inside each fold. When `NULL` (default), standardisation
#'   is enabled automatically for `mode = "cell_cross"` and disabled otherwise.
#' @param parallel Logical; reserved for future parallelism hooks.
#' @param verbose Logical; print progress messages.
#' @param seed Optional random seed applied before permutations.
#' @examples
#' # Simulate toy data with design kernel (includes kernel_info for targets)
#' toy <- dkge_sim_toy(
#'   factors = list(A = list(L = 2), B = list(L = 3)),
#'   active_terms = c("A", "B"), S = 6, P = 20, snr = 5
#' )
#' # Use design_kernel() result to retain factor metadata
#' kern <- design_kernel(
#'   factors = list(A = list(L = 2), B = list(L = 3)),
#'   basis = "effect"
#' )
#' fit <- dkge(toy$B_list, toy$X_list, kernel = kern, rank = 2)
#'
#' # Decode factor A using cell-mode classification
#' \donttest{
#' clf <- dkge_classify(fit, targets = ~A, method = "lda")
#' clf
#' }
#' @export
dkge_classify <- function(fit,
                          targets,
                          y = NULL,
                          method = c("lda", "logit"),
                          folds = NULL,
                          lambda = NULL,
                          metric = c("accuracy", "logloss"),
                          mode = c("auto", "cell", "cell_cross", "delta"),
                          standardize_within_fold = NULL,
                          residualize = TRUE,
                          collapse = NULL,
                          restrict_factors = NULL,
                          n_perm = 0L,
                          scope = NULL,
                          class_weights = c("none", "balanced", "inverse"),
                          ridge = 0,
                          control = NULL,
                          blocks = NULL,
                          parallel = FALSE,
                          verbose = FALSE,
                          seed = NULL) {
  stopifnot(inherits(fit, "dkge"))
  method <- match.arg(method)
  metric <- unique(metric)
  mode <- match.arg(mode)
  class_weights <- match.arg(class_weights)
  if (!is.null(seed)) set.seed(seed)

  control <- control %||% list()
  if (!is.list(control)) {
    stop("`control` must be a list or NULL")
  }
  lambda_grid <- control$lambda_grid %||% NULL
  if (!is.null(lambda_grid)) {
    if (!is.numeric(lambda_grid) || any(lambda_grid <= 0)) {
      stop("`control$lambda_grid` must contain positive numeric values")
    }
    lambda_grid <- sort(unique(lambda_grid))
  }
  lambda_fun <- control$lambda_fun %||% NULL
  if (!is.null(lambda_grid) && !is.null(lambda_fun)) {
    stop("Specify at most one of `control$lambda_grid` or `control$lambda_fun`.")
  }

  if (!is.null(blocks)) {
    if (length(blocks) != length(fit$Btil)) {
      stop("`blocks` must have length equal to the number of subjects in `fit`.")
    }
    blocks <- as.character(blocks)
  }

  target_list <- .dkge_prepare_target_list(fit, targets,
                                           residualize = residualize,
                                           collapse = collapse,
                                           restrict_factors = restrict_factors,
                                           scope_override = scope)
  if (!length(target_list)) {
    stop("No valid targets supplied")
  }

  S <- length(fit$Btil)
  subject_ids <- fit$subject_ids %||% seq_len(S)

  resolve_delta_labels <- function(target_idx, target_obj) {
    if (!is.null(target_obj$y)) {
      return(target_obj$y)
    }
    if (is.null(y)) {
      return(NULL)
    }
    if (is.list(y)) {
      if (length(y) == 0) {
        return(NULL)
      }
      name_vec <- names(y)
      if (!is.null(name_vec) && any(nzchar(name_vec))) {
        if (target_obj$name %in% name_vec) {
          return(y[[target_obj$name]])
        }
      }
      if (target_idx <= length(y)) {
        return(y[[target_idx]])
      }
      return(NULL)
    }
    y
  }

  resolve_standardize <- function(target_idx, target_obj, target_mode) {
    default_val <- identical(target_mode, "cell_cross")
    if (is.null(standardize_within_fold)) {
      return(default_val)
    }
    value <- standardize_within_fold
    if (is.list(value)) {
      if (length(value) == 0) {
        return(default_val)
      }
      name_vec <- names(value)
      if (!is.null(name_vec) && any(nzchar(name_vec)) && target_obj$name %in% name_vec) {
        value <- value[[target_obj$name]]
      } else if (target_idx <= length(value)) {
        value <- value[[target_idx]]
      } else {
        return(default_val)
      }
    } else if (length(value) > 1) {
      if (target_idx <= length(value)) {
        value <- value[[target_idx]]
      } else {
        value <- value[[length(value)]]
      }
    }
    value_logical <- suppressWarnings(as.logical(value))
    if (length(value_logical) == 0 || is.na(value_logical[[1]])) {
      return(default_val)
    }
    isTRUE(value_logical[[1]])
  }

  fold_bundle <- .dkge_prepare_folds(fit, folds)
  fold_info <- .dkge_build_fold_bases(fit,
                                      assignments = fold_bundle$assignments,
                                      ridge = ridge,
                                      align = FALSE,
                                      loader_scope = "all",
                                      verbose = verbose)

  results <- vector("list", length(target_list))
  names(results) <- vapply(target_list, `[[`, character(1), "name")

  for (i in seq_along(target_list)) {
    target <- target_list[[i]]
    target_mode <- if (mode == "auto") .dkge_choose_target_mode(target) else mode
    target_scope <- scope %||% target$scope
    if (target_mode %in% c("cell", "cell_cross")) {
      if (nrow(target$weight_matrix) < 2) {
        warning(sprintf("Target '%s' has fewer than two classes; skipping.", target$name))
        next
      }
      standardize_flag <- resolve_standardize(i, target, target_mode)
      res <- .dkge_classify_cell_target(fit, target, fold_info,
                                        mode = target_mode,
                                        method = method,
                                        lambda = lambda,
                                        metric = metric,
                                        class_weights = class_weights,
                                        n_perm = n_perm,
                                        scope = target_scope,
                                        control = list(lambda_grid = lambda_grid, lambda_fun = lambda_fun),
                                        blocks = blocks,
                                        parallel = parallel,
                                        verbose = verbose,
                                        standardize_within_fold = standardize_flag)
    } else if (target_mode == "delta") {
      delta_labels <- resolve_delta_labels(i, target)
      res <- .dkge_classify_delta_target(fit, target, fold_info,
                                         lambda = lambda,
                                         metric = metric,
                                         n_perm = n_perm,
                                         scope = target_scope,
                                         control = list(lambda_fun = lambda_fun),
                                         verbose = verbose,
                                         y = delta_labels,
                                         subject_ids = subject_ids)
    } else {
      warning(sprintf("Target '%s': unsupported mode '%s'; skipping.", target$name, target_mode))
      next
    }
    results[[i]] <- res
  }

  control_out <- if (length(control) == 0) NULL else control

  structure(list(
    fit = fit,
    targets = target_list,
    results = results,
    method = method,
    metric = metric,
    folds = fold_bundle$folds,
    n_perm = n_perm,
    class_weights = class_weights,
    control = control_out,
    blocks = blocks
  ), class = "dkge_classification")
}

.dkge_choose_target_mode <- function(target) {
  if (!is.null(target$class_labels) && length(target$class_labels) >= 2) {
    "cell"
  } else {
    "delta"
  }
}

.dkge_prepare_target_list <- function(fit, targets,
                                      residualize, collapse,
                                      restrict_factors,
                                      scope_override = NULL) {
  if (inherits(targets, "dkge_target")) {
    target_list <- list(targets)
  } else if (is.list(targets) && all(vapply(targets, inherits, logical(1), "dkge_target"))) {
    target_list <- targets
  } else if (is.matrix(targets)) {
    target_list <- list(.dkge_wrap_direct_target(targets))
  } else {
    if (is.null(fit$kernel_info) || is.null(fit$kernel_info$map)) {
      stop(
        "Target specifications that rely on design formulas require `fit$kernel_info$map`. ",
        "For anchor-based fits, supply a weight matrix (rows = classes, columns = effects) ",
        "or explicit `dkge_target` objects (see `dkge_anchor_targets_from_prototypes()`).")
    }
    target_list <- dkge_targets(fit, targets,
                                residualize = residualize,
                                collapse = collapse,
                                restrict_factors = restrict_factors)
  }

  for (i in seq_along(target_list)) {
    if (!inherits(target_list[[i]], "dkge_target")) {
      stop("All targets must inherit from class 'dkge_target'.")
    }
    if (!is.null(scope_override)) {
      target_list[[i]]$scope <- scope_override
    }
    if (is.null(target_list[[i]]$class_labels)) {
      target_list[[i]]$class_labels <- paste0("class", seq_len(nrow(target_list[[i]]$weight_matrix)))
    }
  }
  target_list
}

.dkge_wrap_direct_target <- function(weight_matrix) {
  stopifnot(is.matrix(weight_matrix))
  nm <- "target1"
  target <- list(
    name = nm,
    factors = character(0),
    labels = data.frame(),
    class_labels = if (!is.null(rownames(weight_matrix))) rownames(weight_matrix) else paste0("class", seq_len(nrow(weight_matrix))),
    weight_matrix = weight_matrix,
    indicator = NULL,
    residualized = NA,
    collapse = NULL,
    scope = "within_subject"
  )
  class(target) <- c("dkge_target", "list")
  target
}

.dkge_prepare_folds <- function(fit, folds) {
  .dkge_normalize_folds(folds, fit)
}

#' @export
print.dkge_classification <- function(x, ...) {
  cat("DKGE Classification\n")
  cat("--------------------\n")
  cat(sprintf("Targets: %d\n", length(x$results)))
  cat(sprintf("Classifier: %s\n", x$method))
  cat(sprintf("Metrics: %s\n", paste(x$metric, collapse = ", ")))
  cat(sprintf("Permutations: %d\n", x$n_perm))
  for (nm in names(x$results)) {
    res <- x$results[[nm]]
    if (is.null(res)) {
      cat(sprintf("  %s: <skipped>\n", nm))
      next
    }
    metric_str <- paste(sprintf("%s=%.3f", names(res$metrics), res$metrics), collapse = ", ")
    cat(sprintf("  %s: %s\n", nm, metric_str))
    if (!is.null(res$p_values)) {
      pv_str <- paste(sprintf("%s p=%.3f", names(res$p_values), res$p_values), collapse = ", ")
      cat(sprintf("    %s\n", pv_str))
    }
  }
  invisible(x)
}

#' @export
as.data.frame.dkge_classification <- function(x, row.names = NULL, optional = FALSE, ...,
                                           what = c("summary", "fold_counts", "confusion", "lambda")) {
  stopifnot(inherits(x, "dkge_classification"))
  what <- match.arg(what)
  target_names <- names(x$results)

  get_fold_diag <- function(res) {
    diag <- res$diagnostics$folds
    if (is.null(diag)) list() else diag
  }

  ensure_counts <- function(counts, classes) {
    if (is.null(counts)) {
      setNames(rep(NA_integer_, length(classes)), classes)
    } else {
      setNames(as.integer(counts), classes)
    }
  }

  if (what == "summary") {
    rows <- list()
    for (nm in target_names) {
      res <- x$results[[nm]]
      if (is.null(res)) next
      metrics <- res$metrics
      for (metric_nm in names(metrics)) {
        rows[[length(rows) + 1L]] <- data.frame(
          target = nm,
          metric = metric_nm,
          value = metrics[[metric_nm]],
          p_value = if (!is.null(res$p_values)) res$p_values[[metric_nm]] else NA_real_,
          n_perm = x$n_perm,
          stringsAsFactors = FALSE
        )
      }
    }
    if (!length(rows)) {
      df <- data.frame(target = character(0), metric = character(0), value = numeric(0),
                       p_value = numeric(0), n_perm = integer(0), stringsAsFactors = FALSE)
    } else {
      df <- do.call(rbind, rows)
    }
    if (!is.null(row.names)) rownames(df) <- row.names else rownames(df) <- NULL
    return(df)
  }

  if (what == "fold_counts") {
    rows <- list()
    for (nm in target_names) {
      res <- x$results[[nm]]
      if (is.null(res)) next
      classes <- res$target$class_labels
      diag_list <- get_fold_diag(res)
      if (!length(diag_list)) next
      for (entry in diag_list) {
        if (is.null(entry)) next
        fold_id <- entry$fold %||% NA_integer_
        train_counts <- ensure_counts(entry$class_counts_train, classes)
        test_counts <- ensure_counts(entry$class_counts_test, classes)
        rows[[length(rows) + 1L]] <- data.frame(
          target = nm,
          fold = fold_id,
          class = classes,
          train = as.integer(train_counts[classes]),
          test = as.integer(test_counts[classes]),
          skipped = isTRUE(entry$skipped),
          reason = entry$reason %||% NA_character_,
          stringsAsFactors = FALSE
        )
      }
    }
    if (!length(rows)) {
      df <- data.frame(target = character(0), fold = integer(0), class = character(0),
                       train = integer(0), test = integer(0), skipped = logical(0),
                       reason = character(0), stringsAsFactors = FALSE)
    } else {
      df <- do.call(rbind, rows)
    }
    if (!is.null(row.names)) rownames(df) <- row.names else rownames(df) <- NULL
    return(df)
  }

  if (what == "confusion") {
    rows <- list()
    for (nm in target_names) {
      res <- x$results[[nm]]
      if (is.null(res)) next
      classes <- res$target$class_labels
      diag_list <- get_fold_diag(res)
      if (!length(diag_list)) next
      for (entry in diag_list) {
        if (is.null(entry) || isTRUE(entry$skipped) || is.null(entry$confusion)) next
        fold_id <- entry$fold %||% NA_integer_
        mat <- entry$confusion
        truth_levels <- rownames(mat) %||% classes
        pred_levels <- colnames(mat) %||% classes
        for (truth in truth_levels) {
          counts_row <- mat[truth, , drop = TRUE]
          counts_row <- setNames(as.integer(counts_row), pred_levels)
          rows[[length(rows) + 1L]] <- data.frame(
            target = nm,
            fold = fold_id,
            truth = truth,
            predicted = pred_levels,
            count = counts_row[pred_levels],
            stringsAsFactors = FALSE
          )
        }
      }
    }
    if (!length(rows)) {
      df <- data.frame(target = character(0), fold = integer(0), truth = character(0),
                       predicted = character(0), count = integer(0), stringsAsFactors = FALSE)
    } else {
      df <- do.call(rbind, rows)
    }
    if (!is.null(row.names)) rownames(df) <- row.names else rownames(df) <- NULL
    return(df)
  }

  # what == "lambda"
  rows <- list()
  for (nm in target_names) {
    res <- x$results[[nm]]
    if (is.null(res)) next
    diag_list <- get_fold_diag(res)
    if (!length(diag_list)) next
    for (entry in diag_list) {
      if (is.null(entry)) next
      rows[[length(rows) + 1L]] <- data.frame(
        target = nm,
        fold = entry$fold %||% NA_integer_,
        lambda = entry$lambda %||% NA_real_,
        skipped = isTRUE(entry$skipped),
        reason = entry$reason %||% NA_character_,
        standardized = entry$standardized %||% NA,
        stringsAsFactors = FALSE
      )
    }
  }
  if (!length(rows)) {
    df <- data.frame(target = character(0), fold = integer(0), lambda = numeric(0),
                     skipped = logical(0), reason = character(0), standardized = logical(0),
                     stringsAsFactors = FALSE)
  } else {
    df <- do.call(rbind, rows)
  }
  if (!is.null(row.names)) rownames(df) <- row.names else rownames(df) <- NULL
  df
}

#' Fold-wise confusion matrices for DKGE classification
#'
#' @param x A \code{dkge_classification} object.
#' @param target Optional target name or index. When \code{NULL}, all targets are returned.
#' @param fold Optional fold identifier (index or numeric label). When \code{NULL},
#'   confusion matrices are summed across folds.
#' @return A confusion matrix (when a single target is requested) or a named list of matrices.
#' @export
dkge_confusion <- function(x, target = NULL, fold = NULL) {
  stopifnot(inherits(x, "dkge_classification"))
  target_names <- names(x$results)
  if (is.null(target)) {
    target_idx <- seq_along(target_names)
  } else if (is.numeric(target)) {
    if (any(target < 1) || any(target > length(target_names))) {
      stop("`target` indices out of range.")
    }
    target_idx <- unique(as.integer(target))
  } else {
    match_idx <- match(as.character(target), target_names)
    if (any(is.na(match_idx))) {
      stop("Unknown target specified.")
    }
    target_idx <- unique(match_idx)
  }

  get_fold_diag <- function(res) {
    diag <- res$diagnostics$folds
    if (is.null(diag)) list() else diag
  }

  select_matrices <- function(diag_list, fold) {
    diag_list <- diag_list[!vapply(diag_list, is.null, logical(1))]
    if (!length(diag_list)) {
      return(list())
    }
    available <- Filter(function(entry) {
      !is.null(entry) && !isTRUE(entry$skipped) && !is.null(entry$confusion)
    }, diag_list)
    if (is.null(fold)) {
      return(available)
    }
    fold_vec <- as.vector(fold)
    fold_ids <- vapply(diag_list, function(entry) entry$fold %||% NA_integer_, integer(1))
    select_idx <- integer(0)
    if (is.numeric(fold_vec)) {
      valid_pos <- fold_vec[fold_vec >= 1 & fold_vec <= length(diag_list)]
      select_idx <- c(select_idx, as.integer(valid_pos))
    }
    select_idx <- unique(c(select_idx, which(fold_ids %in% fold_vec)))
    select_idx <- select_idx[select_idx >= 1 & select_idx <= length(diag_list)]
    chosen <- diag_list[select_idx]
    Filter(function(entry) {
      !is.null(entry) && !isTRUE(entry$skipped) && !is.null(entry$confusion)
    }, chosen)
  }

  output <- vector("list", length(target_idx))
  names(output) <- target_names[target_idx]

  for (j in seq_along(target_idx)) {
    res <- x$results[[target_idx[j]]]
    if (is.null(res)) {
      output[[j]] <- NULL
      next
    }
    diag_list <- get_fold_diag(res)
    chosen <- select_matrices(diag_list, fold)
    if (!length(chosen)) {
      output[[j]] <- NULL
      next
    }
    mats <- lapply(chosen, `[[`, "confusion")
    agg <- Reduce(`+`, mats)
    output[[j]] <- agg
  }

  if (length(output) == 1) {
    return(output[[1]])
  }
  output
}
