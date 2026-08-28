# dkge-contrast.R
# Unified contrast engine for DKGE with multiple cross-fitting strategies

#' Compute DKGE contrasts with cross-fitting
#'
#' Main entry point for computing design contrasts after DKGE fit, with support
#' for multiple cross-fitting strategies to ensure unbiased estimation.
#'
#' @param fit A `dkge` object from [dkge_fit()] or [dkge()]
#' @param contrasts Either a q-length numeric vector for a single contrast,
#'   a named list of contrasts, or a qxk matrix where columns are contrasts
#' @param method Cross-fitting strategy: "loso" (leave-one-subject-out),
#'   "kfold" (K-fold cross-validation), or "analytic" (first-order approximation)
#' @param folds For method="kfold", either an integer K for random folds,
#'   or a list defining custom fold assignments (see [dkge_define_folds()])
#' @param ridge Optional ridge parameter added when recomputing held-out basis
#' @param parallel Logical; if TRUE uses `future.apply::future_lapply()` for
#'   per-subject work (requires the future.apply package)
#' @param verbose Logical; print progress messages
#' @param align Logical; if TRUE (default) align LOSO/K-fold bases to a common reference and compute a consensus basis for reporting.
#' @param transport Optional list describing how to transport subject-level
#'   contrasts to a shared reference parcellation. Supply either explicit
#'   transport matrices via `transforms`/`matrices`, or configuration for the
#'   medoid/atlas transport helpers (e.g., `method`, `centroids`, `medoid`). When
#'   provided, the resulting transport bundle is stored under
#'   `metadata$transport` for downstream reuse.
#' @param collinearity_tol Relative angular tolerance used to flag distinct
#'   input contrasts whose kernel-transformed queries are practically
#'   proportional. The default `0.05` flags absolute query correlations of at
#'   least `0.95`; use `NULL` to disable this warning. This is deliberately
#'   separate from the much smaller null-space tolerance used for estimability.
#' @param ... Additional arguments passed to method-specific functions
#'
#' @return A list with class `dkge_contrasts` containing:
#'   - `values`: Named list of contrast values (one P_s vector per subject per contrast)
#'   - `method`: Cross-fitting method used
#'   - `contrasts`: Input contrast specifications
#'   - `metadata`: Method-specific metadata (fold assignments, bases, etc.)
#'
#' @details
#' This function provides a unified interface to three cross-fitting strategies:
#'
#' 1. **LOSO** (`method = "loso"`): Recomputes the basis excluding each subject,
#'    then projects that subject's data. This is the gold standard for unbiased
#'    estimation but requires S eigen-decompositions.
#'
#' 2. **K-fold** (`method = "kfold"`): Splits data into K folds, recomputes basis
#'    excluding each fold, projects held-out data. More efficient than LOSO while
#'    maintaining good bias properties. Supports time-based, run-based, or custom
#'    fold definitions.
#'
#' 3. **Analytic** (`method = "analytic"`): Uses first-order eigenvalue perturbation
#'    theory to approximate the LOSO solution without full recomputation. Fast but
#'    may be less accurate when subjects have high leverage.
#'
#' All methods work entirely in the qxq design space and respect the K-metric
#' throughout. A semidefinite kernel defines a quotient effect space:
#' `dkge_contrast()` errors when a contrast lies wholly in `null(K)`, warns when
#' only part of a contrast is represented, and reports transformed-query
#' collisions in `metadata$kernel_query_pairs`. Multiple contrasts can be
#' evaluated simultaneously for efficiency.
#'
#' Exact fold replay currently supports fits made with `solver = "pooled"` and
#' `cpca_part = "none"`. CPCA and joint-diagonalization fits fail closed for all
#' three methods; DKGE does not replace their fitted estimator with an ordinary
#' pooled eigensolve while calling the result cross-fitted.
#'
#' @examples
#' # Simulate and fit
#' toy <- dkge_sim_toy(
#'   factors = list(A = list(L = 2), B = list(L = 3)),
#'   active_terms = c("A", "B"), S = 4, P = 20, snr = 5
#' )
#' fit <- dkge(toy$B_list, toy$X_list, kernel = toy$K, rank = 2)
#'
#' # Single contrast with LOSO cross-fitting
#' c1 <- c(1, rep(0, 4))
#' result <- dkge_contrast(fit, c1, method = "loso")
#' result
#'
#' # Fast analytic approximation
#' result_fast <- dkge_contrast(fit, c1, method = "analytic")
#'
#' @seealso [dkge_loso_contrast()], [dkge_define_folds()], [dkge_infer()]
#' @export
dkge_contrast <- function(fit, contrasts,
                         method = c("loso", "kfold", "analytic"),
                         folds = NULL,
                         ridge = 0,
                         parallel = FALSE,
                         verbose = FALSE,
                         align = TRUE,
                         transport = NULL,
                         collinearity_tol = 0.05,
                         ...) {
  stopifnot(inherits(fit, "dkge"))
  method <- match.arg(method)
  .dkge_assert_crossfit_estimator_supported(fit, "`dkge_contrast()`")

  # Normalize contrast input
  contrast_list <- .normalize_contrasts(contrasts, fit)
  kernel_contrast_info <- .dkge_validate_kernel_contrasts(
    contrast_list, fit, collinearity_tol = collinearity_tol
  )
  contrast_info <- .dkge_classify_contrasts(contrast_list, fit)
  .dkge_warn_contrast_inference(contrast_info, method)

  # Dispatch to method
  result <- switch(method,
    loso = .dkge_contrast_loso(fit, contrast_list, ridge, parallel, verbose, align = align, ...),
    kfold = .dkge_contrast_kfold(fit, contrast_list, folds, ridge, parallel, verbose, align = align, ...),
    analytic = .dkge_contrast_analytic(fit, contrast_list, ridge, parallel, verbose, align = align, ...)
  )

  if (is.null(result$metadata)) {
    result$metadata <- list()
  }
  result$metadata$contrast_estimability <- contrast_info
  result$metadata$kernel_estimability <- kernel_contrast_info$table
  result$metadata$kernel_query_pairs <- kernel_contrast_info$pairs
  result$metadata$kernel_query_summary <- kernel_contrast_info$summary
  if (is.null(result$metadata$provenance) && !is.null(fit$provenance)) {
    result$metadata$provenance <- fit$provenance
  }
  if (length(result$values) > 0) {
    first_values <- result$values[[1]]
    subject_ids <- names(first_values)
    if (is.null(subject_ids)) {
      subject_ids <- as.character(seq_along(first_values))
    }
    cluster_dims <- vapply(first_values, length, integer(1))
    names(cluster_dims) <- subject_ids
    result$metadata$cluster_dims <- cluster_dims
  } else {
    result$metadata$cluster_dims <- integer(0)
  }

  if (!is.null(transport)) {
    warning("`transport` argument to dkge_contrast() is deprecated; use `dkge_transport_contrasts_to_reference()`.",
            call. = FALSE)
  }

  structure(result, class = "dkge_contrasts")
}

#' Fill in default names for an unnamed or partly named contrast collection
#'
#' Every downstream consumer indexes contrasts by name (`.dkge_classify_contrasts()`
#' builds a `contrast` column, `dkge_contrast_validated()` builds one summary row
#' per name), so a `NULL`/blank name is filled with the same `contrast<j>` label
#' `dkge_contrast()` already uses for unnamed vectors and unnamed matrix columns.
#'
#' @noRd
.dkge_default_contrast_names <- function(nms, n) {
  default <- paste0("contrast", seq_len(n))
  if (is.null(nms)) return(default)
  nms <- as.character(nms)
  length(nms) <- n
  blank <- is.na(nms) | !nzchar(nms)
  nms[blank] <- default[blank]
  nms
}

#' Convert one contrast to a plain numeric vector, keeping its element names
#'
#' `as.numeric()` drops names, which silently discards the effect labels a
#' matrix contrast carries in its row names; those labels are the only handle
#' name-based matching downstream has on the effect ordering.
#'
#' @noRd
.dkge_contrast_vector <- function(x, nms = names(x)) {
  y <- as.numeric(x)
  if (!is.null(nms) && length(nms) == length(y)) names(y) <- as.character(nms)
  y
}

#' Normalize contrast specifications
#'
#' @param contrasts Various input formats
#' @param fit dkge object for validation
#' @return Named list of q-length numeric vectors
#' @keywords internal
#' @noRd
.normalize_contrasts <- function(contrasts, fit) {
  q <- nrow(fit$U)

  if (is.numeric(contrasts) && is.null(dim(contrasts))) {
    # Single vector
    stopifnot(length(contrasts) == q)
    out <- list(contrast1 = .dkge_contrast_vector(contrasts))
    attr(out[[1]], "dkge_scope") <- attr(contrasts, "dkge_scope", exact = TRUE)
    attr(out[[1]], "dkge_term") <- attr(contrasts, "dkge_term", exact = TRUE)
    .dkge_validate_scope_attr(attr(out[[1]], "dkge_scope"))
    return(out)
  }

  if (is.matrix(contrasts)) {
    # Matrix: columns are contrasts
    stopifnot(nrow(contrasts) == q)
    cn <- .dkge_default_contrast_names(colnames(contrasts), ncol(contrasts))
    rn <- rownames(contrasts)
    scope_attr <- attr(contrasts, "dkge_scope", exact = TRUE)
    term_attr <- attr(contrasts, "dkge_term", exact = TRUE)
    pick <- function(a, j) {
      if (is.null(a)) return(NULL)
      if (length(a) == 1L) return(a[[1]])
      if (length(a) >= j) return(a[[j]])
      NULL
    }
    contrast_list <- lapply(seq_len(ncol(contrasts)), function(j) {
      # Carry the matrix row names onto each column: they name the design
      # effects the weights refer to.
      y <- .dkge_contrast_vector(contrasts[, j], rn)
      attr(y, "dkge_scope") <- pick(scope_attr, j)
      attr(y, "dkge_term") <- pick(term_attr, j)
      y
    })
    names(contrast_list) <- cn
    .dkge_validate_scope_attr(scope_attr)
    return(contrast_list)
  }

  if (is.list(contrasts)) {
    # List, named or not
    stopifnot(all(vapply(contrasts, length, integer(1)) == q))
    out <- lapply(contrasts, function(x) {
      y <- .dkge_contrast_vector(x)
      attr(y, "dkge_scope") <- attr(x, "dkge_scope", exact = TRUE)
      attr(y, "dkge_term") <- attr(x, "dkge_term", exact = TRUE)
      .dkge_validate_scope_attr(attr(y, "dkge_scope"))
      y
    })
    names(out) <- .dkge_default_contrast_names(names(contrasts), length(out))
    return(out)
  }

  stop("contrasts must be a numeric vector, matrix, or named list")
}

#' Diagnose whether a kernel can represent planned contrasts
#'
#' Performs a read-only preflight check in the same pooled-design and kernel
#' geometry used by [dkge_contrast()]. This is useful before running
#' cross-fitting or inference, especially when `K` is singular or strongly
#' structured.
#'
#' @param fit A fitted `dkge` object.
#' @param contrasts A contrast vector, matrix, or list accepted by
#'   [dkge_contrast()].
#' @param tol Numerical tolerance used only for null-space estimability checks.
#' @param collinearity_tol Relative angular tolerance used to flag practically
#'   proportional transformed queries. The default `0.05` corresponds to
#'   `abs(query_correlation) >= 0.95`; use `NULL` to disable collision flags.
#' @return A list with `estimability` (support and null fractions per contrast),
#'   `pairs` (pairwise transformed-query correlations and collision flags), a
#'   compact `summary` naming the maximally correlated pair, and scalar `kernel`
#'   rank and spectral-concentration diagnostics.
#' @export
#' @examples
#' toy <- dkge_sim_toy(
#'   factors = list(cond = list(L = 3)), active_terms = "cond",
#'   S = 3, P = 10, snr = 4
#' )
#' fit <- dkge(toy$B_list, toy$X_list, K = toy$K, rank = 2)
#' c_vec <- c(1, -1, rep(0, nrow(fit$U) - 2))
#' dkge_contrast_diagnostics(fit, c_vec)$estimability
dkge_contrast_diagnostics <- function(fit, contrasts, tol = 1e-8,
                                      collinearity_tol = 0.05) {
  stopifnot(inherits(fit, "dkge"))
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) ||
      tol <= 0 || tol >= 1) {
    stop("`tol` must be one finite scalar in (0, 1).", call. = FALSE)
  }
  .dkge_validate_collinearity_tol(collinearity_tol)
  contrast_list <- .normalize_contrasts(contrasts, fit)
  diagnostics <- .dkge_kernel_contrast_diagnostics(
    contrast_list, fit, tol = tol, collinearity_tol = collinearity_tol
  )
  list(
    estimability = diagnostics$table,
    pairs = diagnostics$pairs,
    summary = diagnostics$summary,
    kernel = diagnostics$kernel
  )
}

.dkge_validate_collinearity_tol <- function(collinearity_tol) {
  if (is.null(collinearity_tol)) return(invisible(NULL))
  if (!is.numeric(collinearity_tol) || length(collinearity_tol) != 1L ||
      !is.finite(collinearity_tol) || collinearity_tol <= 0 ||
      collinearity_tol >= 1) {
    stop("`collinearity_tol` must be NULL or one finite scalar in (0, 1).",
         call. = FALSE)
  }
  invisible(NULL)
}

#' Diagnose contrast estimability in a semidefinite kernel geometry
#'
#' A contrast is first moved through the pooled design ruler because that is
#' the coordinate system in which `K` acts. Its Euclidean projection onto
#' image(K) is the estimable part; the orthogonal remainder lies in null(K).
#' The full query available to any DKGE basis is `K^(1/2) R^(-1) c`.
#'
#' @keywords internal
#' @noRd
.dkge_kernel_contrast_diagnostics <- function(contrast_list, fit, tol = 1e-8,
                                               collinearity_tol = 0.05) {
  .dkge_validate_collinearity_tol(collinearity_tol)
  geometry <- .dkge_kernel_geometry(fit$K)
  P <- fit$kernel_support_projector %||% geometry$support_projector
  contrast_names <- names(contrast_list) %||%
    paste0("contrast", seq_along(contrast_list))

  transformed <- lapply(contrast_list, function(ct) {
    backsolve(fit$R, as.numeric(ct), transpose = FALSE)
  })
  support_parts <- lapply(transformed, function(z) as.numeric(P %*% z))
  queries <- lapply(transformed, function(z) as.numeric(geometry$Khalf %*% z))

  rows <- lapply(seq_along(transformed), function(i) {
    total_norm <- sqrt(sum(transformed[[i]]^2))
    support_norm <- sqrt(sum(support_parts[[i]]^2))
    null_norm <- sqrt(sum((transformed[[i]] - support_parts[[i]])^2))
    support_fraction <- if (total_norm > 0) (support_norm / total_norm)^2 else 0
    null_fraction <- if (total_norm > 0) (null_norm / total_norm)^2 else 1
    status <- if (total_norm == 0 || support_norm <= tol * total_norm) {
      "null"
    } else if (null_norm > tol * total_norm) {
      "partially_estimable"
    } else {
      "estimable"
    }
    data.frame(
      contrast = contrast_names[[i]],
      status = status,
      support_fraction = support_fraction,
      null_fraction = null_fraction,
      query_norm = sqrt(sum(queries[[i]]^2)),
      stringsAsFactors = FALSE
    )
  })
  table <- do.call(rbind, rows)

  pair_rows <- list()
  if (length(queries) >= 2L) {
    pairs <- utils::combn(seq_along(queries), 2L)
    pair_rows <- lapply(seq_len(ncol(pairs)), function(j) {
      i1 <- pairs[1L, j]
      i2 <- pairs[2L, j]
      cosine <- function(a, b) {
        denom <- sqrt(sum(a^2) * sum(b^2))
        if (denom == 0) NA_real_ else sum(a * b) / denom
      }
      query_cor <- cosine(queries[[i1]], queries[[i2]])
      raw_cor <- cosine(transformed[[i1]], transformed[[i2]])
      query_defined <- table$status[[i1]] != "null" && table$status[[i2]] != "null"
      if (!query_defined) query_cor <- NA_real_
      collision <- !is.null(collinearity_tol) && query_defined &&
        is.finite(query_cor) && is.finite(raw_cor) &&
        (1 - abs(query_cor)) <= collinearity_tol &&
        (1 - abs(raw_cor)) > collinearity_tol
      data.frame(
        contrast1 = contrast_names[[i1]],
        contrast2 = contrast_names[[i2]],
        query_correlation = query_cor,
        input_correlation = raw_cor,
        collision = collision,
        stringsAsFactors = FALSE
      )
    })
  }
  pair_table <- if (length(pair_rows)) {
    do.call(rbind, pair_rows)
  } else {
    data.frame(
      contrast1 = character(0), contrast2 = character(0),
      query_correlation = numeric(0), input_correlation = numeric(0),
      collision = logical(0), stringsAsFactors = FALSE
    )
  }

  finite_pairs <- which(is.finite(pair_table$query_correlation))
  max_idx <- if (length(finite_pairs)) {
    finite_pairs[[which.max(abs(pair_table$query_correlation[finite_pairs]))]]
  } else {
    NA_integer_
  }
  query_summary <- list(
    n_contrasts = nrow(table),
    n_pairs = nrow(pair_table),
    n_collisions = sum(pair_table$collision, na.rm = TRUE),
    max_abs_query_correlation = if (is.na(max_idx)) NA_real_ else
      abs(pair_table$query_correlation[[max_idx]]),
    max_query_correlation = if (is.na(max_idx)) NA_real_ else
      pair_table$query_correlation[[max_idx]],
    max_query_pair = if (is.na(max_idx)) character(0) else
      c(contrast1 = pair_table$contrast1[[max_idx]],
        contrast2 = pair_table$contrast2[[max_idx]]),
    min_query_norm = if (nrow(table)) min(table$query_norm) else NA_real_,
    max_query_norm = if (nrow(table)) max(table$query_norm) else NA_real_,
    collinearity_tol = collinearity_tol
  )

  list(table = table, pairs = pair_table, summary = query_summary, queries = queries,
       kernel = .dkge_kernel_diagnostics(geometry))
}

#' Enforce the kernel contrast contract
#'
#' @keywords internal
#' @noRd
.dkge_validate_kernel_contrasts <- function(contrast_list, fit, tol = 1e-8,
                                            collinearity_tol = 0.05) {
  diagnostics <- .dkge_kernel_contrast_diagnostics(
    contrast_list, fit, tol = tol, collinearity_tol = collinearity_tol
  )
  null_names <- diagnostics$table$contrast[diagnostics$table$status == "null"]
  if (length(null_names)) {
    .dkge_abort(
      sprintf(
        paste0(
          "Contrast(s) %s lie entirely in null(K) after pooled-design scaling; ",
          "the selected kernel cannot represent these estimands. Choose a ",
          "fuller-rank kernel or revise the contrasts."
        ),
        paste(shQuote(null_names), collapse = ", ")
      ),
      "dkge_kernel_contrast_error"
    )
  }

  partial_names <- diagnostics$table$contrast[
    diagnostics$table$status == "partially_estimable"
  ]
  if (length(partial_names)) {
    .dkge_warn(
      sprintf(
        paste0(
          "Contrast(s) %s contain directions in null(K); DKGE will evaluate ",
          "only their projection onto image(K). Inspect ",
          "`metadata$kernel_estimability`."
        ),
        paste(shQuote(partial_names), collapse = ", ")
      ),
      "dkge_kernel_contrast_warning"
    )
  }

  collisions <- diagnostics$pairs[diagnostics$pairs$collision, , drop = FALSE]
  if (nrow(collisions)) {
    labels <- apply(collisions[c("contrast1", "contrast2")], 1L, function(x) {
      paste(shQuote(x), collapse = " / ")
    })
    .dkge_warn(
      sprintf(
        paste0(
          "Distinct contrast pairs produce nearly proportional kernel queries: %s. ",
          "Their DKGE maps may be practically indistinguishable up to scale or ",
          "sign; inspect ",
          "`metadata$kernel_query_pairs`."
        ),
        paste(labels, collapse = "; ")
      ),
      "dkge_kernel_contrast_collision_warning"
    )
  }
  diagnostics
}

.dkge_contrast_recommendation <- function(estimability) {
  if (identical(estimability, "within")) {
    return("LOSO/k-fold cross-fitting")
  }
  if (estimability %in% c("between", "mixed")) {
    return("subject-label permutation or dkge_between_* inference")
  }
  "unknown; inspect contrast scope"
}

#' Structural estimability scope for a contrast over design cells
#'
#' Classifies a contrast by *which design factors it actually varies over*,
#' rather than by matching its name against kernel term names. For each factor
#' the contrast weights are grouped by the levels of all other factors; the
#' contrast depends on that factor when the weights differ within at least one
#' such group. The scope is then `"between"` when only between-scope factors are
#' involved, `"within"` when only within-scope factors are, and `"mixed"` when
#' both are. A contrast with constant weights (a grand mean) depends on no
#' factor and is reported as `"within"`, since every subject can estimate it.
#'
#' Returns `NULL` when the fit carries no usable cell metadata (for example an
#' effect-basis kernel, whose coordinates are not design cells), leaving the
#' caller to fall back to other evidence.
#'
#' @keywords internal
#' @noRd
.dkge_contrast_structural_scope <- function(contrast, fit, tol = 1e-10) {
  info <- fit$kernel_info
  if (is.null(info)) return(NULL)
  factor_scope <- info$factor_scope
  if (is.null(factor_scope) || is.null(names(factor_scope))) {
    return(NULL)
  }
  cells <- .dkge_match_kernel_cells(fit, info)
  cvec <- as.numeric(contrast)
  if (is.null(cells) || !nrow(cells) || nrow(cells) != length(cvec) || anyNA(cvec)) {
    return(NULL)
  }
  factor_names <- intersect(names(factor_scope), names(cells))
  if (!length(factor_names)) return(NULL)

  scale <- max(abs(cvec), 1)
  depends <- vapply(factor_names, function(nm) {
    others <- setdiff(factor_names, nm)
    key <- if (length(others)) {
      do.call(paste, c(lapply(others, function(o) as.character(cells[[o]])), sep = "\r"))
    } else {
      rep("", nrow(cells))
    }
    spreads <- tapply(cvec, key, function(z) max(z) - min(z))
    any(spreads > tol * scale)
  }, logical(1))
  names(depends) <- factor_names

  scopes <- as.character(factor_scope[factor_names])
  between_factors <- factor_names[scopes == "between"]
  within_factors <- setdiff(factor_names, between_factors)
  uses_between <- any(depends[between_factors])
  uses_within <- any(depends[within_factors])

  if (uses_between && uses_within) return("mixed")
  if (uses_between) return("between")
  "within"
}

.dkge_allowed_scope <- c("within", "between", "mixed")

.dkge_validate_scope_attr <- function(scope) {
  if (is.null(scope)) return(invisible(NULL))
  scope <- as.character(scope)
  bad <- setdiff(unique(scope[!is.na(scope)]), .dkge_allowed_scope)
  if (length(bad)) {
    stop(sprintf(
      "`dkge_scope` must be one of %s; got %s.",
      paste(shQuote(.dkge_allowed_scope), collapse = ", "),
      paste(shQuote(bad), collapse = ", ")
    ), call. = FALSE)
  }
  invisible(scope)
}

.dkge_classify_one_contrast <- function(name, contrast, fit, tol = 1e-10) {
  scope_attr <- attr(contrast, "dkge_scope", exact = TRUE)
  if (!is.null(scope_attr)) {
    .dkge_validate_scope_attr(scope_attr)
    return(as.character(scope_attr)[[1]])
  }

  term_attr <- attr(contrast, "dkge_term", exact = TRUE)
  term_scope <- fit$kernel_info$term_scope %||% NULL
  if (!is.null(term_attr) && !is.null(term_scope) && term_attr %in% names(term_scope)) {
    return(unname(term_scope[[term_attr]]))
  }

  # Structural evidence (does the contrast vary across levels of a
  # between-subject factor?) is authoritative when available. A contrast NAME
  # that happens to match a kernel term is only a hint: names are arbitrary
  # user labels, so a between-subject vector named "task" must not be
  # classified as within. Where the two disagree, take the structural answer.
  structural <- .dkge_contrast_structural_scope(contrast, fit, tol = tol)
  if (!is.null(structural)) {
    return(structural)
  }
  if (!is.null(term_scope) && !is.null(name) && name %in% names(term_scope)) {
    return(unname(term_scope[[name]]))
  }

  blocks <- fit$kernel_info$blocks %||% NULL
  if (!is.null(blocks) && !is.null(term_scope)) {
    active <- names(blocks)[vapply(blocks, function(idx) {
      any(abs(contrast[idx]) > tol)
    }, logical(1))]
    active <- intersect(active, names(term_scope))
    if (length(active)) {
      scopes <- unname(term_scope[active])
      if (all(scopes == "within")) {
        return("within")
      }
      if (all(scopes == "between")) {
        return("between")
      }
      return("mixed")
    }
  }

  "unknown"
}

.dkge_classify_contrasts <- function(contrast_list, fit) {
  contrast_names <- names(contrast_list) %||% paste0("contrast", seq_along(contrast_list))
  estimability <- vapply(seq_along(contrast_list), function(i) {
    .dkge_classify_one_contrast(contrast_names[[i]], contrast_list[[i]], fit)
  }, character(1))
  data.frame(
    contrast = contrast_names,
    estimability = unname(estimability),
    recommended_inference = unname(vapply(estimability, .dkge_contrast_recommendation, character(1))),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

.dkge_warn_contrast_inference <- function(contrast_info, method) {
  if (!method %in% c("loso", "kfold", "analytic")) {
    return(invisible(NULL))
  }
  needs_between <- contrast_info$estimability %in% c("between", "mixed")
  if (any(needs_between)) {
    bad <- contrast_info$contrast[needs_between]
    warning(
      sprintf(
        "Contrast(s) %s are between/mixed effects; %s results are descriptive. Use subject-label permutation or dkge_between_* inference for group-effect testing.",
        paste(shQuote(bad), collapse = ", "),
        method
      ),
      call. = FALSE
    )
  }
  invisible(NULL)
}

#' LOSO contrast implementation
#'
#' @inheritParams dkge_contrast
#' @param contrast_list Normalized list of contrasts
#' @keywords internal
#' @noRd
.dkge_contrast_loso <- function(fit, contrast_list, ridge, parallel, verbose, align, ...) {
  S <- length(fit$Btil)
  n_contrasts <- length(contrast_list)
  verbose_flag <- .dkge_verbose(verbose)

  if (verbose_flag) {
    message(sprintf("Computing %d contrast(s) via LOSO for %d subjects", n_contrasts, S))
  }

  subject_labels <- fit$subject_ids %||% paste0("subject", seq_len(S))
  assignments <- lapply(seq_len(S), function(s) s)

  fold_info <- .dkge_build_fold_bases(
    fit,
    assignments = assignments,
    ridge = ridge,
    align = align,
    loader_scope = "heldout",
    verbose = verbose
  )
  folds <- fold_info$folds

  c_tilde_list <- lapply(contrast_list, function(ct) backsolve(fit$R, ct, transpose = FALSE))

  values <- vector("list", n_contrasts)
  names(values) <- names(contrast_list)
  alphas <- vector("list", n_contrasts)
  names(alphas) <- names(contrast_list)

  r <- ncol(fit$U)
  fold_row_names <- vapply(folds, function(fold) paste(subject_labels[fold$subjects], collapse = ","), character(1))

  for (i in seq_along(contrast_list)) {
    values[[i]] <- vector("list", S)
    names(values[[i]]) <- subject_labels
    alpha_mat <- matrix(NA_real_, nrow = length(folds), ncol = r)
    rownames(alpha_mat) <- fold_row_names

    for (fold in folds) {
      U_fold <- fold$basis
      alpha_vec <- as.numeric(t(U_fold) %*% fit$K %*% c_tilde_list[[i]])
      alpha_mat[fold$index, seq_along(alpha_vec)] <- alpha_vec

      holdout <- fold$subjects
      loaders <- fold$loaders
      subject_scores <- .dkge_apply(
        holdout,
        function(s) {
          loader <- loaders[[as.character(s)]]
          as.numeric(loader$A %*% alpha_vec)
        },
        parallel = parallel
      )

      for (pos in seq_along(holdout)) {
        s <- holdout[[pos]]
        values[[i]][[subject_labels[[s]]]] <- subject_scores[[pos]]
      }
    }

    alphas[[i]] <- alpha_mat
  }

  for (i in seq_along(values)) {
    missing <- which(vapply(values[[i]], is.null, logical(1)))
    if (length(missing) > 0) {
      warning(sprintf("Subjects %s not processed for contrast %s",
                      paste(subject_labels[missing], collapse = ","),
                      names(values)[i] %||% paste0("contrast", i)))
    }
  }

  metadata <- list(
    bases = lapply(folds, `[[`, "basis"),
    aligned_bases = lapply(folds, `[[`, "basis_aligned"),
    rotations = lapply(folds, `[[`, "rotation"),
    alphas = alphas,
    alignment_receipts = .dkge_alignment_receipts_from_folds(
      fit, fold_info, alphas, method = "loso"
    ),
    ridge = ridge,
    procrustes = if (align) list(alignment = fold_info$alignment, consensus = fold_info$consensus) else NULL
  )

  list(
    values = values,
    method = "loso",
    contrasts = contrast_list,
    metadata = metadata
  )
}

#' Analytic LOSO contrast implementation
#'
#' @inheritParams .dkge_contrast_loso
#' @keywords internal
#' @noRd
.dkge_contrast_analytic <- function(fit, contrast_list, ridge, parallel, verbose, align = TRUE, ...) {
  .dkge_contrast_analytic_impl(fit, contrast_list, ridge, parallel, verbose, align = align, ...)
}

#' Print method for dkge_contrasts
#'
#' @param x A dkge_contrasts object
#' @param ... Additional arguments (unused)
#' @export
print.dkge_contrasts <- function(x, ...) {
  n_contrasts <- length(x$contrasts)
  n_subjects <- length(x$values[[1]])

  cat("DKGE Contrasts\n")
  cat("--------------\n")
  cat(sprintf("Method: %s\n", x$method))
  cat(sprintf("Contrasts: %d\n", n_contrasts))
  cat(sprintf("Subjects: %d\n", n_subjects))

  if (n_contrasts <= 5) {
    cat("Contrast names:", paste(names(x$contrasts), collapse = ", "), "\n")
  } else {
    cat("Contrast names:", paste(names(x$contrasts)[1:5], collapse = ", "), "...\n")
  }

  if (!is.null(x$metadata$ridge) && x$metadata$ridge > 0) {
    cat(sprintf("Ridge: %g\n", x$metadata$ridge))
  }

  detail <- x$metadata$fallback_detail
  if (!is.null(detail) && nrow(detail)) {
    fallback_rows <- detail[is.na(detail$reason) | detail$reason != "analytic", , drop = FALSE]
    fallback_count <- nrow(fallback_rows)
    if (fallback_count > 0) {
      cat(sprintf("Fallback triggered for %d/%d subject-contrast pairs (see metadata$fallback_detail).\n",
                  fallback_count, nrow(detail)))
      top_rows <- head(fallback_rows, 5)
      for (k in seq_len(nrow(top_rows))) {
        row <- top_rows[k, , drop = FALSE]
        cat(sprintf("    subject=%s, contrast=%s, reason=%s\n",
                    row$subject, row$contrast, row$reason))
      }
      if (fallback_count > nrow(top_rows)) {
        cat("    ...\n")
      }
    }
  }

  invisible(x)
}

#' Extract contrast values as matrix
#'
#' @param x A dkge_contrasts object
#' @param contrast Name or index of contrast to extract
#' @param ... Additional arguments (not used)
#' @return SxP matrix of contrast values
#' @export
as.matrix.dkge_contrasts <- function(x, contrast = 1, ...) {
  if (is.character(contrast)) {
    contrast <- match(contrast, names(x$contrasts))
    if (is.na(contrast)) stop("Contrast not found")
  }

  value_list <- x$values[[contrast]]
  dims <- vapply(value_list, length, integer(1))
  unique_dims <- unique(dims)

  if (length(unique_dims) != 1) {
    msg <- "Subject cluster counts differ; use dkge_transport_contrasts_to_reference() before stacking."
    cond <- structure(list(message = msg, call = sys.call()),
                     class = c("dkge_transport_needed", "error", "condition"))
    stop(cond)
  }

  do.call(rbind, value_list)
}

#' @export
as.data.frame.dkge_contrasts <- function(x, row.names = NULL, optional = FALSE, ...,
                                         stringsAsFactors = FALSE) {
  contrast_names <- names(x$contrasts)
  if (is.null(contrast_names) || any(!nzchar(contrast_names))) {
    contrast_names <- paste0("contrast", seq_along(x$contrasts))
  }

  rows <- vector("list", length(x$values))
  for (i in seq_along(x$values)) {
    subj_values <- x$values[[i]]
    if (!length(subj_values)) {
      next
    }
    subj_names <- names(subj_values)
    if (is.null(subj_names) || any(!nzchar(subj_names))) {
      subj_names <- paste0("subject", seq_along(subj_values))
    }

    entries <- Map(function(vals, subj) {
      cluster_ids <- names(vals)
      if (is.null(cluster_ids) || any(!nzchar(cluster_ids))) {
        cluster_ids <- paste0("cluster", seq_along(vals))
      }
      data.frame(
        contrast = contrast_names[[i]],
        subject = subj,
        component = cluster_ids,
        value = as.numeric(vals),
        method = x$method,
        stringsAsFactors = stringsAsFactors
      )
    }, subj_values, subj_names)
    rows[[i]] <- do.call(rbind, entries)
  }

  rows <- Filter(Negate(is.null), rows)
  result <- if (length(rows)) do.call(rbind, rows) else NULL
  if (is.null(result)) {
    result <- data.frame(
      contrast = character(0),
      subject = character(0),
      component = character(0),
      value = numeric(0),
      method = character(0),
      stringsAsFactors = stringsAsFactors
    )
  }
  if (!is.null(row.names)) {
    rownames(result) <- row.names
  } else {
    rownames(result) <- NULL
  }
  result
}
