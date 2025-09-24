# dkge-targets.R
# Target construction utilities for DKGE classification.

#' Build classification targets from a DKGE fit
#'
#' Generates a list of target specifications that map cell-level design
#' structure to the DKGE effect space. Each target provides a weight matrix that
#' can be multiplied with subject effect coefficients to obtain class-specific
#' patterns for decoding.
#'
#' @param fit dkge object containing `kernel_info$map` metadata.
#' @param spec Target specification. Accepts a formula (e.g. `~ A + B + A:B`),
#'   a character vector of term labels, the string "fullcell", or an existing
#'   list of `dkge_target` objects (returned unchanged).
#' @param residualize Logical; if `TRUE` (default) residualise higher-order
#'   targets against previously constructed lower-order targets.
#' @param collapse Optional named list describing how to collapse factors that do
#'   not appear in a target. Each entry may be `"mean"`,
#'   `list(method = "mean", window = 3:8)`, or a numeric vector of length equal
#'   to the number of levels providing custom weights (automatically normalised).
#' @param scope Permutation/exchangeability scope stored with each target
#'   (default "within_subject").
#' @param restrict_factors Optional character vector restricting the factors used
#'   when `spec = "fullcell"`. When `NULL`, all factors are used.
#'
#' @return List of objects with class `dkge_target`. Each target contains
#'   `name`, `factors`, `labels`, `weight_matrix`, `indicator`, `residualized`,
#'   `collapse`, and `scope` fields.
#' @export
#' @examples
#' \dontrun{
#' tg <- dkge_targets(fit, ~ A + B + A:B, collapse = list(time = list(method = "mean", window = 3:8)))
#' }
dkge_targets <- function(fit,
                         spec = NULL,
                         residualize = TRUE,
                         collapse = NULL,
                         scope = "within_subject",
                         restrict_factors = NULL) {
  stopifnot(inherits(fit, "dkge"))

  if (is.list(spec) && all(vapply(spec, inherits, logical(1), what = "dkge_target"))) {
    return(spec)
  }

  info <- fit$kernel_info
  if (is.null(info) || is.null(info$map)) {
    stop("fit does not contain kernel mapping information; supply custom targets or rebuild with design_kernel(..., basis = 'effect').")
  }

  factor_levels <- info$levels
  if (is.null(factor_levels) || is.null(names(factor_levels))) {
    stop("kernel_info$levels not available; cannot derive targets from design factors.")
  }
  factor_names <- names(factor_levels)
  map <- info$map

  if (is.null(spec)) {
    spec <- stats::as.formula(paste("~", paste(factor_names, collapse = " + ")))
  }

  target_terms <- .dkge_parse_target_spec(spec, factor_names, restrict_factors)
  if (!length(target_terms)) {
    stop("No target terms recognised in `spec`.")
  }

  collapse_weights <- .dkge_parse_collapse_rules(factor_levels, collapse)

  targets <- list()
  lower_basis <- NULL

  for (term in target_terms) {
    indicator <- .dkge_indicator_from_term(term$factors, factor_levels, factor_names, collapse_weights)
    indicator <- .dkge_normalize_rows(indicator)

    if (residualize && !is.null(lower_basis) && nrow(lower_basis) > 0 && length(term$factors) > 0) {
      indicator <- .dkge_residualize_rows(indicator, lower_basis)
    }

    keep <- which(.dkge_row_norm2(indicator) > 1e-10)
    if (!length(keep)) {
      warning(sprintf("Target '%s' eliminated because all rows vanished after residualisation.", term$name))
      next
    }
    indicator <- indicator[keep, , drop = FALSE]

    labels_df <- .dkge_target_labels(term$factors, factor_levels)[keep, , drop = FALSE]
    class_labels <- .dkge_collapse_label_rows(labels_df)

    W <- indicator %*% map
    keep_w <- which(.dkge_row_norm2(W) > 1e-10)
    if (!length(keep_w)) {
      warning(sprintf("Target '%s' eliminated because weight rows are numerically zero.", term$name))
      next
    }
    if (length(keep_w) < nrow(W)) {
      indicator <- indicator[keep_w, , drop = FALSE]
      labels_df <- labels_df[keep_w, , drop = FALSE]
      class_labels <- class_labels[keep_w]
      W <- W[keep_w, , drop = FALSE]
    }

    target <- list(
      name = term$name,
      factors = term$factors,
      labels = labels_df,
      class_labels = class_labels,
      weight_matrix = W,
      indicator = indicator,
      residualized = residualize,
      collapse = collapse,
      scope = scope
    )
    class(target) <- c("dkge_target", "list")
    targets[[length(targets) + 1L]] <- target

    lower_basis <- rbind(lower_basis, indicator)
  }

  targets
}

#' Target helper for a single factor
#'
#' @inheritParams dkge_targets
#' @param factor Factor name.
#' @export
dkge_target_factor <- function(fit, factor,
                               residualize = TRUE,
                               collapse = NULL,
                               scope = "within_subject") {
  stopifnot(length(factor) == 1L)
  term <- stats::as.formula(paste("~", factor))
  tg <- dkge_targets(fit, term, residualize = residualize,
                     collapse = collapse, scope = scope)
  if (length(tg)) tg[[1]] else tg
}

#' Target helper for an interaction
#'
#' @inheritParams dkge_targets
#' @param factors Character vector of factor names to include in the interaction.
#' @export
dkge_target_interaction <- function(fit, factors,
                                    residualize = TRUE,
                                    collapse = NULL,
                                    scope = "within_subject") {
  stopifnot(is.character(factors), length(factors) >= 2)
  term <- stats::as.formula(paste("~", paste(factors, collapse = ":")))
  tg <- dkge_targets(fit, term, residualize = residualize,
                     collapse = collapse, scope = scope)
  if (length(tg)) tg[[1]] else tg
}

# Helper: parse target specification into ordered terms
.dkge_parse_target_spec <- function(spec, factor_names, restrict_factors = NULL) {
  if (is.character(spec) && length(spec) == 1 && identical(spec, "fullcell")) {
    if (is.null(restrict_factors)) {
      restrict_factors <- factor_names
    }
    restrict_factors <- intersect(restrict_factors, factor_names)
    if (!length(restrict_factors)) {
      stop("`restrict_factors` did not match any design factors.")
    }
    term <- list(name = paste(restrict_factors, collapse = ":"),
                 factors = restrict_factors,
                 order = length(restrict_factors))
    return(list(term))
  }

  term_labels <- NULL
  if (inherits(spec, "formula")) {
    term_obj <- stats::terms(spec, simplify = TRUE)
    term_labels <- attr(term_obj, "term.labels")
  } else if (is.character(spec)) {
    term_labels <- spec
  } else {
    stop("Unsupported target specification. Use a formula, character vector, or 'fullcell'.")
  }

  term_labels <- unique(term_labels)
  terms <- vector("list", length(term_labels))
  for (i in seq_along(term_labels)) {
    lbl <- term_labels[[i]]
    parts <- strsplit(lbl, ":", fixed = TRUE)[[1]]
    parts <- trimws(parts)
    if (!all(parts %in% factor_names)) {
      stop(sprintf("Target term '%s' references unknown factors.", lbl))
    }
    terms[[i]] <- list(name = lbl, factors = parts, order = length(parts))
  }

  # Sort by interaction order (ascending), preserve user order within ties
  order_idx <- order(vapply(terms, `[[`, integer(1), "order"))
  terms[order_idx]
}

# Helper: build collapse weights per factor
.dkge_parse_collapse_rules <- function(levels, collapse) {
  weights <- lapply(levels, function(L) rep(1, L))
  if (is.null(collapse)) {
    return(weights)
  }
  for (nm in names(collapse)) {
    if (!nm %in% names(levels)) {
      warning(sprintf("Collapse rule provided for unknown factor '%s'; ignoring.", nm))
      next
    }
    L <- levels[[nm]]
    rule <- collapse[[nm]]
    vec <- rep(1, L)
    if (is.character(rule) && length(rule) == 1) {
      if (identical(rule, "mean")) {
        vec <- rep(1, L)
      } else {
        warning(sprintf("Unsupported collapse method '%s' for factor '%s'.", rule, nm))
      }
    } else if (is.list(rule)) {
      method <- tolower(rule$method %||% "mean")
      if (method == "mean") {
        window <- rule$window %||% seq_len(L)
        window <- intersect(window, seq_len(L))
        vec <- rep(0, L)
        vec[window] <- 1
      } else {
        warning(sprintf("Unsupported collapse method '%s' for factor '%s'.", method, nm))
      }
    } else if (is.numeric(rule)) {
      if (length(rule) != L) {
        stop(sprintf("Collapse weights for factor '%s' must have length %d.", nm, L))
      }
      vec <- as.numeric(rule)
    } else {
      stop(sprintf("Unsupported collapse specification for factor '%s'.", nm))
    }
    if (sum(vec) > 0) vec <- vec / sum(vec)
    weights[[nm]] <- vec
  }
  weights
}

# Helper: build Kronecker indicator matrix for a term
.dkge_indicator_from_term <- function(term_factors, levels, factor_names, collapse_weights) {
  mats <- vector("list", length(factor_names))
  names(mats) <- factor_names
  for (nm in factor_names) {
    L <- levels[[nm]]
    if (nm %in% term_factors) {
      mats[[nm]] <- diag(L)
    } else {
      w <- collapse_weights[[nm]]
      if (is.null(w)) {
        w <- rep(1 / L, L)
      } else if (abs(sum(w)) > 0) {
        # ensure weights sum to one for averaging
        w <- w / sum(w)
      }
      mats[[nm]] <- matrix(w, nrow = 1)
    }
  }
  Reduce(kronecker, mats)
}

# Helper: expand class labels for a term
.dkge_target_labels <- function(term_factors, levels) {
  if (!length(term_factors)) {
    return(data.frame())
  }
  lvl_list <- lapply(term_factors, function(nm) seq_len(levels[[nm]]))
  labels <- expand.grid(lvl_list, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  colnames(labels) <- term_factors
  labels
}

# Helper: residualise rows of M against row space of L
.dkge_residualize_rows <- function(M, L) {
  if (is.null(L) || nrow(L) == 0) {
    return(M)
  }
  gram <- L %*% t(L)
  diag(gram) <- diag(gram) + 1e-8
  coeff <- M %*% t(L)
  adjust <- coeff %*% solve(gram, L)
  M - adjust
}

.dkge_normalize_rows <- function(M) {
  if (!nrow(M)) return(M)
  rs <- rowSums(M)
  idx <- which(rs > 0)
  if (length(idx)) {
    M[idx, ] <- M[idx, , drop = FALSE] / rs[idx]
  }
  M
}

.dkge_row_norm2 <- function(M) rowSums(M * M)

.dkge_collapse_label_rows <- function(labels_df) {
  if (!nrow(labels_df)) return(character(0))
  apply(labels_df, 1, function(row) paste(row, collapse = ":"))
}
