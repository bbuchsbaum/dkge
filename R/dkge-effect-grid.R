# dkge-effect-grid.R
# Lightweight metadata for global effect grids with within/between scopes.

#' Construct a DKGE global effect grid
#'
#' @param factors Named list of factor specifications. Each entry may be a
#'   design-kernel factor spec, a vector of level labels, or a scalar level
#'   count.
#' @param scope Named character vector or list assigning each factor to
#'   `"within"` or `"between"`. Defaults to `"within"` for all factors. Scope
#'   is term metadata used by contrast routing; it does not alter the kernel.
#' @param block_factors Optional factor names that should define independent
#'   blocks when passed to [design_kernel()].
#' @param sep Separator used when constructing cell labels.
#' @return A `dkge_effect_grid` object with factor specs, scopes, cells, and
#'   labels.
#' @details Passing the resulting grid to [design_kernel()] preserves its
#'   globally named cell order. For three or more factors, note that
#'   `design_kernel(terms = NULL)` includes the main effects and the single full
#'   interaction only. Supply `terms` explicitly when two-way interactions are
#'   part of the scientific model.
#' @export
dkge_effect_grid <- function(factors,
                             scope = NULL,
                             block_factors = NULL,
                             sep = ":") {
  stopifnot(is.list(factors), length(factors) > 0)
  factor_names <- names(factors)
  if (is.null(factor_names) || any(!nzchar(factor_names)) || any(duplicated(factor_names))) {
    stop("`factors` must be a named list with unique names.", call. = FALSE)
  }

  specs <- Map(function(x, nm) {
    if (is.list(x) && !is.null(x$L)) {
      spec <- x
      if (is.null(spec$levels)) {
        spec$levels <- paste0(nm, seq_len(as.integer(spec$L)))
      }
      spec$L <- as.integer(spec$L)
      return(spec)
    }
    if (is.numeric(x) && length(x) == 1L) {
      L <- as.integer(x)
      return(list(L = L, type = "nominal", levels = paste0(nm, seq_len(L))))
    }
    lev <- as.character(x)
    list(L = length(lev), type = "nominal", levels = lev)
  }, factors, factor_names)
  names(specs) <- factor_names

  scopes <- setNames(rep("within", length(factor_names)), factor_names)
  if (!is.null(scope)) {
    scope <- unlist(scope, use.names = TRUE)
    if (is.null(names(scope)) || !all(names(scope) %in% factor_names)) {
      stop("`scope` must be named with factor names.", call. = FALSE)
    }
    scopes[names(scope)] <- as.character(scope)
  }
  if (!all(scopes %in% c("within", "between"))) {
    stop("`scope` entries must be either 'within' or 'between'.", call. = FALSE)
  }

  if (!is.null(block_factors)) {
    block_factors <- as.character(block_factors)
    missing <- setdiff(block_factors, factor_names)
    if (length(missing)) {
      stop("Unknown block factor(s): ", paste(missing, collapse = ", "), call. = FALSE)
    }
  }

  level_list <- lapply(specs, `[[`, "levels")
  rev_cells <- do.call(expand.grid, c(rev(level_list),
                                      KEEP.OUT.ATTRS = FALSE,
                                      stringsAsFactors = FALSE))
  cells <- rev_cells[rev(seq_along(rev_cells))]
  names(cells) <- factor_names
  labels <- do.call(paste, c(cells, sep = sep))

  out <- list(
    factors = specs,
    factor_names = factor_names,
    scope = scopes,
    block_factors = block_factors,
    cells = cells,
    cell_labels = labels,
    sep = sep
  )
  class(out) <- "dkge_effect_grid"
  out
}

.dkge_effect_grid_cell_labels <- function(factors, sep = ":") {
  levels <- lapply(names(factors), function(nm) {
    f <- factors[[nm]]
    if (!is.null(f$levels)) as.character(f$levels) else paste0(nm, seq_len(f$L))
  })
  names(levels) <- names(factors)
  rev_cells <- do.call(expand.grid, c(rev(levels),
                                      KEEP.OUT.ATTRS = FALSE,
                                      stringsAsFactors = FALSE))
  cells <- rev_cells[rev(seq_along(rev_cells))]
  names(cells) <- names(factors)
  list(cells = cells, labels = do.call(paste, c(cells, sep = sep)))
}

.dkge_term_scope <- function(term, factor_scope) {
  if (!length(term)) {
    return("within")
  }
  scopes <- factor_scope[term]
  if (all(scopes == "between")) {
    "between"
  } else if (any(scopes == "between")) {
    "mixed"
  } else {
    "within"
  }
}
