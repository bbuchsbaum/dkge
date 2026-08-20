# dkge-effect-grid.R
# Lightweight metadata for global effect grids with within/between scopes.

# --- Shared factor-specification validators -------------------------------
# Used by both dkge_effect_grid() and design_kernel() so the two constructors
# accept and reject exactly the same specs. Every guard is length-safe: a
# vector-valued `L` or `type` used to reach an `if()` and abort with R's
# "the condition has length > 1" instead of a message naming the factor.

#' Coerce a factor's `L` field to a single positive integer
#'
#' @noRd
.dkge_factor_level_count <- function(L, nm) {
  if (is.null(L) || length(L) != 1L) {
    stop(sprintf(
      "Factor '%s': `L` must be a single level count (got length %d).",
      nm, length(L)
    ), call. = FALSE)
  }
  if (!is.numeric(L) && !is.logical(L) && !is.character(L)) {
    stop(sprintf("Factor '%s': `L` must be a single level count.", nm),
         call. = FALSE)
  }
  L <- suppressWarnings(as.integer(L))
  if (is.na(L) || L < 1L) {
    stop(sprintf("Factor '%s' must have at least one level.", nm), call. = FALSE)
  }
  L
}

#' Validate a factor's `type` field
#'
#' @noRd
.dkge_factor_type <- function(type, nm) {
  type <- type %||% "nominal"
  if (!is.character(type) || length(type) != 1L || is.na(type)) {
    stop(sprintf(
      "Factor '%s': `type` must be a single character string (one of nominal, ordinal, circular, continuous).",
      nm
    ), call. = FALSE)
  }
  type <- tolower(type)
  if (!type %in% c("nominal", "ordinal", "circular", "continuous")) {
    stop(sprintf(
      "Factor '%s': unsupported type '%s' (use nominal, ordinal, circular, or continuous).",
      nm, type
    ), call. = FALSE)
  }
  type
}

#' Validate a factor's level labels
#'
#' Duplicated labels within a factor make the expanded cell labels ambiguous,
#' so the kernel could not be indexed by name afterwards.
#'
#' @noRd
.dkge_factor_levels <- function(levels, nm) {
  levels <- as.character(levels)
  if (anyNA(levels) || any(!nzchar(levels))) {
    stop(sprintf("Factor '%s': level labels must be non-empty strings.", nm),
         call. = FALSE)
  }
  if (anyDuplicated(levels)) {
    dup <- unique(levels[duplicated(levels)])
    stop(sprintf(
      "Factor '%s' has duplicated level label(s): %s. Cell labels must be unique.",
      nm, paste(dup, collapse = ", ")
    ), call. = FALSE)
  }
  levels
}

#' Construct a DKGE global effect grid
#'
#' @param factors Named list of factor specifications. Each entry may be a
#'   design-kernel factor spec, a vector of level labels, or a scalar level
#'   count. List specs accept `type`, `L`, `levels`, `values`, and `l`
#'   (positive finite length-scale for ordinal/circular/continuous factors).
#' @param scope Named character vector or list assigning each factor to
#'   `"within"` or `"between"`. Defaults to `"within"` for all factors.
#' @param block_factors Optional factor names that should define independent
#'   blocks when passed to [design_kernel()].
#' @param sep Separator used when constructing cell labels.
#' @return A `dkge_effect_grid` object with factor specs, scopes, cells, labels,
#'   and the default kernel terms implied by the factor order. A one-factor
#'   grid has one default main-effect term because its full interaction is the
#'   same term.
#' @export
#' @examples
#' grid <- dkge_effect_grid(
#'   factors = list(
#'     cue = c("valid", "invalid"),
#'     load = list(L = 3, type = "ordinal"),
#'     group = 2
#'   ),
#'   scope = c(group = "between"),
#'   block_factors = "group"
#' )
#' grid$cell_labels
#' grid$scope
#'
#' # Feed the grid to design_kernel() so the kernel indexes the same cells
#' kern <- design_kernel(grid, basis = "cell")
#' identical(rownames(kern$K), grid$cell_labels)
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
    if (is.list(x)) {
      spec <- x
      known_fields <- c("type", "L", "levels", "values", "l")
      unknown <- setdiff(names(spec), known_fields)
      if (length(unknown) || (length(spec) && is.null(names(spec)))) {
        stop(sprintf(
          "Factor '%s' has unrecognised field(s): %s. Known fields: %s.",
          nm,
          if (length(unknown)) paste(unknown, collapse = ", ") else "(unnamed)",
          paste(known_fields, collapse = ", ")
        ), call. = FALSE)
      }
      # `spec[["..."]]` throughout: `$` partially matches, so `spec$l` would
      # return the `levels` element of any spec that carries one.
      spec[["type"]] <- .dkge_factor_type(spec[["type"]], nm)
      if (identical(spec[["type"]], "continuous")) {
        if (is.null(spec[["values"]])) {
          stop(sprintf("Continuous factor '%s' must provide `values`.", nm),
               call. = FALSE)
        }
        spec[["L"]] <- length(spec[["values"]])
        if (is.null(spec[["levels"]])) {
          spec[["levels"]] <- as.character(spec[["values"]])
        }
      } else {
        if (is.null(spec[["L"]]) && !is.null(spec[["levels"]])) {
          spec[["L"]] <- length(spec[["levels"]])
        }
        if (is.null(spec[["L"]])) {
          stop(sprintf(
            "Factor '%s' must provide `L` (number of levels), `levels`, or `values` for a continuous factor.",
            nm
          ), call. = FALSE)
        }
        spec[["L"]] <- .dkge_factor_level_count(spec[["L"]], nm)
        if (is.null(spec[["levels"]])) {
          spec[["levels"]] <- paste0(nm, seq_len(spec[["L"]]))
        }
      }
      spec[["L"]] <- .dkge_factor_level_count(spec[["L"]], nm)
      if (length(spec[["levels"]]) != spec[["L"]]) {
        stop(sprintf("Factor '%s': `levels` has %d entries but `L` is %d.",
                     nm, length(spec[["levels"]]), spec[["L"]]), call. = FALSE)
      }
      spec[["levels"]] <- .dkge_factor_levels(spec[["levels"]], nm)
      if (!is.null(spec[["l"]])) {
        l <- spec[["l"]]
        if (!is.numeric(l) || length(l) != 1L || !is.finite(l) || l <= 0) {
          stop(sprintf(
            "Factor '%s': length-scale `l` must be a positive finite scalar.",
            nm
          ), call. = FALSE)
        }
      }
      return(spec)
    }
    if (is.numeric(x) && length(x) == 1L) {
      L <- .dkge_factor_level_count(x, nm)
      return(list(L = L, type = "nominal", levels = paste0(nm, seq_len(L))))
    }
    lev <- as.character(x)
    if (!length(lev)) {
      stop(sprintf("Factor '%s' must have at least one level.", nm), call. = FALSE)
    }
    list(L = length(lev), type = "nominal",
         levels = .dkge_factor_levels(lev, nm))
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

  grid_labels <- .dkge_effect_grid_cell_labels(specs, sep = sep)

  out <- list(
    factors = specs,
    factor_names = factor_names,
    default_terms = .dkge_default_kernel_terms(factor_names),
    scope = scopes,
    block_factors = block_factors,
    cells = grid_labels$cells,
    cell_labels = grid_labels$labels,
    sep = sep
  )
  class(out) <- "dkge_effect_grid"
  out
}

#' Derive unique default kernel terms from an ordered factor set
#'
#' Main effects and the full interaction are distinct only when at least two
#' factors are present. Keeping this derivation next to dkge_effect_grid() makes
#' the grid and design_kernel() constructors share one ordering contract.
#'
#' @noRd
.dkge_default_kernel_terms <- function(factor_names) {
  main_effects <- as.list(as.character(factor_names))
  if (length(factor_names) > 1L) {
    c(main_effects, list(as.character(factor_names)))
  } else {
    main_effects
  }
}

#' Expand factor levels into cell labels
#'
#' Shared by [dkge_effect_grid()] and [design_kernel()] so that both produce the
#' same cell ordering (first factor varying slowest).
#'
#' @noRd
.dkge_effect_grid_cell_labels <- function(factors, sep = ":") {
  levels <- lapply(names(factors), function(nm) {
    f <- factors[[nm]]
    if (!is.null(f[["levels"]])) {
      as.character(f[["levels"]])
    } else {
      paste0(nm, seq_len(f[["L"]]))
    }
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
  unknown <- setdiff(as.character(term), names(factor_scope))
  if (length(unknown)) {
    stop("Term references unknown factor(s): ",
         paste(unknown, collapse = ", "), call. = FALSE)
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
