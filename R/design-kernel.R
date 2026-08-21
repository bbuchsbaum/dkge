
# design-kernel.R
# Flexible constructors for design similarity kernels.

#' Build a flexible design-similarity kernel
#'
#' Constructs a PSD kernel that captures factorial similarity across design
#' effects, optionally mapping from cell space to effect space using contrasts.
#'
#' @param factors Named list of factor specifications. Each factor is described
#'   by a list containing:
#'   \describe{
#'     \item{type}{"nominal" | "ordinal" | "circular" | "continuous" (default "nominal").}
#'     \item{L}{Number of levels (for discrete types).}
#'     \item{values}{Numeric coordinates for continuous factors.}
#'     \item{l}{Optional length-scale for ordinal/circular/continuous factors.}
#'   }
#' @param terms List of character vectors describing which factors appear in each
#'   kernel term (e.g. list("A","B", c("A","B"))). Defaults to all main
#'   effects plus the full interaction when there is more than one factor. For
#'   a one-factor design, the main effect and full interaction coincide and are
#'   included once.
#' @param rho Named numeric weights per term (names like "A", "A:B"). Defaults
#'   to 1 for each term if omitted. Must be non-negative.
#' @param include_intercept Logical; if TRUE adds a small identity ridge
#'   (controlled by `rho0`) to keep the kernel full rank (default TRUE).
#' @param rho0 Non-negative scalar ridge weight added when `include_intercept` is
#'   TRUE (default 1e-8).
#' @param basis Either "cell" (kernel over all design cells) or "effect" (kernel
#'   over regressors/effects). Default "cell".
#' @param contrasts Optional named list of per-factor contrast matrices used when
#'   `basis="effect"`. Defaults to sum-to-zero contrasts for discrete factors and
#'   a single column of ones for continuous factors.
#' @param block_structure Optional ordering of effect blocks (names matching
#'   terms). If NULL, uses the order of `terms`.
#' @param block_factors Optional factor names forced to use their factor kernel
#'   in every cell-space term. This gives block-diagonal replication across
#'   between-subject factors such as group instead of accidental all-ones
#'   coupling through terms that omit the block factor.
#' @param normalize One of "unit_trace", "none", "unit_fro", "max_diag".
#'   Controls how the kernel is scaled after construction (default "unit_trace").
#' @param jitter Small diagonal jitter added to `K_cell` for numerical stability
#'   (default 1e-8).
#' @return A list with elements `K` (kernel in requested basis), `K_cell` (always
#'   returned), and `info` containing metadata such as factor/term names, mapping
#'   matrix, block indices, explicit cell/effect coordinate spaces, and the
#'   labels on each coordinate axis.
#' @examples
#' # Simple 2x3 factorial design
#' kern <- design_kernel(
#'   factors = list(A = list(L = 2), B = list(L = 3)),
#'   basis = "effect"
#' )
#' dim(kern$K)  # 5x5 effect-space kernel
#'
#' # Ordinal factor with RBF smoothing
#' kern_ord <- design_kernel(
#'   factors = list(time = list(L = 5, type = "ordinal", l = 1.5)),
#'   basis = "effect"
#' )
#' @export

design_kernel <- function(factors,
                          terms = NULL,
                          rho = NULL,
                          include_intercept = TRUE,
                          rho0 = 1e-8,
                          basis = c("cell", "effect"),
                          contrasts = NULL,
                          block_structure = NULL,
                          block_factors = NULL,
                          normalize = c("unit_trace", "none", "unit_fro", "max_diag"),
                          jitter = 1e-8) {
  basis <- match.arg(basis)
  normalize <- match.arg(normalize)

  effect_grid <- NULL
  if (inherits(factors, "dkge_effect_grid")) {
    effect_grid <- factors
    factors <- effect_grid$factors
    if (is.null(block_factors)) {
      block_factors <- effect_grid$block_factors
    }
  }

  stopifnot(is.list(factors), length(factors) > 0)
  fact_names <- names(factors)
  if (is.null(fact_names) || any(!nzchar(fact_names)) || any(duplicated(fact_names))) {
    stop("`factors` must be a named list with unique names.")
  }
  if (!is.null(block_factors)) {
    block_factors <- as.character(block_factors)
    missing_blocks <- setdiff(block_factors, fact_names)
    if (length(missing_blocks)) {
      stop("Unknown block factor(s): ", paste(missing_blocks, collapse = ", "))
    }
  }

  # Normalise factor specifications
  for (nm in fact_names) {
    f <- factors[[nm]]
    if (!is.list(f)) {
      stop("Factor '", nm, "' must be a list specification (see ?design_kernel).")
    }
    known_fields <- c("type", "L", "levels", "values", "l")
    unknown <- setdiff(names(f), known_fields)
    if (length(unknown) || (length(f) && is.null(names(f)))) {
      stop(sprintf(
        "Factor '%s' has unrecognised field(s): %s. Known fields: %s.",
        nm,
        if (length(unknown)) paste(unknown, collapse = ", ") else "(unnamed)",
        paste(known_fields, collapse = ", ")
      ), call. = FALSE)
    }
    # `f[["..."]]` throughout: `$` partially matches, so `f$l` would return
    # `levels` and `f$val` would return `values` for specs that carry them.
    f[["type"]] <- .dkge_factor_type(f[["type"]], nm)
    if (identical(f[["type"]], "continuous")) {
      if (is.null(f[["values"]])) stop("Continuous factor '", nm, "' must provide `values`.")
      f[["L"]] <- length(f[["values"]])
    } else {
      if (is.null(f[["L"]]) && !is.null(f[["levels"]])) f[["L"]] <- length(f[["levels"]])
      if (is.null(f[["L"]])) stop("Discrete factor '", nm, "' must provide `L` (number of levels).")
    }
    f[["L"]] <- .dkge_factor_level_count(f[["L"]], nm)
    if (!is.null(f[["levels"]])) {
      if (length(f[["levels"]]) != f[["L"]]) {
        stop(sprintf(
          "Factor '%s': `levels` has %d entries but `L` is %d.",
          nm, length(f[["levels"]]), f[["L"]]
        ))
      }
      f[["levels"]] <- .dkge_factor_levels(f[["levels"]], nm)
    }
    factors[[nm]] <- f
  }

  Ls <- vapply(factors, function(f) f[["L"]], integer(1))
  types <- vapply(factors, function(f) f[["type"]], character(1))
  Qcell <- prod(Ls)
  grid_labels <- if (is.null(effect_grid)) {
    .dkge_effect_grid_cell_labels(factors)
  } else {
    list(cells = effect_grid$cells, labels = effect_grid$cell_labels)
  }
  factor_scope <- if (is.null(effect_grid)) {
    setNames(rep("within", length(fact_names)), fact_names)
  } else {
    effect_grid$scope
  }

  .rbf1d <- function(x, l) {
    x <- as.numeric(x)
    D2 <- outer(x, x, function(a, b) (a - b)^2)
    exp(- D2 / (2 * l * l))
  }

  k_factor <- function(f) {
    L <- f[["L"]]
    type <- f[["type"]]
    # `f[["l"]]`, not `f$l`: `$` partially matches, so a spec carrying `levels`
    # (as every dkge_effect_grid() spec does) would hand the level labels back
    # as the length-scale.
    l <- f[["l"]] %||% 1.0
    if (type != "nominal" &&
        (!is.numeric(l) || length(l) != 1L || !is.finite(l) || l <= 0)) {
      stop("Factor length-scale `l` must be a positive finite scalar.")
    }
    if (type == "nominal") {
      diag(L)
    } else if (type == "ordinal") {
      .rbf1d(seq_len(L), l)
    } else if (type == "circular") {
      idx <- seq_len(L)
      D2 <- outer(idx, idx, function(i, j) {
        d <- abs(i - j)
        dd <- pmin(d, L - d)  # pmin is vectorized, min is not
        dd * dd
      })
      exp(- D2 / (2 * l * l))
    } else if (type == "continuous") {
      values <- as.numeric(f[["values"]])
      l <- f[["l"]] %||% (stats::IQR(values) / 1.349 + 1e-8)
      .rbf1d(values, l)
    } else {
      stop("Unsupported factor type '", type, "'.")
    }
  }

  K_fac <- lapply(factors, k_factor)
  J_fac <- lapply(Ls, function(L) matrix(1, L, L))

  terms_supplied <- !is.null(terms)
  if (!terms_supplied) {
    terms <- if (!is.null(effect_grid$default_terms)) {
      effect_grid$default_terms
    } else {
      .dkge_default_kernel_terms(fact_names)
    }
  }
  if (!is.list(terms)) terms <- as.list(terms)
  terms <- lapply(seq_along(terms), function(k) {
    S <- as.character(terms[[k]])
    unknown <- setdiff(S, fact_names)
    if (length(unknown)) {
      stop(sprintf(
        "Term %d references unknown factor(s): %s. Known factors: %s.",
        k, paste(unknown, collapse = ", "), paste(fact_names, collapse = ", ")
      ))
    }
    if (anyDuplicated(S)) {
      stop(sprintf("Term %d repeats a factor: %s.", k,
                   paste(unique(S[duplicated(S)]), collapse = ", ")))
    }
    S
  })
  term_name <- function(S) paste(S, collapse = ":")
  tnames <- vapply(terms, term_name, character(1))
  if (terms_supplied && anyDuplicated(tnames)) {
    stop("Duplicate kernel terms: ",
         paste(unique(tnames[duplicated(tnames)]), collapse = ", "))
  }

  # `rho` is resolved positionally into one weight per term. Indexing by name
  # (`rho[[tnames[k]]]`) failed with "subscript out of bounds" for a partially
  # named `rho`, and could not express duplicated term names at all.
  rho_vec <- rep(1, length(tnames))
  if (!is.null(rho)) {
    if (!is.numeric(rho) || !length(rho)) {
      stop("`rho` must be a non-empty numeric vector of term weights.")
    }
    if (anyNA(rho) || any(!is.finite(rho))) {
      stop("`rho` must be finite.")
    }
    if (any(rho < 0)) stop("`rho` must be non-negative.")
    rho_names <- names(rho)
    if (is.null(rho_names)) {
      if (length(rho) != length(tnames)) {
        stop(sprintf("Unnamed `rho` must have one entry per term (got %d, expected %d).",
                     length(rho), length(tnames)))
      }
      rho_vec <- as.numeric(rho)
    } else {
      if (any(!nzchar(rho_names))) {
        stop("`rho` must be either fully named or fully unnamed.")
      }
      if (anyDuplicated(rho_names)) {
        stop("`rho` repeats term name(s): ",
             paste(unique(rho_names[duplicated(rho_names)]), collapse = ", "))
      }
      unknown_rho <- setdiff(rho_names, tnames)
      if (length(unknown_rho)) {
        stop("`rho` names must be term names. Unknown: ",
             paste(unknown_rho, collapse = ", "),
             ". Known terms: ", paste(tnames, collapse = ", "))
      }
      # Terms `rho` does not mention keep the default weight of 1.
      idx <- match(tnames, rho_names)
      rho_vec[!is.na(idx)] <- as.numeric(rho[idx[!is.na(idx)]])
    }
  }
  names(rho_vec) <- tnames

  .kron_all <- function(mats) Reduce(kronecker, mats)

  per_term_kron <- function(S) {
    mats <- Map(function(nm, i) if (nm %in% S || nm %in% block_factors) K_fac[[nm]] else J_fac[[i]],
                fact_names, seq_along(fact_names))
    .kron_all(mats)
  }

  K_cell <- matrix(0, Qcell, Qcell)
  for (k in seq_along(terms)) {
    K_cell <- K_cell + rho_vec[[k]] * per_term_kron(terms[[k]])
  }
  if (include_intercept && rho0 > 0) K_cell <- K_cell + rho0 * diag(Qcell)
  if (jitter > 0)                     K_cell <- K_cell + jitter * diag(Qcell)
  dimnames(K_cell) <- list(grid_labels$labels, grid_labels$labels)

  if (normalize == "unit_trace") {
    tr <- sum(diag(K_cell)); if (tr > 0) K_cell <- K_cell / tr
  } else if (normalize == "unit_fro") {
    fn <- sqrt(sum(K_cell * K_cell)); if (fn > 0) K_cell <- K_cell / fn
  } else if (normalize == "max_diag") {
    md <- max(diag(K_cell)); if (md > 0) K_cell <- K_cell / md
  }

  info <- list(levels = Ls,
               factor_names = fact_names,
               term_names = tnames,
               factor_scope = factor_scope,
               term_scope = setNames(vapply(terms, .dkge_term_scope, character(1),
                                             factor_scope = factor_scope), tnames),
               block_factors = block_factors,
               cells = grid_labels$cells,
               cell_labels = grid_labels$labels,
               effect_labels = NULL,
               kernel_labels = NULL,
               basis = basis,
               coordinate_space = list(
                 kernel = basis,
                 kernel_labels = basis,
                 cell_labels = "cell",
                 cells = "cell",
                 map_rows = "cell",
                 map_columns = "effect",
                 blocks = basis
               ),
               map = NULL,
               blocks = NULL,
               dims = list(Qcell = Qcell))

  if (basis == "cell") {
    info$kernel_labels <- grid_labels$labels
    info$blocks <- list(cells = seq_len(Qcell))
    return(list(K = K_cell, K_cell = K_cell, info = info))
  }

  if (is.null(contrasts)) {
    contrasts <- lapply(factors, function(f) {
      if (f[["type"]] %in% c("nominal", "ordinal", "circular")) {
        cm <- contr.sum(f[["L"]])
        storage.mode(cm) <- "double"
        cm
      } else {
        matrix(1, f[["L"]], 1)
      }
    })
  } else {
    if (!all(names(contrasts) %in% fact_names)) {
      stop("All `contrasts` must be named and match factor names.")
    }
  }

  one_column <- function(L) matrix(1, L, 1)
  term_map <- function(S) {
    mats <- lapply(fact_names, function(nm) {
      if (nm %in% S) contrasts[[nm]] else one_column(factors[[nm]][["L"]])
    })
    T <- .kron_all(mats)
    attr(T, "out_dim") <- ncol(T)
    T
  }

  if (is.null(block_structure)) block_structure <- tnames
  T_blocks <- setNames(lapply(terms, term_map), tnames)
  missing <- setdiff(block_structure, names(T_blocks))
  if (length(missing)) stop("Blocks requested but not present in terms: ",
                             paste(missing, collapse = ", "))

  T_list <- T_blocks[block_structure]
  out_cols <- vapply(T_list, function(T) attr(T, "out_dim"), integer(1))
  T <- do.call(cbind, T_list)
  rownames(T) <- grid_labels$labels

  effect_names <- unlist(Map(function(block, n) {
    if (n == 1L) block else paste0(block, seq_len(n))
  }, block_structure, out_cols), use.names = FALSE)
  colnames(T) <- effect_names

  block_idx <- split(seq_len(sum(out_cols)),
                     rep(block_structure, times = out_cols))

  info$map    <- T
  info$blocks <- block_idx
  info$dims$q <- ncol(T)
  info$effect_labels <- effect_names
  info$kernel_labels <- effect_names

  K <- crossprod(T, K_cell %*% T)
  dimnames(K) <- list(effect_names, effect_names)
  list(K = K, K_cell = K_cell, info = info)
}

#' Sum-to-zero contrasts for a set of factors
#' @param Ls Numeric vector of factor levels (named or unnamed).
#' @return Named list of Lx(L-1) contrast matrices.
#' @export
sum_contrasts <- function(Ls) {
  lapply(Ls, function(L) { cm <- contr.sum(L); storage.mode(cm) <- "double"; cm })
}

#' Helmert contrasts for a set of factors
#' @param Ls Numeric vector of factor levels (named or unnamed).
#' @return Named list of orthonormal Helmert contrast matrices.
#' @export
helmert_contrasts <- function(Ls) {
  lapply(Ls, function(L) {
    cm <- contr.helmert(L)
    Q  <- qr.Q(qr(cm))
    Q[, seq_len(ncol(cm)), drop = FALSE]
  })
}

#' Robust kernel roots
#' @param K Positive semi-definite kernel matrix.
#' @param jitter Small diagonal jitter added before inversion.
#' @return List with `Khalf`, `Kihalf`, eigenvalues, eigenvectors, and basic diagnostics.
#' @export
kernel_roots <- function(K, jitter = 1e-10) {
  if (!isTRUE(isSymmetric(K, tol = 1e-8))) {
    warning("Kernel matrix not symmetric; applying symmetrization via (K + t(K))/2.")
  }
  Ks <- (K + t(K)) / 2
  ee <- eigen(Ks, symmetric = TRUE)
  vals <- ee$values
  V <- ee$vectors

  if (is.null(jitter)) {
    pinned <- vals <= 0
    tiny <- 1e-10
    vals[pinned] <- tiny
  } else {
    if (!is.numeric(jitter) || length(jitter) != 1 || jitter < 0) {
      stop("`jitter` must be NULL or a non-negative numeric scalar")
    }
    pinned <- vals < jitter
    vals <- pmax(vals, jitter)
  }

  n_clamped <- sum(pinned)
  n_total <- length(vals)
  if (n_clamped > 0 && n_clamped / max(n_total, 1) > 0.1) {
    warning(sprintf("%.0f%% of eigenvalues were clamped during kernel root computation.",
                    100 * n_clamped / n_total))
  }

  sqrt_vals <- sqrt(vals)
  inv_sqrt_vals <- rep(0, length(vals))
  pos_idx <- vals > 0
  inv_sqrt_vals[pos_idx] <- 1 / sqrt_vals[pos_idx]

  Khalf  <- V %*% diag(sqrt_vals, length(sqrt_vals)) %*% t(V)
  Kihalf <- V %*% diag(inv_sqrt_vals, length(inv_sqrt_vals)) %*% t(V)
  dimnames(Khalf) <- dimnames(Ks)
  dimnames(Kihalf) <- dimnames(Ks)

  rank_est <- sum(vals > .Machine$double.eps)

  list(Khalf = Khalf,
       Kihalf = Kihalf,
       evals = vals,
       evecs = V,
       rank = rank_est,
       n_clamped = n_clamped)
}

#' Kernel alignment score
#' @param A,B Matrices of identical dimensions.
#' @return Cosine similarity between the flattened matrices.
#' @export
kernel_alignment <- function(A, B) {
  stopifnot(all(dim(A) == dim(B)))
  num <- sum(A * B)
  den <- sqrt(sum(A * A) * sum(B * B) + 1e-24)
  if (den == 0) 0 else as.numeric(num / den)
}

#' @rdname design_kernel
#' @export
dkge_design_kernel <- design_kernel

#' @rdname kernel_roots
#' @export
dkge_kernel_roots <- kernel_roots

#' @rdname kernel_alignment
#' @export
dkge_kernel_alignment <- kernel_alignment
