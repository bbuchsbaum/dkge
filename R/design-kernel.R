
# design-kernel.R
# Flexible constructors for design similarity kernels.

"%||%" <- function(a, b) if (is.null(a)) b else a

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
#'   effects plus the full interaction.
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
#' @param normalize One of "unit_trace", "none", "unit_fro", "max_diag".
#'   Controls how the kernel is scaled after construction (default "unit_trace").
#' @param jitter Small diagonal jitter added to `K_cell` for numerical stability
#'   (default 1e-8).
#' @return A list with elements `K` (kernel in requested basis), `K_cell` (always
#'   returned), and `info` containing metadata such as factor/term names, mapping
#'   matrix, and block indices.
#' @export

design_kernel <- function(factors,
                          terms = NULL,
                          rho = NULL,
                          include_intercept = TRUE,
                          rho0 = 1e-8,
                          basis = c("cell", "effect"),
                          contrasts = NULL,
                          block_structure = NULL,
                          normalize = c("unit_trace", "none", "unit_fro", "max_diag"),
                          jitter = 1e-8) {
  basis <- match.arg(basis)
  normalize <- match.arg(normalize)

  stopifnot(is.list(factors), length(factors) > 0)
  fact_names <- names(factors)
  if (is.null(fact_names) || any(!nzchar(fact_names)) || any(duplicated(fact_names))) {
    stop("`factors` must be a named list with unique names.")
  }

  # Normalise factor specifications
  for (nm in fact_names) {
    f <- factors[[nm]]
    f$type <- tolower(f$type %||% "nominal")
    if (f$type == "continuous") {
      if (is.null(f$values)) stop("Continuous factor '", nm, "' must provide `values`.")
      f$L <- length(f$values)
    } else {
      if (is.null(f$L)) stop("Discrete factor '", nm, "' must provide `L` (number of levels).")
    }
    f$L <- as.integer(f$L)
    factors[[nm]] <- f
  }

  Ls <- vapply(factors, function(f) f$L, integer(1))
  types <- vapply(factors, function(f) f$type, character(1))
  Qcell <- prod(Ls)

  .rbf1d <- function(x, l) {
    x <- as.numeric(x)
    D2 <- outer(x, x, function(a, b) (a - b)^2)
    exp(- D2 / (2 * l * l))
  }

  k_factor <- function(f) {
    L <- f$L
    type <- f$type
    l <- f$l %||% 1.0
    if (type == "nominal") {
      diag(L)
    } else if (type == "ordinal") {
      .rbf1d(seq_len(L), l)
    } else if (type == "circular") {
      idx <- seq_len(L)
      D2 <- outer(idx, idx, function(i, j) {
        d <- abs(i - j)
        dd <- min(d, L - d)
        dd * dd
      })
      exp(- D2 / (2 * l * l))
    } else if (type == "continuous") {
      values <- as.numeric(f$values)
      l <- f$l %||% (stats::IQR(values) / 1.349 + 1e-8)
      .rbf1d(values, l)
    } else {
      stop("Unsupported factor type '", type, "'.")
    }
  }

  K_fac <- lapply(factors, k_factor)
  J_fac <- lapply(Ls, function(L) matrix(1, L, L))

  if (is.null(terms)) {
    terms <- c(as.list(fact_names), list(fact_names))
  }
  term_name <- function(S) paste(S, collapse = ":")
  tnames <- vapply(terms, term_name, character(1))

  if (is.null(rho)) {
    rho <- setNames(rep(1, length(terms)), tnames)
  } else {
    if (is.null(names(rho))) names(rho) <- tnames
    if (any(rho < 0)) stop("`rho` must be non-negative.")
  }

  .kron_all <- function(mats) Reduce(kronecker, mats)

  per_term_kron <- function(S) {
    mats <- Map(function(nm, i) if (nm %in% S) K_fac[[nm]] else J_fac[[i]],
                fact_names, seq_along(fact_names))
    .kron_all(mats)
  }

  K_cell <- matrix(0, Qcell, Qcell)
  for (k in seq_along(terms)) {
    K_cell <- K_cell + (rho[[tnames[k]]] %||% 1) * per_term_kron(terms[[k]])
  }
  if (include_intercept && rho0 > 0) K_cell <- K_cell + rho0 * diag(Qcell)
  if (jitter > 0)                     K_cell <- K_cell + jitter * diag(Qcell)

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
               basis = basis,
               map = NULL,
               blocks = NULL,
               dims = list(Qcell = Qcell))

  if (basis == "cell") {
    return(list(K = K_cell, K_cell = K_cell, info = info))
  }

  if (is.null(contrasts)) {
    contrasts <- lapply(factors, function(f) {
      if (f$type %in% c("nominal", "ordinal", "circular")) {
        cm <- contr.sum(f$L)
        storage.mode(cm) <- "double"
        cm
      } else {
        matrix(1, f$L, 1)
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
      if (nm %in% S) contrasts[[nm]] else one_column(factors[[nm]]$L)
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

  block_idx <- split(seq_len(sum(out_cols)),
                     rep(block_structure, times = out_cols))

  info$map    <- T
  info$blocks <- block_idx
  info$dims$q <- ncol(T)

  K <- crossprod(T, K_cell %*% T)
  list(K = K, K_cell = K_cell, info = info)
}

#' Sum-to-zero contrasts for a set of factors
#' @param Ls Numeric vector of factor levels (named or unnamed).
#' @return Named list of L×(L-1) contrast matrices.
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
#' @return List with `Khalf`, `Kihalf`, eigenvalues, and eigenvectors.
#' @export
kernel_roots <- function(K, jitter = 1e-10) {
  Ks <- (K + t(K)) / 2
  ee <- eigen(Ks, symmetric = TRUE)
  vals <- pmax(ee$values, jitter)
  V <- ee$vectors
  Khalf  <- V %*% diag(sqrt(vals)) %*% t(V)
  Kihalf <- V %*% diag(1 / sqrt(vals)) %*% t(V)
  list(Khalf = Khalf, Kihalf = Kihalf, evals = vals, evecs = V)
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
