#' @import ggplot2
NULL

# -----------------------------------------------------------------------------
# Shared aesthetics -----------------------------------------------------------
# -----------------------------------------------------------------------------

#' DKGE minimal theme for ggplot2 outputs
#'
#' Produces a light-weight theme used by the DKGE plotting helpers. Adjust
#' `base_size` / `base_family` to customize font size or typeface.
#'
#' @param base_size Base font size.
#' @param base_family Base font family.
#' @return A ggplot2 theme object.
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   ggplot2::ggplot(mtcars, ggplot2::aes(mpg, wt)) +
#'     ggplot2::geom_point() +
#'     theme_dkge()
#' }
#' @export
theme_dkge <- function(base_size = 12, base_family = "") {
 ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
 ggplot2::theme(
 panel.grid.minor = ggplot2::element_blank(),
 panel.grid.major = ggplot2::element_line(linewidth = 0.25, colour = "grey85"),
 axis.title = ggplot2::element_text(face = "bold"),
 plot.title = ggplot2::element_text(face = "bold", hjust = 0, margin = ggplot2::margin(b = 6)),
 plot.subtitle = ggplot2::element_text(hjust = 0, margin = ggplot2::margin(b = 10)),
 legend.position = "right",
 strip.text = ggplot2::element_text(face = "bold")
 )
}

.dkge_scale_fill_diverging <- function() {
 ggplot2::scale_fill_gradient2(low = "#3B4CC0", mid = "white", high = "#B40426", midpoint = 0)
}
.dkge_scale_colour_diverging <- function() {
 ggplot2::scale_colour_gradient2(low = "#3B4CC0", mid = "white", high = "#B40426", midpoint = 0)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

.dkge_effect_indexer <- function(fit) {
 q <- nrow(fit$U)
 eff_names <- fit$effects %||% paste0("effect", seq_len(q))
 blocks <- fit$kernel_info$blocks %||% list(seq_len(q))
 if (!is.list(blocks) || length(unlist(blocks)) != q) {
 blocks <- list(seq_len(q))
 }
 list(names = eff_names, blocks = blocks)
}

.dkge_component_indexer <- function(fit, comps = NULL) {
 r <- ncol(fit$U)
 if (is.null(comps)) comps <- seq_len(min(r, 6L))
 comps <- suppressWarnings(as.integer(comps))
 if (!length(comps) || anyNA(comps)) {
 stop("`comps` must be a non-empty vector of component indices.", call. = FALSE)
 }
 bad <- unique(comps[comps < 1L | comps > r])
 if (length(bad)) {
 stop(sprintf("Component index out of range (fit has rank %d): %s",
 r, paste(bad, collapse = ", ")), call. = FALSE)
 }
 comps
}

.dkge_component_labels <- function(comps) {
 paste0("LV", comps)
}

.dkge_energy_by_subject_component <- function(fit) {
 r <- ncol(fit$U)
 S <- length(fit$Btil)
 subjects <- fit$subject_ids %||% paste0("sub", seq_len(S))
 energy <- matrix(NA_real_, nrow = S, ncol = r)
 for (s in seq_len(S)) {
 B_t <- fit$Btil[[s]]
 w_s <- .dkge_fit_subject_voxel_weights(fit, s, B_t, subject_id = subjects[[s]])
 B_t <- .dkge_apply_voxel_weights(B_t, w_s)
 A_s <- t(B_t) %*% fit$K %*% fit$U
 energy[s, ] <- colSums(A_s * A_s)
 }
 rownames(energy) <- subjects
 colnames(energy) <- paste0("comp", seq_len(r))
 energy
}

# Euclidean orthonormal basis of the image of `U` under K^{1/2}, truncated to
# the numerical rank. Working in the K^{1/2} coordinates turns K-metric
# principal angles into ordinary ones, and the truncation is what keeps a
# rank-deficient input honest: `qr.Q()` would pad the deficient directions with
# an arbitrary complement, inventing extra "shared" directions between the two
# subspaces (span{e1,e2} against a duplicated-column basis of span{e3} then
# reports angles 0 and 90 instead of a single 90).
.dkge_k_span_basis <- function(U, Khalf, label) {
 W <- Khalf %*% U
 sv <- svd(W)
 d <- sv$d
 tol <- max(dim(W)) * .Machine$double.eps * max(d, 0)
 rank <- sum(d > tol)
 if (rank == 0L) {
 stop(sprintf("`%s` spans nothing in the K metric (numerical rank 0).", label),
 call. = FALSE)
 }
 sv$u[, seq_len(rank), drop = FALSE]
}

.dkge_principal_angles_K <- function(U1, U2, K, orthonormalize = TRUE) {
 U1 <- as.matrix(U1)
 U2 <- as.matrix(U2)
 K <- as.matrix(K)
 if (nrow(U1) != nrow(K) || nrow(U2) != nrow(K)) {
 stop("`U1` and `U2` must have as many rows as `K`.", call. = FALSE)
 }
 if (orthonormalize) {
 Khalf <- .dkge_kernel_roots(K)$Khalf
 Q1 <- .dkge_k_span_basis(U1, Khalf, "U1")
 Q2 <- .dkge_k_span_basis(U2, Khalf, "U2")
 # min(rank(U1), rank(U2)) angles, not min(ncol(U1), ncol(U2)).
 sv <- svd(crossprod(Q1, Q2), nu = 0, nv = 0)$d
 } else {
 sv <- svd(t(U1) %*% K %*% U2, nu = 0, nv = 0)$d
 }
 # Singular values are non-negative by construction; only the upper clamp is
 # needed to absorb floating-point overshoot past 1.
 acos(pmin(sv, 1))
}

#' Principal angles between two DKGE bases in the K metric
#'
#' Computes subspace principal angles between two bases using the DKGE design
#' kernel as the inner-product metric. This is useful for comparing stratified
#' or resampled DKGE fits after K-Procrustes alignment diagnostics.
#'
#' @details
#' Principal angles are a property of the *subspaces* spanned by `U1` and `U2`,
#' not of the particular bases used to represent them. Both inputs are therefore
#' reduced to an ordinary Euclidean basis of \eqn{K^{1/2} U} via
#' \code{svd(K^{1/2} U)} (the internal helper `.dkge_k_span_basis()`) before
#' the cosines \eqn{\sigma_i(Q_1^\top Q_2)} are formed, so any two bases of
#' the same subspace (for example `U` and `2 * U`) return angles of zero. Set
#' `orthonormalize = FALSE` only when the inputs are already known to satisfy
#' \eqn{U^\top K U = I}; the singular values are otherwise not cosines and the
#' clamp to \eqn{[0, 1]} would silently report perfect alignment.
#'
#' The returned vector has one entry per shared dimension, ordered from the
#' smallest angle (most aligned direction) to the largest. Rank-deficient inputs
#' are truncated to the subspace they actually span, so the length is
#' `min(rank_K(U1), rank_K(U2))`, which can be smaller than
#' `min(ncol(U1), ncol(U2))`; a basis of numerical rank 0 is an error rather
#' than a fabricated subspace. With `orthonormalize = FALSE` the inputs are
#' taken at face value and all `min(ncol(U1), ncol(U2))` singular values are
#' returned.
#'
#' @param U1,U2 Basis matrices with the same row dimension as `K`. Columns need
#'   not be K-orthonormal.
#' @param K Design kernel defining the inner product.
#' @param orthonormalize Logical; K-orthonormalize `U1` and `U2` first
#'   (default `TRUE`).
#' @return Numeric vector of principal angles in radians.
#' @examples
#' K <- diag(c(2, 1, 1))
#' U1 <- cbind(c(1, 0, 0), c(0, 1, 0))
#' # A rescaled basis of the same subspace has zero angles.
#' dkge_principal_angles_K(U1, 3 * U1, K)
#' # A K-orthogonal direction sits at 90 degrees.
#' dkge_principal_angles_K(U1[, 1, drop = FALSE], U1[, 2, drop = FALSE], K)
#' @export
dkge_principal_angles_K <- function(U1, U2, K, orthonormalize = TRUE) {
 .dkge_principal_angles_K(U1, U2, K, orthonormalize = orthonormalize)
}

#' DKGE component saliences in effect space
#'
#' Extracts the design/effect-space saliences \eqn{K U}. These are the raw
#' coordinates most users should inspect first: each row is an effect or design
#' cell, each column is a DKGE latent variable.
#'
#' @details
#' The salience matrix is exactly \eqn{K U[, comps]}. Because `fit$U` is
#' K-orthonormal (\eqn{U^\top K U = I}), each salience column already has unit
#' norm in the metric dual to `K`, i.e.
#' \eqn{s_j^\top K^{-1} s_j = u_j^\top K u_j = 1}. `scale = "unit"` therefore
#' divides column `j` by \eqn{\sqrt{u_j^\top K u_j}} and is a no-op for a
#' well-formed fit; it is retained so that hand-assembled or perturbed bases are
#' put back on a comparable footing. It deliberately does **not** divide by the
#' Euclidean norm, which would distort the K geometry.
#'
#' @param fit Fitted `dkge` object.
#' @param comps Components to include (defaults to first min(rank, 6)).
#' @param scale Optional within-component display scaling. `"raw"` leaves
#'   saliences on their original scale, `"unit"` rescales each component to unit
#'   K-norm of the underlying latent vector (see Details), and `"zscore"`
#'   z-scores each component across effects.
#' @param long Logical; return a tidy long data frame when `TRUE`, otherwise a
#'   numeric effects-by-components matrix.
#' @return A data frame or matrix of component saliences.
#' @examples
#' toy <- dkge_sim_toy(
#'   factors = list(A = list(L = 2), B = list(L = 3)),
#'   active_terms = c("A", "B"), S = 3, P = 15, snr = 5
#' )
#' fit <- dkge(toy$B_list, toy$X_list, K = toy$K, rank = 2)
#' head(dkge_component_saliences(fit, comps = 1:2))
#' dkge_component_saliences(fit, comps = 1, long = FALSE)
#' @export
dkge_component_saliences <- function(fit,
                                     comps = NULL,
                                     scale = c("raw", "unit", "zscore"),
                                     long = TRUE) {
 stopifnot(inherits(fit, "dkge"))
 scale <- match.arg(scale)
 comps <- .dkge_component_indexer(fit, comps)
 idx <- .dkge_effect_indexer(fit)

 U_sub <- fit$U[, comps, drop = FALSE]
 sal <- fit$K %*% U_sub
 rownames(sal) <- idx$names
 colnames(sal) <- .dkge_component_labels(comps)

 if (scale == "unit") {
 norms <- sqrt(pmax(colSums(U_sub * sal), 0))
 norms[!is.finite(norms) | norms < 1e-12] <- 1
 sal <- sweep(sal, 2L, norms, "/")
 } else if (scale == "zscore") {
 sal <- apply(sal, 2L, function(x) {
 s <- stats::sd(x)
 if (!is.finite(s) || s < 1e-12) return(rep(0, length(x)))
 (x - mean(x)) / s
 })
 sal <- as.matrix(sal)
 rownames(sal) <- idx$names
 colnames(sal) <- .dkge_component_labels(comps)
 }

 if (!long) return(sal)

 data.frame(
 effect = rep(idx$names, times = length(comps)),
 component = factor(rep(.dkge_component_labels(comps), each = length(idx$names)),
 levels = .dkge_component_labels(comps)),
 component_id = rep(comps, each = length(idx$names)),
 salience = as.vector(sal),
 stringsAsFactors = FALSE
 )
}

.dkge_identity_design_basis <- function(fit, source = "identity",
                                        normalize = "none") {
 idx <- .dkge_effect_indexer(fit)
 q <- length(idx$names)
 C <- diag(q)
 dimnames(C) <- list(idx$names, idx$names)
 # The identity fallback is a *default* basis like any other, so it honours
 # `normalize`. Skipping it left the documented `unit_K` default silently
 # inapplicable to every fit without cell metadata (an effect-basis kernel
 # gives identity columns K-norms of sqrt(diag(K)), not 1).
 C <- .dkge_normalize_basis_columns(C, fit$K, normalize)
 attr(C, "term") <- idx$names
 attr(C, "source") <- source
 C
}

# Level order for factor `nm`. `design_kernel()` stores cells built by
# `expand.grid()` over the declared level labels, in which each column's
# first-appearance order *is* the declared order; it stores no per-factor level
# labels of its own (`info$levels` holds level *counts*). Reading the order off
# the cell table is therefore both the declared order and the only source that
# actually exists.
.dkge_kernel_factor_levels <- function(values) {
 unique(as.character(values))
}

# Reorder cell-basis `kernel_info$cells` so rows follow `fit$effects`, matching
# `cell_labels` by name. Effect-basis cell metadata is deliberately not treated
# as fit coordinates, even when Qcell == q. Missing factorial metadata returns
# NULL; a declared cell mapping that is incomplete or ambiguous fails closed.
.dkge_match_kernel_cells <- function(fit, info = NULL) {
  info <- info %||% fit$kernel_info %||% list()
  labels <- info$cell_labels
  cells <- info$cells
  effects <- fit$effects
  if (is.null(effects)) {
    return(NULL)
  }
  kernel_labels <- rownames(fit$K) %||% colnames(fit$K) %||% effects
  space <- .dkge_kernel_info_space(info, kernel_labels)
  if (identical(space, "effect")) return(NULL)
  if (is.null(cells) && is.null(labels)) return(NULL)
  if (is.null(cells) || is.null(labels)) {
    stop(
      "Cell-basis kernel metadata requires both `cell_labels` and `cells`.",
      call. = FALSE
    )
  }
  if (is.null(space)) {
    stop(
      paste0(
        "kernel_info cell coordinates are ambiguous; declare `basis = 'cell'` ",
        "or `coordinate_space$kernel = 'cell'`."
      ),
      call. = FALSE
    )
  }
  cells <- as.data.frame(cells, stringsAsFactors = FALSE)
  labels <- as.character(labels)
  effects <- as.character(effects)
  if (!nrow(cells) || nrow(cells) != length(labels)) {
    stop("kernel_info$cells must have one row per cell_labels entry.",
         call. = FALSE)
  }
  if (anyDuplicated(labels) || anyDuplicated(effects)) {
    stop("kernel_info$cell_labels and fit$effects must be unique for name matching.",
         call. = FALSE)
  }
  idx <- match(effects, labels)
  if (length(labels) != length(effects) || anyNA(idx) ||
      !setequal(labels, effects)) {
    stop(
      "kernel_info$cell_labels cannot be aligned uniquely to fit$effects.",
      call. = FALSE
    )
  }
  cells[idx, , drop = FALSE]
}

.dkge_quote_factor_names <- function(nms) {
  sprintf("`%s`", nms)
}

.dkge_normalize_basis_columns <- function(C, K, normalize) {
 if (identical(normalize, "none")) return(C)
 norms <- if (identical(normalize, "unit_K")) {
 sqrt(pmax(colSums(C * (K %*% C)), 0))
 } else {
 sqrt(colSums(C * C))
 }
 norms[!is.finite(norms) | norms < 1e-12] <- 1
 sweep(C, 2L, norms, "/")
}

.dkge_default_design_basis <- function(fit,
                                       include_intercept = TRUE,
                                       coding = c("sum", "helmert", "poly"),
                                       normalize = c("unit_K", "unit_l2", "none")) {
 coding <- match.arg(coding)
 normalize <- match.arg(normalize)
 idx <- .dkge_effect_indexer(fit)
 q <- length(idx$names)
 info <- fit$kernel_info %||% list()
 factor_names <- info$factor_names
 cells <- .dkge_match_kernel_cells(fit, info)

 if (is.null(cells) || is.null(factor_names) || nrow(cells) != q) {
 return(.dkge_identity_design_basis(fit, normalize = normalize))
 }

 missing_factors <- setdiff(factor_names, names(cells))
 if (length(missing_factors)) {
 return(.dkge_identity_design_basis(fit, normalize = normalize))
 }

 factor_names <- factor_names[factor_names %in% names(cells)]
 for (nm in factor_names) {
 cells[[nm]] <- factor(cells[[nm]],
 levels = .dkge_kernel_factor_levels(cells[[nm]]))
 }
 factor_names <- factor_names[vapply(factor_names, function(nm) {
   nlevels(cells[[nm]]) > 1L
 }, logical(1))]
 if (!length(factor_names)) {
   return(.dkge_identity_design_basis(fit, normalize = normalize))
 }

 contrast_fun <- switch(coding,
 sum = stats::contr.sum,
 helmert = stats::contr.helmert,
 poly = stats::contr.poly)
 contrasts_arg <- lapply(factor_names, function(nm) contrast_fun(nlevels(cells[[nm]])))
 names(contrasts_arg) <- factor_names

 rhs <- paste(.dkge_quote_factor_names(factor_names), collapse = " * ")
 form <- stats::as.formula(paste("~", rhs))
 mm <- stats::model.matrix(form, data = cells, contrasts.arg = contrasts_arg)
 assign <- attr(mm, "assign")
 term_labels <- gsub("`", "", attr(stats::terms(form), "term.labels"),
                     fixed = TRUE)

 colnames(mm)[assign == 0L] <- "grand_mean"
 for (a in unique(assign[assign > 0L])) {
 loc <- which(assign == a)
 if (length(loc) == 1L) colnames(mm)[loc] <- term_labels[[a]]
 }

 # The model-matrix `assign` attribute maps every column to the term that
 # generated it, so terms are read off exactly instead of being recovered by
 # prefix-matching column names (which is ambiguous when one factor name is a
 # prefix of another, e.g. "task" and "task2").
 term <- character(ncol(mm))
 term[assign == 0L] <- "grand_mean"
 has_term <- assign > 0L
 term[has_term] <- term_labels[assign[has_term]]

 if (!include_intercept) {
 keep <- assign != 0L
 mm <- mm[, keep, drop = FALSE]
 term <- term[keep]
 }

 if (!ncol(mm)) {
 return(.dkge_identity_design_basis(fit, normalize = normalize))
 }
 rownames(mm) <- idx$names

 mm <- .dkge_normalize_basis_columns(mm, fit$K, normalize)
 attr(mm, "assign") <- NULL
 attr(mm, "contrasts") <- NULL

 attr(mm, "term") <- as.character(term)
 attr(mm, "source") <- "factorial"
 mm
}

.dkge_validate_design_basis <- function(fit,
                                        basis = NULL,
                                        include_intercept = TRUE,
                                        coding = c("sum", "helmert", "poly"),
                                        normalize = c("unit_K", "unit_l2", "none")) {
 coding <- match.arg(coding)
 normalize <- match.arg(normalize)
 idx <- .dkge_effect_indexer(fit)
 q <- length(idx$names)

 if (is.null(basis)) {
 return(.dkge_default_design_basis(fit,
 include_intercept = include_intercept,
 coding = coding,
 normalize = normalize))
 }

 if (is.list(basis) && !is.null(basis[["basis"]])) basis <- basis[["basis"]]
 if (is.list(basis) && !is.null(basis[["C"]])) basis <- basis[["C"]]
 if (is.list(basis) && is.null(dim(basis))) {
   stop("`basis` must have one row per DKGE effect/design cell.")
 }
 term_attr <- attr(basis, "term", exact = TRUE)
 source_attr <- attr(basis, "source", exact = TRUE)
 basis <- as.matrix(basis)
 if (nrow(basis) != q) {
 stop("`basis` must have one row per DKGE effect/design cell.")
 }

 if (!is.null(rownames(basis))) {
 missing <- setdiff(idx$names, rownames(basis))
 if (length(missing)) {
 stop("`basis` row names do not cover fit effects: ",
 paste(missing, collapse = ", "))
 }
 basis <- basis[idx$names, , drop = FALSE]
 } else {
 rownames(basis) <- idx$names
 }
 if (is.null(colnames(basis))) colnames(basis) <- paste0("contrast", seq_len(ncol(basis)))
 if (is.null(term_attr)) term_attr <- colnames(basis)
 if (is.null(source_attr)) source_attr <- "custom"
 attr(basis, "term") <- as.character(term_attr)
 attr(basis, "source") <- as.character(source_attr)
 basis
}

#' Build a plotting contrast basis for DKGE effects
#'
#' Returns a q-by-m matrix whose columns define interpretable coordinates for
#' plotting DKGE saliences. With cell-space design metadata, the default is a
#' factorial basis with a grand-mean column plus main effects/interactions. With
#' declared effect-space metadata or no recoverable design metadata, the default
#' is the identity basis over effect rows. Declared cell metadata is matched by
#' unique names; incomplete or ambiguous mappings are errors.
#'
#' @param fit Fitted `dkge` object.
#' @param basis Optional custom q-by-m matrix. Row names, when present, are
#'   matched to `fit$effects`.
#' @param include_intercept Include a grand-mean column for automatically built
#'   factorial bases.
#' @param coding Contrast coding for automatically built factorial bases:
#'   `"sum"`, `"helmert"`, or `"poly"`.
#' @param normalize Column scaling applied to automatically built bases,
#'   including the identity fallback used when the fit carries no design-cell
#'   metadata: `"unit_K"` (default) gives every column unit K-norm,
#'   \eqn{\sqrt{c^\top K c} = 1}; `"unit_l2"` gives unit Euclidean norm; and
#'   `"none"` leaves the coding matrix unscaled. User-supplied `basis` matrices
#'   are never rescaled.
#' @details
#' Basis columns are scored against saliences through \eqn{C^\top K U}, so the
#' comparability of coordinates across columns is governed by the K metric, not
#' by Euclidean length. Normalizing to unit Euclidean norm (`"unit_l2"`) leaves
#' columns with different K-norms — a grand-mean column is typically inflated
#' relative to interaction columns purely by the metric — which is why
#' `"unit_K"` is the default.
#' @return A numeric matrix with attributes `term` and `source`.
#' @examples
#' dk <- design_kernel(list(task = list(L = 2), measure = list(L = 3)),
#'                     basis = "cell", normalize = "none")
#' labels <- rownames(dk$K)
#' q <- length(labels)
#' set.seed(1)
#' B_list <- replicate(3, matrix(rnorm(q * 6), q, 6,
#'                               dimnames = list(labels, NULL)), simplify = FALSE)
#' X_list <- replicate(3, matrix(diag(q), q, q,
#'                               dimnames = list(labels, labels)), simplify = FALSE)
#' fit <- dkge(B_list, X_list, K = dk, rank = 2)
#' C <- dkge_design_basis(fit)
#' attr(C, "term")
#' # Columns have unit K-norm by default.
#' round(diag(t(C) %*% fit$K %*% C), 12)
#' @export
dkge_design_basis <- function(fit,
                              basis = NULL,
                              include_intercept = TRUE,
                              coding = c("sum", "helmert", "poly"),
                              normalize = c("unit_K", "unit_l2", "none")) {
 stopifnot(inherits(fit, "dkge"))
 coding <- match.arg(coding)
 normalize <- match.arg(normalize)
 .dkge_validate_design_basis(fit,
 basis = basis,
 include_intercept = include_intercept,
 coding = coding,
 normalize = normalize)
}

#' Express DKGE component saliences in a contrast basis
#'
#' Computes \eqn{C' K U}, where columns of `C` are user-supplied or inferred
#' design contrasts. This is a coordinate summary of the raw saliences, not a
#' replacement for inspecting `K %*% U` directly.
#'
#' @inheritParams dkge_component_saliences
#' @inheritParams dkge_design_basis
#' @return Tidy data frame with one row per component and contrast-basis column.
#' @examples
#' dk <- design_kernel(list(task = list(L = 2), measure = list(L = 3)),
#'                     basis = "cell", normalize = "none")
#' labels <- rownames(dk$K)
#' q <- length(labels)
#' set.seed(1)
#' B_list <- replicate(3, matrix(rnorm(q * 6), q, 6,
#'                               dimnames = list(labels, NULL)), simplify = FALSE)
#' X_list <- replicate(3, matrix(diag(q), q, q,
#'                               dimnames = list(labels, labels)), simplify = FALSE)
#' fit <- dkge(B_list, X_list, K = dk, rank = 2)
#' head(dkge_component_contrast_scores(fit, comps = 1:2))
#' @export
dkge_component_contrast_scores <- function(fit,
                                           basis = NULL,
                                           comps = NULL,
                                           include_intercept = TRUE,
                                           coding = c("sum", "helmert", "poly"),
                                           normalize = c("unit_K", "unit_l2", "none")) {
 stopifnot(inherits(fit, "dkge"))
 comps <- .dkge_component_indexer(fit, comps)
 C <- dkge_design_basis(fit,
 basis = basis,
 include_intercept = include_intercept,
 coding = coding,
 normalize = normalize)
 S <- dkge_component_saliences(fit, comps = comps, scale = "raw", long = FALSE)
 scores <- crossprod(C, S)
 terms <- attr(C, "term") %||% colnames(C)
 terms <- rep(as.character(terms), length.out = nrow(scores))

 data.frame(
 term = rep(terms, times = length(comps)),
 contrast = rep(rownames(scores), times = length(comps)),
 component = factor(rep(.dkge_component_labels(comps), each = nrow(scores)),
 levels = .dkge_component_labels(comps)),
 component_id = rep(comps, each = nrow(scores)),
 score = as.vector(scores),
 basis_source = attr(C, "source") %||% "custom",
 stringsAsFactors = FALSE
 )
}

#' DKGE scree plot with cumulative curve
#'
#' @param fit Fitted `dkge` object.
#' @param one_se_pick Optional integer component chosen by one-SE rule.
#' @return A ggplot object.
#' @examples
#' toy <- dkge_sim_toy(
#'   factors = list(A = list(L = 2), B = list(L = 3)),
#'   active_terms = c("A", "B"), S = 3, P = 15, snr = 5
#' )
#' fit <- dkge(toy$B_list, toy$X_list, K = toy$K, rank = 2)
#' dkge_plot_scree(fit)
#' @export
dkge_plot_scree <- function(fit, one_se_pick = NULL) {
 stopifnot(inherits(fit, "dkge"))
 tab <- dkge_variance_explained(fit)
 tab$cumulative <- cumsum(tab$prop_var)

 p <- ggplot2::ggplot(tab, ggplot2::aes(x = component)) +
 ggplot2::geom_col(ggplot2::aes(y = prop_var), fill = "#4C78A8", width = 0.7) +
 ggplot2::geom_point(ggplot2::aes(y = cumulative), size = 2.1) +
 ggplot2::geom_line(ggplot2::aes(y = cumulative), linewidth = 0.7) +
 ggplot2::labs(title = "DKGE scree", subtitle = "Bars: variance proportion; line: cumulative",
 x = "Component", y = "Proportion of variance") +
 theme_dkge()

 if (!is.null(one_se_pick) && is.finite(one_se_pick) && one_se_pick >= 1) {
 ymax <- max(tab$cumulative)
 p <- p +
 ggplot2::annotate("rect", xmin = one_se_pick - 0.5, xmax = one_se_pick + 0.5,
 ymin = -Inf, ymax = Inf, alpha = 0.08, fill = "#E45756") +
 ggplot2::annotate("segment", x = one_se_pick, xend = one_se_pick,
 y = 0, yend = ymax, colour = "#E45756",
 linewidth = 0.8, linetype = 2) +
 ggplot2::annotate("text", x = one_se_pick, y = ymax,
 label = "one-SE pick", vjust = -0.6, colour = "#E45756")
 }
 p
}

#' Effect-space loadings heatmap
#'
#' Displays the design-kernel-weighted component basis \eqn{K U}.
#'
#' @param fit Fitted `dkge` object.
#' @param comps Components to include (defaults to first min(rank,6)).
#' @param zscore Logical; z-score loadings within each effect.
#' @return A ggplot object.
#' @examples
#' toy <- dkge_sim_toy(
#'   factors = list(A = list(L = 2), B = list(L = 3)),
#'   active_terms = c("A", "B"), S = 3, P = 15, snr = 5
#' )
#' fit <- dkge(toy$B_list, toy$X_list, K = toy$K, rank = 2)
#' dkge_plot_effect_loadings(fit, comps = 1:2)
#' @export
dkge_plot_effect_loadings <- function(fit, comps = NULL, zscore = FALSE) {
 stopifnot(inherits(fit, "dkge"))
 r <- ncol(fit$U)
 if (is.null(comps)) comps <- seq_len(min(r, 6L))
 comps <- comps[comps >= 1 & comps <= r]
 if (!length(comps)) stop("No valid components selected.")

 idx <- .dkge_effect_indexer(fit)
 eff_names <- idx$names
 blocks <- idx$blocks

 load_mat <- fit$K %*% fit$U[, comps, drop = FALSE]
 df <- data.frame(
 effect = rep(eff_names, times = length(comps)),
 component = factor(rep(paste0("comp", comps), each = length(eff_names)),
 levels = paste0("comp", comps)),
 loading = as.vector(load_mat)
 )

 if (zscore) {
 df <- within(df, {
 effect_f <- effect
 loading <- stats::ave(loading, effect_f, FUN = function(x) {
 s <- stats::sd(x)
 if (!is.finite(s) || s < 1e-12) return(rep(0, length(x)))
 (x - mean(x)) / s
 })
 })
 }

 df$effect <- factor(df$effect, levels = eff_names)

 p <- ggplot2::ggplot(df, ggplot2::aes(x = component, y = effect, fill = loading)) +
 ggplot2::geom_tile() +
 .dkge_scale_fill_diverging() +
 ggplot2::labs(title = "Design-space loadings (K U)",
 x = "Component", y = "Effect / contrast") +
 theme_dkge()

 if (length(blocks) > 1) {
 block_order <- order(vapply(blocks, function(b) min(as.integer(b)), integer(1)))
 blocks <- blocks[block_order]
 sep <- cumsum(vapply(blocks, length, integer(1)))
 sep <- sep[-length(sep)]
 if (length(sep)) {
 p <- p + ggplot2::geom_hline(yintercept = sep + 0.5, colour = "grey70", linewidth = 0.3)
 }
 }
 p
}

#' Plot raw DKGE component saliences
#'
#' @inheritParams dkge_component_saliences
#' @param type Plot type. `"heatmap"` is the most general view; `"profile"`
#'   draws one line per component across ordered effect rows.
#' @return A ggplot object.
#' @examples
#' toy <- dkge_sim_toy(
#'   factors = list(A = list(L = 2), B = list(L = 3)),
#'   active_terms = c("A", "B"), S = 3, P = 15, snr = 5
#' )
#' fit <- dkge(toy$B_list, toy$X_list, K = toy$K, rank = 2)
#' dkge_plot_component_saliences(fit, comps = 1:2)
#' @export
dkge_plot_component_saliences <- function(fit,
                                          comps = NULL,
                                          scale = c("raw", "unit", "zscore"),
                                          type = c("heatmap", "profile")) {
 scale <- match.arg(scale)
 type <- match.arg(type)
 comps <- .dkge_component_indexer(fit, comps)
 df <- dkge_component_saliences(fit, comps = comps, scale = scale, long = TRUE)
 idx <- .dkge_effect_indexer(fit)
 df$effect <- factor(df$effect, levels = idx$names)

 if (type == "profile") {
 return(ggplot2::ggplot(df, ggplot2::aes(x = effect, y = salience,
 group = component, colour = component)) +
 ggplot2::geom_hline(yintercept = 0, colour = "grey70", linetype = 2) +
 ggplot2::geom_line(linewidth = 0.7, alpha = 0.85) +
 ggplot2::geom_point(size = 1.8) +
 ggplot2::scale_colour_brewer(palette = "Dark2") +
 ggplot2::labs(title = "Raw component saliences",
 subtitle = if (scale == "raw") "K U" else paste("Display scale:", scale),
 x = "Effect / design cell", y = "Salience") +
 theme_dkge() +
 ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1)))
 }

 ggplot2::ggplot(df, ggplot2::aes(x = component, y = effect, fill = salience)) +
 ggplot2::geom_tile() +
 .dkge_scale_fill_diverging() +
 ggplot2::labs(title = "Raw component saliences",
 subtitle = if (scale == "raw") "K U" else paste("Display scale:", scale),
 x = "Component", y = "Effect / design cell") +
 theme_dkge()
}

#' Plot DKGE component saliences in a contrast basis
#'
#' @inheritParams dkge_component_contrast_scores
#' @param type Plot type: `"heatmap"` for compact comparison, `"bars"` for a
#'   component-by-component coordinate plot.
#' @return A ggplot object.
#' @examples
#' dk <- design_kernel(list(task = list(L = 2), measure = list(L = 3)),
#'                     basis = "cell", normalize = "none")
#' labels <- rownames(dk$K)
#' q <- length(labels)
#' set.seed(1)
#' B_list <- replicate(3, matrix(rnorm(q * 6), q, 6,
#'                               dimnames = list(labels, NULL)), simplify = FALSE)
#' X_list <- replicate(3, matrix(diag(q), q, q,
#'                               dimnames = list(labels, labels)), simplify = FALSE)
#' fit <- dkge(B_list, X_list, K = dk, rank = 2)
#' dkge_plot_component_contrast_scores(fit, comps = 1:2)
#' @export
dkge_plot_component_contrast_scores <- function(fit,
                                                basis = NULL,
                                                comps = NULL,
                                                include_intercept = TRUE,
                                                coding = c("sum", "helmert", "poly"),
                                                normalize = c("unit_K", "unit_l2", "none"),
                                                type = c("heatmap", "bars")) {
 type <- match.arg(type)
 scores <- dkge_component_contrast_scores(fit,
 basis = basis,
 comps = comps,
 include_intercept = include_intercept,
 coding = coding,
 normalize = normalize)
 scores$contrast <- factor(scores$contrast, levels = unique(scores$contrast))

 if (type == "bars") {
 return(ggplot2::ggplot(scores, ggplot2::aes(x = contrast, y = score, fill = component)) +
 ggplot2::geom_hline(yintercept = 0, colour = "grey70", linetype = 2) +
 ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8), width = 0.72) +
 ggplot2::scale_fill_brewer(palette = "Dark2") +
 ggplot2::facet_wrap(stats::as.formula("~ term"), scales = "free_x") +
 ggplot2::labs(title = "Component saliences in contrast basis",
 x = "Contrast basis column", y = "Coordinate") +
 theme_dkge() +
 ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1)))
 }

 ggplot2::ggplot(scores, ggplot2::aes(x = contrast, y = component, fill = score)) +
 ggplot2::geom_tile() +
 .dkge_scale_fill_diverging() +
 ggplot2::facet_grid(stats::as.formula(". ~ term"), scales = "free_x", space = "free_x") +
 ggplot2::labs(title = "Component saliences in contrast basis",
 x = "Contrast basis column", y = "Component") +
 theme_dkge() +
 ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1))
}

.dkge_subject_projection_groups <- function(groups, subjects) {
 subjects <- as.character(subjects)
 if (is.null(groups)) {
 out <- rep("all", length(subjects))
 names(out) <- subjects
 return(out)
 }
 if (is.data.frame(groups)) {
 subject_col <- intersect(c("subject", "subject_id", "id"), names(groups))[1]
 group_col <- intersect(c("group", "condition", "label"), names(groups))[1]
 if (is.na(subject_col) || is.na(group_col)) {
 stop("`groups` data frame must contain subject/subject_id/id and group/condition/label columns.",
 call. = FALSE)
 }
 keys <- as.character(groups[[subject_col]])
 if (anyDuplicated(keys)) {
 stop("`groups` contains duplicate subject identifiers: ",
 paste(unique(keys[duplicated(keys)]), collapse = ", "), call. = FALSE)
 }
 vals <- as.character(groups[[group_col]])
 names(vals) <- keys
 missing <- setdiff(subjects, keys)
 if (length(missing)) {
 stop("`groups` does not cover subject(s): ", paste(missing, collapse = ", "),
 call. = FALSE)
 }
 return(vals[subjects])
 }
 nms <- names(groups)
 vals <- as.character(groups)
 names(vals) <- nms
 if (!is.null(nms)) {
 if (anyDuplicated(nms)) {
 stop("`groups` contains duplicate subject names: ",
 paste(unique(nms[duplicated(nms)]), collapse = ", "), call. = FALSE)
 }
 missing <- setdiff(subjects, nms)
 if (length(missing)) {
 stop("`groups` names do not cover subject(s): ", paste(missing, collapse = ", "),
 call. = FALSE)
 }
 return(vals[subjects])
 }
 if (length(vals) != length(subjects)) {
 stop("`groups` must be NULL, a vector with one value per subject, or a subject/group data frame.",
 call. = FALSE)
 }
 names(vals) <- subjects
 vals
}

# Resolve the stored voxel weights for one subject. This defers entirely to
# `.dkge_subject_voxel_weights()` (R/dkge-moments.R), which is what
# `dkge_fit()` itself uses, so the plotting path and the fitting path agree:
# a per-subject list is indexed, a cohort vector of the right length is used
# as-is, a constant cohort vector is re-expanded, and a *varying* vector of the
# wrong length is an error rather than being recycled across a subject whose
# parcels it does not describe.
.dkge_fit_subject_voxel_weights <- function(fit, s, Bts, subject_id = NULL) {
 spec <- fit$voxel_weights_subject %||% fit$voxel_weights
 .dkge_subject_voxel_weights(spec, s, ncol(Bts), subject_id = subject_id)
}

# Voxel weights enter the projection as a per-cluster reweighting of the mean,
# so their overall scale must not matter: a uniform weight of 4 would otherwise
# double every score. Weights produced by `dkge_fit()` are already mean-1 (see
# `.dkge_combine_weights()`); normalising here makes hand-set weights obey the
# same convention.
.dkge_apply_voxel_weights <- function(Bts, w_s) {
 if (is.null(w_s) || length(w_s) == 0L) return(Bts)
 w <- pmax(as.numeric(w_s), 0)
 scale <- mean(w)
 if (is.finite(scale) && scale > 0) w <- w / scale
 sweep(Bts, 2L, sqrt(w), "*")
}

#' Project subjects onto DKGE components
#'
#' Computes one **signed** scalar per subject and component: the subject's mean
#' cluster brain score along that component.
#'
#' @details
#' For subject \eqn{s} with whitened, voxel-weighted betas \eqn{\tilde B_s}
#' (\eqn{q \times P_s}) and group component \eqn{u_j}, the per-cluster brain
#' scores are the columns of \eqn{a_{s,j} = \tilde B_s^\top K u_j \in R^{P_s}}.
#' The reported projection is their average over clusters,
#' \deqn{\pi_{s,j} = \frac{1}{P_s} \sum_{p} a_{s,j}[p]
#'                 = \langle K u_j, \bar b_s \rangle,
#'       \qquad \bar b_s = \frac{1}{P_s}\sum_p \tilde B_s[, p].}
#' Equivalently it is the K-inner product between the component salience
#' \eqn{K u_j} and the subject's cluster-averaged effect-space profile.
#'
#' This quantity is **linear** in the subject's betas, so it is signed and
#' \eqn{\pi_{s,j}} flips sign when \eqn{\tilde B_s} is negated. That is the
#' point of the function: it complements the sign-blind quadratic participation
#' measure \eqn{\lVert \tilde B_s^\top K u_j \rVert^2} reported by
#' [dkge_plot_subject_contrib()], which cannot distinguish a subject expressing
#' a component from a subject expressing its mirror image. Averaging over
#' clusters is what makes the score subject-independent in cluster space and
#' therefore comparable across subjects with different parcel counts; a subject
#' whose cluster scores are large but evenly split in sign will have a small
#' projection and a large energy.
#'
#' In `"loso"` mode the same definition is applied with the basis
#' \eqn{U^{(-s)}} estimated without subject \eqn{s}, using the fold loaders'
#' voxel weighting. When `align = TRUE` each fold basis is rotated onto the
#' *pooled* basis `fit$U` by K-Procrustes ([dkge_procrustes_K()]), so component
#' \eqn{j} names the same pooled axis for every subject. (The fold builder's own
#' alignment uses fold 1 as its reference, not `fit$U`, and is deliberately not
#' reused here.) With `align = FALSE` component signs and order are only defined
#' up to each fold's own eigen-solve, so cross-subject comparison is not
#' meaningful.
#'
#' Voxel weights are normalized to mean 1 and applied as \eqn{\sqrt{w}} column
#' scales — the same convention as the training blocks — so a uniform weight
#' of any magnitude leaves the projection unchanged. This is not the
#' column-wise weighted mean \eqn{\sum_p w_p b_p / \sum_p w_p}; that would
#' change the scale relative to [dkge_plot_subject_contrib()] energy, which
#' uses the same \eqn{\sqrt{w}} reweighting. Weights produced by
#' [dkge_fit()] already have mean 1.
#'
#' @param fit Fitted `dkge` object.
#' @param groups Optional vector of group labels, or a data frame with subject
#'   and group columns. Named vectors are matched by subject name and must cover
#'   every subject; unnamed vectors are matched positionally.
#' @param mode `"loso"` for held-out supplementary projections or `"pooled"` for
#'   descriptive projections on the pooled fit.
#' @param comps Components to include.
#' @param align Rotate each LOSO fold basis onto the pooled component axes
#'   (`fit$U`) by K-Procrustes before scoring.
#' @param ridge Ridge passed to the LOSO fold basis builder.
#' @return Tidy data frame with subject, group, component, component_id,
#'   projection, and mode.
#' @examples
#' toy <- dkge_sim_toy(
#'   factors = list(A = list(L = 2), B = list(L = 3)),
#'   active_terms = c("A", "B"), S = 4, P = 15, snr = 5
#' )
#' fit <- dkge(toy$B_list, toy$X_list, K = toy$K, rank = 2)
#' dkge_subject_component_projections(fit, mode = "pooled", comps = 1:2)
#' @export
dkge_subject_component_projections <- function(fit,
                                               groups = NULL,
                                               mode = c("loso", "pooled"),
                                               comps = NULL,
                                               align = TRUE,
                                               ridge = 0) {
 stopifnot(inherits(fit, "dkge"))
 mode <- match.arg(mode)
 comps <- .dkge_component_indexer(fit, comps)
 S <- length(fit$Btil)
 if (mode == "loso" && S < 2L) {
 stop("LOSO subject projections require at least two subjects.")
 }
 subjects <- as.character(fit$subject_ids %||% paste0("sub", seq_len(S)))
 group_vals <- .dkge_subject_projection_groups(groups, subjects)

 n_comp <- length(comps)
 proj <- matrix(NA_real_, nrow = S, ncol = n_comp)

 if (mode == "pooled") {
 saliences <- fit$K %*% fit$U[, comps, drop = FALSE]
 bbar <- matrix(NA_real_, nrow = nrow(fit$U), ncol = S)
 for (s in seq_len(S)) {
 Bts <- fit$Btil[[s]]
 w_s <- .dkge_fit_subject_voxel_weights(fit, s, Bts, subject_id = subjects[[s]])
 bbar[, s] <- rowMeans(.dkge_apply_voxel_weights(Bts, w_s))
 }
 proj <- crossprod(bbar, saliences)
 } else {
 assignments <- lapply(seq_len(S), function(s) s)
 fold_info <- .dkge_build_fold_bases(fit,
 assignments = assignments,
 ridge = ridge,
 align = FALSE,
 loader_scope = "heldout",
 verbose = FALSE,
 missingness = fit$missingness %||% "none",
 miss_args = fit$miss_args %||% list())
 for (fold in fold_info$folds) {
 s <- fold$subjects[[1]]
 loader <- fold$loaders[[as.character(s)]]
 # `.dkge_build_fold_bases(align = TRUE)` aligns folds to fold 1, not to
 # the pooled basis, so its rotation would make component j mean
 # "whatever fold 1 called j". Rotate onto `fit$U` directly instead, which
 # is what these projections are documented to be expressed in.
 R_align <- if (isTRUE(align)) {
 dkge_procrustes_K(fit$U, fold$basis, fit$K)$R
 } else {
 diag(ncol(fit$U))
 }
 A_aligned <- loader$A %*% R_align
 proj[s, ] <- colMeans(A_aligned[, comps, drop = FALSE])
 }
 }

 comp_labels <- .dkge_component_labels(comps)
 out <- data.frame(
 subject = rep(subjects, each = n_comp),
 group = rep(unname(group_vals), each = n_comp),
 component = factor(rep(comp_labels, times = S), levels = comp_labels),
 component_id = rep(comps, times = S),
 projection = as.vector(t(proj)),
 mode = mode,
 stringsAsFactors = FALSE
 )
 out$group <- factor(out$group)
 out
}

#' Plot subject projections onto DKGE components
#'
#' @inheritParams dkge_subject_component_projections
#' @param projections Optional data frame previously returned by
#'   [dkge_subject_component_projections()]. Supplying it skips recomputation,
#'   which matters for `mode = "loso"` because that path refits one basis per
#'   subject. It must contain `subject`, `group`, `component`, and a numeric
#'   `projection`. When supplied, `fit`, `groups`, `comps`, `align`, `ridge`,
#'   and `mode` are all ignored: the panel label is read from the frame's own
#'   `mode` column, and a frame without one (or mixing modes) is left
#'   unlabeled rather than being captioned with the `mode` default.
#' @return A ggplot object.
#' @examples
#' toy <- dkge_sim_toy(
#'   factors = list(A = list(L = 2), B = list(L = 3)),
#'   active_terms = c("A", "B"), S = 4, P = 15, snr = 5
#' )
#' fit <- dkge(toy$B_list, toy$X_list, K = toy$K, rank = 2)
#' proj <- dkge_subject_component_projections(fit, mode = "pooled", comps = 1:2)
#' dkge_plot_subject_component_projections(fit, projections = proj)
#' @export
dkge_plot_subject_component_projections <- function(fit = NULL,
                                                    groups = NULL,
                                                    mode = c("loso", "pooled"),
                                                    comps = NULL,
                                                    align = TRUE,
                                                    ridge = 0,
                                                    projections = NULL) {
 mode <- match.arg(mode)
 if (is.null(fit) && is.null(projections)) {
 stop("Supply either `fit` or a precomputed `projections` data frame.", call. = FALSE)
 }
 if (!is.null(projections)) {
 df <- as.data.frame(projections)
 required <- c("subject", "group", "component", "projection")
 missing_cols <- setdiff(required, names(df))
 if (length(missing_cols)) {
 stop("`projections` must contain columns: ", paste(missing_cols, collapse = ", "),
 call. = FALSE)
 }
 if (!nrow(df)) {
 stop("`projections` must have at least one row.", call. = FALSE)
 }
 if (!is.numeric(df$projection)) {
 stop("`projections$projection` must be numeric.", call. = FALSE)
 }
 # The label is a claim about how the numbers were produced, so it is read
 # from the frame rather than from the ignored `mode` argument: a frame with
 # no `mode` column was previously captioned "LOSO" (the `mode` default) even
 # when it held pooled projections.
 mode <- .dkge_projection_frame_mode(df)
 } else {
 df <- dkge_subject_component_projections(fit,
 groups = groups,
 mode = mode,
 comps = comps,
 align = align,
 ridge = ridge)
 }
 ggplot2::ggplot(df, ggplot2::aes(x = component, y = projection, colour = group)) +
 ggplot2::geom_hline(yintercept = 0, colour = "grey70", linetype = 2) +
 ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.12, height = 0, seed = 1),
 size = 2.2, alpha = 0.9) +
 ggplot2::scale_colour_brewer(palette = "Set1", na.value = "grey35") +
 ggplot2::labs(title = .dkge_projection_mode_title(mode),
 subtitle = .dkge_projection_mode_subtitle(mode),
 x = "Component", y = "Mean cluster brain score", colour = "Group") +
 theme_dkge()
}

# A supplied projections frame is labelled by what it actually records. A frame
# with no `mode` column, or one mixing modes, cannot be claimed to be either.
.dkge_projection_frame_mode <- function(df) {
 if (is.null(df$mode)) return("unknown")
 modes <- unique(as.character(df$mode))
 modes <- modes[!is.na(modes) & nzchar(modes)]
 if (length(modes) == 1L) modes else "unknown"
}

.dkge_projection_mode_title <- function(mode) {
 if (identical(mode, "unknown")) {
 return("Subject projections on DKGE components")
 }
 paste0("Subject projections on DKGE components (", toupper(mode), ")")
}

.dkge_projection_mode_subtitle <- function(mode) {
 switch(mode,
 loso = "Signed mean cluster brain score, held-out basis aligned to the pooled axes",
 pooled = "Signed mean cluster brain score on the pooled component axes",
 "Signed mean cluster brain score; component basis not recorded in the supplied projections")
}

#' Subject weights and per-component energy heatmap
#'
#' @param fit Fitted `dkge` object.
#' @param comps Components to display (default first min(rank,6)).
#' @return List with `weights` and `energy` ggplots.
#' @examples
#' toy <- dkge_sim_toy(
#'   factors = list(A = list(L = 2), B = list(L = 3)),
#'   active_terms = c("A", "B"), S = 3, P = 15, snr = 5
#' )
#' fit <- dkge(toy$B_list, toy$X_list, K = toy$K, rank = 2)
#' panels <- dkge_plot_subject_contrib(fit, comps = 1:2)
#' panels$weights
#' @export
dkge_plot_subject_contrib <- function(fit, comps = NULL) {
 stopifnot(inherits(fit, "dkge"))
 r <- ncol(fit$U)
 if (is.null(comps)) comps <- seq_len(min(r, 6L))
 comps <- comps[comps >= 1 & comps <= r]
 if (!length(comps)) stop("No valid components selected.")

 weights <- fit$weights %||% rep(1, length(fit$Btil))
 subjects <- fit$subject_ids %||% paste0("sub", seq_along(weights))
 df_weights <- data.frame(subject = factor(subjects, levels = subjects), weight = weights)

 energy <- .dkge_energy_by_subject_component(fit)
 energy <- energy[, comps, drop = FALSE]
 energy_df <- data.frame(
 subject = factor(rep(subjects, times = length(comps)), levels = subjects),
 component = factor(rep(paste0("comp", comps), each = length(subjects)),
 levels = paste0("comp", comps)),
 energy = as.vector(energy)
 )

 p_weights <- ggplot2::ggplot(df_weights, ggplot2::aes(x = subject, y = weight)) +
 ggplot2::geom_col(fill = "#72B7B2") +
 ggplot2::labs(title = "Subject weights", x = NULL, y = "Weight") +
 theme_dkge() +
 ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, hjust = 1, vjust = 1))

 p_energy <- ggplot2::ggplot(energy_df, ggplot2::aes(x = component, y = subject, fill = energy)) +
 ggplot2::geom_tile() +
 ggplot2::scale_fill_viridis_c() +
 ggplot2::labs(title = "Component participation (energy)", x = "Component", y = NULL) +
 theme_dkge()

 list(weights = p_weights, energy = p_energy)
}

#' Subspace stability via principal angles
#'
#' @param bases List of basis matrices.
#' @param K Design kernel.
#' @param consensus Optional consensus basis.
#' @param labels Optional labels.
#' @return A ggplot object.
#' @examples
#' toy <- dkge_sim_toy(
#'   factors = list(A = list(L = 2), B = list(L = 3)),
#'   active_terms = c("A", "B"), S = 3, P = 15, snr = 5
#' )
#' fit <- dkge(toy$B_list, toy$X_list, K = toy$K, rank = 2)
#' bases <- list(fit$U, fit$U)
#' dkge_plot_subspace_stability(bases, K = fit$K)
#' @export
dkge_plot_subspace_stability <- function(bases, K, consensus = NULL, labels = NULL) {
 stopifnot(is.list(bases), is.matrix(K))
 if (is.null(consensus)) {
 consensus <- dkge_consensus_basis_K(bases, K)$U
 } else if (is.list(consensus) && !is.null(consensus$U)) {
 consensus <- consensus$U
 }
 if (is.null(labels)) labels <- paste0("base", seq_along(bases))

 angle_df <- do.call(rbind, lapply(seq_along(bases), function(i) {
 theta <- .dkge_principal_angles_K(bases[[i]], consensus, K) * 180 / pi
 # `seq_along(theta)`, not `seq_len(r)`: a rank-deficient basis yields fewer
 # angles than it has columns.
 data.frame(base = labels[i], component = seq_along(theta), angle_deg = theta)
 }))

 ggplot2::ggplot(angle_df, ggplot2::aes(x = component, y = angle_deg, group = base, colour = base)) +
 ggplot2::geom_line(alpha = 0.7) +
 ggplot2::geom_point(size = 1.8) +
 ggplot2::scale_colour_brewer(palette = 'Set2') +
 ggplot2::labs(title = "Subspace stability", y = "Angle (degrees)", x = "Component") +
 theme_dkge() +
 ggplot2::guides(colour = "none")
}

#' Anchor-level information plots
#'
#' @param info_haufe Result from `dkge_info_map_haufe()`.
#' @param info_loco Result from `dkge_info_map_loco()`.
#' @param top Number of anchors to annotate.
#' @return List of ggplot objects.
#' @examples
#' info_haufe <- list(mean_anchor = rnorm(12))
#' panels <- dkge_plot_info_anchor(info_haufe = info_haufe, top = 3)
#' panels$haufe
#' @export
dkge_plot_info_anchor <- function(info_haufe = NULL, info_loco = NULL, top = 20) {
 if (is.null(info_haufe) && is.null(info_loco)) {
 stop("Provide at least one of info_haufe or info_loco.")
 }

 panels <- list()
 annotate_top <- function(p, df, mapping) {
 if (is.null(top) || top <= 0) return(p)
 df$rank <- rank(-abs(df$value), ties.method = "first")
 lab <- df[df$rank <= top, , drop = FALSE]
 if (requireNamespace("ggrepel", quietly = TRUE)) {
 p + ggrepel::geom_text_repel(data = lab, mapping, size = 3, max.overlaps = 100)
 } else {
 p + ggplot2::geom_text(data = lab, mapping, size = 3)
 }
 }

 if (!is.null(info_haufe)) {
 values <- as.numeric(info_haufe$mean_anchor %||% info_haufe$anchor %||% info_haufe$y)
 df <- data.frame(anchor = seq_along(values), value = values)
 p <- ggplot2::ggplot(df, ggplot2::aes(x = anchor, y = value, colour = value)) +
 ggplot2::geom_segment(ggplot2::aes(xend = anchor, y = 0, yend = value), linewidth = 0.6, alpha = 0.8) +
 .dkge_scale_colour_diverging() +
 ggplot2::labs(title = "Haufe anchor weights", x = "Anchor", y = "Weight") +
 theme_dkge()
 p <- annotate_top(p, df, ggplot2::aes(x = anchor, y = value, label = anchor))
 panels$haufe <- p
 }

 if (!is.null(info_loco)) {
 values <- as.numeric(info_loco$loco_anchor %||% info_loco$delta %||% info_loco$y)
 df <- data.frame(anchor = seq_along(values), value = values)
 p <- ggplot2::ggplot(df, ggplot2::aes(x = anchor, y = value, fill = value)) +
 ggplot2::geom_col(width = 0.8) +
 .dkge_scale_fill_diverging() +
 ggplot2::labs(title = "LOCO importance (anchors)", x = "Anchor", y = "Delta score") +
 theme_dkge()
 p <- annotate_top(p, df, ggplot2::aes(x = anchor, y = value, label = anchor))
 panels$loco <- p
 }

 panels
}
