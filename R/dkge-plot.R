#' @import ggplot2
NULL

# -----------------------------------------------------------------------------
# Shared aesthetics -----------------------------------------------------------
# -----------------------------------------------------------------------------

#' DKGE minimal theme for ggplot2 outputs
#'
#' Produces a light-weight theme used by the DKGE plotting helpers. Adjust
#' `base_size` / `base_family` to customise font size or typeface.
#'
#' @param base_size Base font size.
#' @param base_family Base font family.
#' @return A ggplot2 theme object.
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

.dkge_energy_by_subject_component <- function(fit) {
 r <- ncol(fit$U)
 S <- length(fit$Btil)
 energy <- matrix(NA_real_, nrow = S, ncol = r)
 for (s in seq_len(S)) {
 B_t <- fit$Btil[[s]]
 A_s <- t(B_t) %*% fit$K %*% fit$U
 energy[s, ] <- colSums(A_s * A_s)
 }
 rownames(energy) <- fit$subject_ids %||% paste0("sub", seq_len(S))
 colnames(energy) <- paste0("comp", seq_len(r))
 energy
}

.dkge_principal_angles_K <- function(U1, U2, K) {
 sv <- svd(t(U1) %*% K %*% U2, nu = 0, nv = 0)$d
 sv <- pmin(pmax(sv, -1), 1)
 acos(sv)
}

#' DKGE scree plot with cumulative curve
#'
#' @param fit Fitted `dkge` object.
#' @param one_se_pick Optional integer component chosen by one-SE rule.
#' @return A ggplot object.
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

#' Effect-space loadings heatmap (K \%*\% U)
#'
#' @param fit Fitted `dkge` object.
#' @param comps Components to include (defaults to first min(rank,6)).
#' @param zscore Logical; z-score loadings within each effect.
#' @return A ggplot object.
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
 ggplot2::labs(title = "Design-space loadings (K %*% U)",
 x = "Component", y = "Effect / contrast") +
 theme_dkge()

 if (length(blocks) > 1) {
 sep <- cumsum(vapply(blocks, length, integer(1)))
 sep <- sep[-length(sep)]
 if (length(sep)) {
 p <- p + ggplot2::geom_hline(yintercept = sep + 0.5, colour = "grey70", linewidth = 0.3)
 }
 }
 p
}

#' Subject weights and per-component energy heatmap
#'
#' @param fit Fitted `dkge` object.
#' @param comps Components to display (default first min(rank,6)).
#' @return List with `weights` and `energy` ggplots.
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
 subject = rep(subjects, times = length(comps)),
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
#' @export
dkge_plot_subspace_stability <- function(bases, K, consensus = NULL, labels = NULL) {
 stopifnot(is.list(bases), is.matrix(K))
 r <- ncol(bases[[1]])
 if (is.null(consensus)) {
 consensus <- dkge_consensus_basis_K(bases, K)$U
 } else if (is.list(consensus) && !is.null(consensus$U)) {
 consensus <- consensus$U
 }
 if (is.null(labels)) labels <- paste0("base", seq_along(bases))

 angle_df <- do.call(rbind, lapply(seq_along(bases), function(i) {
 theta <- .dkge_principal_angles_K(bases[[i]], consensus, K) * 180 / pi
 data.frame(base = labels[i], component = seq_len(r), angle_deg = theta)
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
