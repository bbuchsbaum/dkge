#' DKGE "Five Fundamentals" dashboard
#'
#' Layout helper combining the main DKGE diagnostic plots into a single
#' patchwork canvas. Individual panels degrade gracefully when optional inputs
#' (e.g. LOSO bases or information maps) are missing.
#'
#' @inheritParams dkge_plot_scree
#' @inheritParams dkge_plot_effect_loadings
#' @inheritParams dkge_plot_subject_contrib
#' @param bases Optional list of basis matrices (e.g. LOSO/fold bases) for the
#'   stability panel.
#' @param consensus Optional consensus basis for stability plotting.
#' @param base_labels Optional labels for the supplied bases.
#' @param info_haufe Optional result from [dkge_info_map_haufe()].
#' @param info_loco Optional result from [dkge_info_map_loco()].
#' @param top Number of anchors to annotate in information plots (set to 0 to
#'   disable).
#' @param width,height Dimensions (inches) when saving to disk.
#' @param dpi Resolution (dots per inch) when saving.
#' @param save_path Optional file path (png/pdf/svg) to save the dashboard.
#' @return A patchwork object (invisibly if saved).
#' @export
dkge_plot_suite <- function(fit,
                            one_se_pick = NULL,
                            comps = NULL,
                            zscore = FALSE,
                            bases = NULL,
                            consensus = NULL,
                            base_labels = NULL,
                            info_haufe = NULL,
                            info_loco = NULL,
                            top = 20,
                            width = 12,
                            height = 12,
                            dpi = 300,
                            save_path = NULL) {
  if (!inherits(fit, "dkge")) stop("`fit` must be a dkge object.")
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required for dkge_plot_suite().")
  }

  placeholder <- function(msg) {
    ggplot2::ggplot() + ggplot2::theme_void() +
      ggplot2::annotate("text", x = 0.5, y = 0.5, label = msg, size = 4)
  }

  p_scree <- tryCatch(dkge_plot_scree(fit, one_se_pick = one_se_pick),
                      error = function(e) placeholder("Scree unavailable"))

  p_effect <- tryCatch(dkge_plot_effect_loadings(fit, comps = comps, zscore = zscore),
                       error = function(e) placeholder("Effect loadings unavailable"))

  subj_panels <- tryCatch(dkge_plot_subject_contrib(fit, comps = comps), error = identity)
  if (inherits(subj_panels, "error")) {
    p_weights <- placeholder("Subject weights unavailable")
    p_energy <- placeholder("Participation heatmap unavailable")
  } else {
    p_weights <- subj_panels$weights
    p_energy <- subj_panels$energy
  }

  if (is.null(bases)) {
    p_stab <- placeholder("Provide `bases` for stability plot")
  } else {
    p_stab <- tryCatch(dkge_plot_subspace_stability(bases, K = fit$K,
                                                    consensus = consensus,
                                                    labels = base_labels),
                       error = function(e) placeholder("Stability plotting failed"))
  }

  if (is.null(info_haufe) && is.null(info_loco)) {
    p_info <- placeholder("Information maps not supplied")
  } else {
    info_panels <- tryCatch(dkge_plot_info_anchor(info_haufe = info_haufe,
                                                  info_loco = info_loco,
                                                  top = top),
                            error = identity)
    if (inherits(info_panels, "error")) {
      p_info <- placeholder("Information maps failed")
    } else {
      if (!is.null(info_panels$haufe) && !is.null(info_panels$loco)) {
        p_info <- info_panels$haufe + info_panels$loco + patchwork::plot_layout(ncol = 2)
      } else if (!is.null(info_panels$haufe)) {
        p_info <- info_panels$haufe
      } else {
        p_info <- info_panels$loco
      }
    }
  }

  design <- "
  AABBB
  CCDDD
  EEEEE
  FFFFF
  "

  dashboard <- p_scree + p_effect + p_weights + p_energy + p_stab + p_info +
    patchwork::plot_layout(design = design)

  if (!is.null(save_path)) {
    ggplot2::ggsave(filename = save_path, plot = dashboard,
                    width = width, height = height, dpi = dpi, limitsize = FALSE)
    return(invisible(dashboard))
  }
  dashboard
}
