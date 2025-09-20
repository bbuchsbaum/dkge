# dkge-components.R
# Convenience helpers for component-level inference and transport

#' Component-level consensus statistics
#'
#' Transports each subject's component loadings to a reference parcellation,
#' computes inference statistics, and returns tidy summaries ready for
#' visualisation.
#'
#' @param fit A fitted `dkge` object.
#' @param mapper Mapper strategy (string or [dkge_mapper_spec()]). Defaults to
#'   "sinkhorn".
#' @param centroids Optional list of subject centroid matrices; defaults to
#'   centroids stored in `fit` if available.
#' @param sizes Optional list of cluster masses (one vector per subject).
#' @param inference One of "signflip" or "parametric", or a list providing
#'   `type`, `B`, `tail`, and `alpha`.
#' @param medoid Reference subject index (defaults to 1).
#' @param components Optional vector of component indices or names; default is
#'   all components.
#' @param adjust Method supplied to [stats::p.adjust()] for multiple testing
#'   correction in the tidy summary.
#' @param ... Additional mapper-specific parameters (e.g. `epsilon`).
#'
#' @return A list with fields:
#'   - `summary`: tidy data frame of statistics and p-values.
#'   - `statistics`: per-component statistic vectors.
#'   - `transport`: per-component transported subject matrices.
#' @export
dkge_component_stats <- function(fit,
                                 mapper = "sinkhorn",
                                 centroids = NULL,
                                 sizes = NULL,
                                 inference = "signflip",
                                 medoid = 1L,
                                 components = NULL,
                                 adjust = "fdr",
                                 ...) {
  stopifnot(inherits(fit, "dkge"))

  centroids <- centroids %||% fit$centroids %||% fit$input$centroids %||%
    stop("Centroids must be supplied or stored in the fit object.")

  # Build mapper specification
  mapper_spec <- .dkge_resolve_mapper_spec(mapper, method = NULL, dots = list(...))

  loadings <- lapply(fit$Btil, function(Bts) t(Bts) %*% fit$K %*% fit$U)
  rank <- ncol(loadings[[1]])

  if (is.null(components)) {
    comp_idx <- seq_len(rank)
  } else if (is.numeric(components)) {
    comp_idx <- components
  } else {
    comp_idx <- match(components, colnames(fit$U))
  }

  transport <- dkge_transport_loadings_to_medoid(fit,
                                                 medoid = medoid,
                                                 centroids = centroids,
                                                 loadings = loadings,
                                                 sizes = sizes,
                                                 mapper = mapper_spec)

  subj_mats <- lapply(transport$subjects, function(M) M[, comp_idx, drop = FALSE])

  inference_res <- .dkge_component_inference(subj_mats, inference)
  tidy <- .dkge_component_tidy(inference_res, comp_idx, adjust)

  list(summary = tidy,
       statistics = inference_res$stats,
       transport = subj_mats)
}

#' @rdname dkge_component_stats
#' @param file Path to the CSV file where component statistics will be written.
#' @export
dkge_write_component_stats <- function(fit, file, ...) {
  res <- dkge_component_stats(fit, ...)
  utils::write.csv(res$summary, file = file, row.names = FALSE)
  invisible(res)
}

.dkge_component_inference <- function(subj_mats, inference) {
  if (is.character(inference)) {
    inference <- list(type = inference)
  }
  type <- inference$type %||% "signflip"
  alpha <- inference$alpha %||% 0.05

  stats <- vector("list", length(subj_mats))
  pvals <- vector("list", length(subj_mats))
  for (i in seq_along(subj_mats)) {
    Y <- subj_mats[[i]]
    if (type == "signflip") {
      B <- inference$B %||% 2000
      tail <- inference$tail %||% "two.sided"
      res <- dkge_signflip_maxT(Y, B = B, tail = tail)
      stats[[i]] <- res$stat
      pvals[[i]] <- res$p
    } else if (type == "parametric") {
      mu <- colMeans(Y)
      se <- apply(Y, 2, stats::sd) / sqrt(nrow(Y))
      tstat <- mu / (se + 1e-12)
      stats[[i]] <- tstat
      pvals[[i]] <- 2 * stats::pt(-abs(tstat), df = nrow(Y) - 1)
    } else {
      stop("Unsupported inference type")
    }
  }
  list(stats = stats, pvals = pvals, alpha = alpha)
}

.dkge_component_tidy <- function(inference_res, comp_idx, adjust) {
  out <- Map(function(stat, p, comp) {
    data.frame(component = comp,
               cluster = seq_along(stat),
               stat = stat,
               p = p,
               stringsAsFactors = FALSE)
  }, inference_res$stats, inference_res$pvals, comp_idx)
  df <- do.call(rbind, out)
  df$p_adj <- stats::p.adjust(df$p, method = adjust)
  df$significant <- df$p_adj <= inference_res$alpha
  df
}
