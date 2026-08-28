# User-facing summary for fitted DKGE models.

#' Print a fitted DKGE model
#'
#' Reports the fitted dimensions, spatial-smoothing effectiveness, and the
#' subject-level block weighting used in the pooled moment. When selected,
#' MFA/energy normalization weights are not automatically inverse-variance
#' weights for group inference; `w_method = "none"` reports equal moment
#' weights.
#'
#' @param x A fitted `dkge` object.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.dkge <- function(x, ...) {
  if (!is.null(x$spatial)) {
    .dkge_spatial_fit_payload(x$spatial, validate = TRUE)
  }
  n_subjects <- length(x$Btil %||% x$subject_ids)
  n_effects <- if (!is.null(x$K)) nrow(x$K) else NA_integer_
  fitted_rank <- ncol(x$U %||% matrix(numeric(), 0L, 0L))
  weights <- as.numeric(x$weights %||% rep(1, n_subjects))
  finite <- is.finite(weights) & weights >= 0
  active <- finite & weights > 0

  cat("<dkge>\n", sep = "")
  cat("  Subjects:", n_subjects, "\n")
  cat("  Effects:", n_effects, "\n")
  cat("  Rank:", fitted_rank, "\n")
  if (!is.null(x$spatial)) {
    spatial_status <- x$spatial$status %||%
      if (isTRUE(x$spatial$active)) "active" else if (x$spatial$lambda > 0) {
        "inert"
      } else {
        "inactive"
      }
    edge_counts <- x$spatial$diagnostics$n_edges %||% numeric(0)
    edge_summary <- if (length(edge_counts)) {
      paste(format(edge_counts, trim = TRUE), collapse = ", ")
    } else {
      "unavailable"
    }
    cat(
      "  Spatial smoothing:", spatial_status,
      sprintf("(lambda = %s; edges = %s)",
              format(x$spatial$lambda), edge_summary),
      "\n"
    )
  }
  cat(
    "  Subject weighting:", x$w_method %||% "unknown",
    sprintf("(tau = %s)", format(x$w_tau %||% NA_real_, digits = 4)),
    "\n"
  )

  if (any(active)) {
    w <- weights[active]
    ess <- sum(w)^2 / sum(w^2)
    dispersion <- if (length(w) > 1L && mean(w) > 0) {
      stats::sd(w) / mean(w)
    } else {
      0
    }
    cat(
      "  Weight range:",
      paste(format(range(w), digits = 4), collapse = " to "),
      sprintf("(median = %s, CV = %s)",
              format(stats::median(w), digits = 4),
              format(dispersion, digits = 4)),
      "\n"
    )
    cat(
      "  Effective subject mass:", format(ess, digits = 4),
      sprintf("of %d usable", length(w)), "\n"
    )
  } else {
    cat("  Weight range: unavailable\n")
    cat("  Effective subject mass: 0\n")
  }

  invisible(x)
}
