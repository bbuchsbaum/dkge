# dkge-effect-weights.R
# Subject-by-effect reliability specifications for q-space pooling.

#' Specify effect-level reliability weighting
#'
#' Effect weights control how unequally estimated subject-by-effect rows enter
#' the pooled second moment. They are distinct from [dkge_weights()], which
#' weights voxels or parcels, and from `w_method`, which supplies one scalar per
#' subject.
#'
#' @param method Effect precision source. `"none"` gives every observed effect
#'   equal precision, `"count"` uses `effect_n` stored on each subject, and
#'   `"precision"` uses explicitly supplied `effect_precision` values.
#'
#' @return An object of class `dkge_effect_weights`, suitable for the
#'   `effect_weights` argument of [dkge_fit()]. A precision of zero excludes that
#'   subject-effect row. Pairwise reliability is the geometric mean of the two
#'   effect precisions.
#' @export
#' @examples
#' dkge_effect_weights("none")
#' dkge_effect_weights("count")
dkge_effect_weights <- function(method = c("none", "count", "precision")) {
  method <- match.arg(method)
  structure(list(method = method), class = "dkge_effect_weights")
}

#' @noRd
.dkge_resolve_effect_precision <- function(dataset, spec, obs_masks = NULL,
                                           indices = seq_len(dataset$n_subjects)) {
  if (is.null(spec)) spec <- dkge_effect_weights("none")
  if (!inherits(spec, "dkge_effect_weights")) {
    stop("`effect_weights` must be created by `dkge_effect_weights()`.",
         call. = FALSE)
  }

  method <- spec$method %||% "none"
  q <- dataset$q
  effects <- dataset$effects
  source <- switch(
    method,
    none = NULL,
    count = dataset$effect_n,
    precision = dataset$effect_precision,
    stop("Unknown effect-weighting method.", call. = FALSE)
  )

  out <- vector("list", length(indices))
  for (j in seq_along(indices)) {
    s <- indices[[j]]
    mask <- if (is.null(obs_masks)) rep(TRUE, q) else as.logical(obs_masks[[s]])
    if (length(mask) != q) {
      stop(sprintf(
        "Observation mask for subject '%s' must have one entry per effect (got %d, expected %d).",
        dataset$subject_ids[[s]], length(mask), q
      ), call. = FALSE)
    }

    if (identical(method, "none")) {
      p <- as.numeric(mask)
    } else {
      p <- source[[s]]
      if (is.null(p)) {
        stop(sprintf(
          "Subject '%s' has no `%s` values required by effect_weights='%s'.",
          dataset$subject_ids[[s]],
          if (method == "count") "effect_n" else "effect_precision",
          method
        ), call. = FALSE)
      }
      p <- as.numeric(p)
      if (length(p) != q || any(!is.finite(p)) || any(p < 0)) {
        stop(sprintf(
          "Subject '%s' must have %d finite, non-negative effect precision values.",
          dataset$subject_ids[[s]], q
        ), call. = FALSE)
      }
      p[!mask] <- 0
    }
    names(p) <- effects
    out[[j]] <- p
  }
  names(out) <- dataset$subject_ids[indices]
  out
}
