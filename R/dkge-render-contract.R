# dkge-render-contract.R
# Rendering is deliberately downstream of correspondence and inference.

.dkge_renderer_payload <- function(x) {
  list(
    schema_version = x$schema_version,
    support_id = x$support_id,
    support_hash = x$support_hash,
    decoder_hash = x$decoder_hash,
    provenance = x$provenance
  )
}

.dkge_validate_renderer <- function(x) {
  if (!inherits(x, "dkge_renderer")) {
    .dkge_abort("Expected a `dkge_renderer` object.", "dkge_renderer_error")
  }
  .dkge_validate_reference_support(x$reference_support)
  if (!identical(x$support_id, x$reference_support$support_id) ||
      !identical(x$support_hash, x$reference_support$structural_hash) ||
      !identical(x$decoder_hash,
                 .dkge_object_hash(x$reference_support$decoder))) {
    .dkge_abort("Renderer and reference support identities do not match.",
                "dkge_renderer_error")
  }
  expected <- .dkge_object_hash(.dkge_renderer_payload(x))
  if (!identical(x$structural_hash, expected)) {
    .dkge_abort("Renderer was mutated after construction.",
                "dkge_renderer_error")
  }
  invisible(x)
}

#' Construct a renderer for an identified reference support
#'
#' A renderer never fits subject correspondence. It only attaches display
#' mechanics to an already identified support and consumes aligned maps or an
#' already computed statistic.
#'
#' @param support A [dkge_reference_support()] object.
#' @param provenance Optional rendering provenance.
#' @return An immutable `dkge_renderer` object.
#' @export
dkge_renderer <- function(support, provenance = NULL) {
  .dkge_validate_reference_support(support)
  out <- list(
    schema_version = "1.0.0",
    reference_support = support,
    support_id = support$support_id,
    support_hash = support$structural_hash,
    decoder_hash = .dkge_object_hash(support$decoder),
    provenance = provenance
  )
  out$structural_hash <- .dkge_object_hash(.dkge_renderer_payload(out))
  structure(out, class = c("dkge_renderer", "list"))
}

.dkge_renderer_statistic <- function(x,
                                     statistic = c("weighted_mean", "mean", "median")) {
  statistic <- match.arg(statistic)
  if (identical(statistic, "weighted_mean")) {
    weights <- x$subject_weights / sum(x$subject_weights)
    return(lapply(x$values, function(Y) as.numeric(crossprod(weights, Y))))
  }
  if (identical(statistic, "mean")) {
    return(lapply(x$values, colMeans))
  }
  lapply(x$values, apply, 2L, stats::median)
}

.dkge_normalize_render_values <- function(x, Q) {
  if (is.numeric(x) && is.null(dim(x))) x <- list(statistic = x)
  if (is.matrix(x)) {
    if (nrow(x) == Q) {
      x <- lapply(seq_len(ncol(x)), function(j) x[, j])
      names(x) <- colnames(x) %||% paste0("statistic", seq_along(x))
    } else if (ncol(x) == Q) {
      x <- lapply(seq_len(nrow(x)), function(i) x[i, ])
      names(x) <- rownames(x) %||% paste0("statistic", seq_along(x))
    }
  }
  if (!is.list(x) || !length(x)) {
    .dkge_abort("Rendered statistics must be a vector, matrix, or named list.",
                "dkge_renderer_error")
  }
  x <- lapply(x, as.numeric)
  if (any(vapply(x, length, integer(1)) != Q) ||
      any(!vapply(x, function(v) all(is.finite(v)), logical(1)))) {
    .dkge_abort("Every rendered statistic must be finite and match the support size.",
                "dkge_renderer_error")
  }
  if (is.null(names(x))) names(x) <- paste0("statistic", seq_along(x))
  x
}

#' Render aligned maps or an existing support-level statistic
#'
#' @param renderer A [dkge_renderer()] object.
#' @param x A [dkge_aligned_maps()] object or an already computed support-level
#'   vector, matrix, or named list.
#' @param statistic Aggregation used only when `x` contains aligned subject
#'   rows. `"weighted_mean"` uses the explicit weights recorded on `x`.
#' @param decode Logical; apply the support decoder when one exists.
#' @return A `dkge_rendered_alignment` with support values and optional decoded
#'   values. No correspondence is learned by this function.
#' @export
dkge_render_aligned <- function(renderer,
                                x,
                                statistic = c("weighted_mean", "mean", "median"),
                                decode = TRUE) {
  .dkge_validate_renderer(renderer)
  weighting <- "precomputed_statistic"
  estimand <- NULL
  alignment_hash <- NULL
  if (inherits(x, "dkge_aligned_maps")) {
    .dkge_validate_aligned_maps(x)
    if (!identical(x$support_id, renderer$support_id) ||
        !identical(x$support_hash, renderer$support_hash)) {
      .dkge_abort("Aligned maps and renderer identify different supports.",
                  "dkge_renderer_error")
    }
    statistic <- match.arg(statistic)
    values <- .dkge_renderer_statistic(x, statistic)
    weighting <- if (identical(statistic, "weighted_mean")) {
      x$subject_weighting
    } else {
      statistic
    }
    estimand <- x$estimand
    alignment_hash <- x$fitted_alignment_hash
  } else {
    values <- x
  }
  values <- .dkge_normalize_render_values(
    values, renderer$reference_support$n_locations
  )
  values <- lapply(values, function(v) {
    names(v) <- as.character(renderer$reference_support$labels)
    v
  })
  decoded <- NULL
  if (isTRUE(decode) && !is.null(renderer$reference_support$decoder)) {
    decoded <- lapply(values, function(v) {
      dkge_anchor_to_voxel_apply(renderer$reference_support$decoder, v)
    })
  }
  structure(
    list(
      support_id = renderer$support_id,
      support_hash = renderer$support_hash,
      support_values = values,
      decoded_values = decoded,
      statistic = if (inherits(x, "dkge_aligned_maps")) statistic else "supplied",
      subject_weighting = weighting,
      estimand = estimand,
      fitted_alignment_hash = alignment_hash,
      renderer_hash = renderer$structural_hash
    ),
    class = c("dkge_rendered_alignment", "list")
  )
}

#' @export
print.dkge_renderer <- function(x, ...) {
  cat("<dkge_renderer>\n")
  cat("  support  :", x$support_id, "\n")
  cat("  decoder  :", !is.null(x$reference_support$decoder), "\n")
  cat("  learns correspondence: FALSE\n")
  invisible(x)
}
