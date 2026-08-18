# dkge-diagnostics.R
# Compact, auditable summaries of fitted estimator contracts.

.dkge_moment_estimator_name <- function(fit) {
  debias <- fit$debias %||% "none"
  if (identical(debias, "none")) return("observed_second_moment")

  if (identical(debias, "analytic")) {
    dependencies <- vapply(fit$error_models %||% list(), function(model) {
      model$trial_dependence %||% "not_recorded"
    }, character(1))
    known <- dependencies[dependencies != "not_recorded"]
    if (length(known) && all(known == "iid")) {
      return("analytic_iid_error_correction")
    }
    if (any(known %in% c("covariance", "precision", "precision_operator"))) {
      return("analytic_covariance_aware_correction")
    }
    return("analytic_supplied_error_correction")
  }

  subjects <- fit$subjects %||% list()
  independent <- vapply(subjects, function(subject) {
    isTRUE(subject$split_provenance$independent)
  }, logical(1))
  if (length(independent) && all(independent)) {
    "independent_split_cross_moment"
  } else {
    "exploratory_split_cross_moment"
  }
}

.dkge_named_counts <- function(values, missing_label = "not_recorded") {
  values <- as.character(values)
  values[is.na(values) | !nzchar(values)] <- missing_label
  result <- table(values)
  stats::setNames(as.integer(result), names(result))
}

.dkge_error_model_summary <- function(fit) {
  models <- fit$error_models %||% list()
  if (!length(models)) {
    return(list(
      estimator = c(not_recorded = length(fit$subject_ids %||% list())),
      trial_dependence = c(not_recorded = length(fit$subject_ids %||% list())),
      effect_covariance_source = c(not_recorded = length(fit$subject_ids %||% list()))
    ))
  }
  field_counts <- function(field) {
    .dkge_named_counts(vapply(models, function(model) {
      if (is.null(model)) "not_recorded" else model[[field]] %||% "not_recorded"
    }, character(1)))
  }
  list(
    estimator = field_counts("estimator"),
    trial_dependence = field_counts("trial_dependence"),
    effect_covariance_source = field_counts("effect_covariance_source")
  )
}

.dkge_range_summary <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (!length(x)) return(c(min = NA_real_, median = NA_real_, max = NA_real_))
  c(min = min(x), median = stats::median(x), max = max(x))
}

#' Summarize a DKGE fit contract and diagnostics
#'
#' Reports effect coverage, pair effective sample size, negative spectral mass,
#' error-model provenance, the active moment estimator, and whether the fitted
#' object supports a physical block reconstruction or only q-space operations.
#'
#' @param object A fitted `dkge` object.
#' @param ... Unused.
#' @return An object of class `summary.dkge`.
#' @export
summary.dkge <- function(object, ...) {
  if (!inherits(object, "dkge")) {
    stop("`object` must inherit from 'dkge'.", call. = FALSE)
  }
  S <- length(object$subject_ids %||% object$Btil)
  q <- nrow(object$K)
  pair_counts <- object$pair_counts %||% matrix(S, q, q)
  pair_ess <- object$pair_ess %||% pair_counts
  effect_observed <- diag(pair_counts)
  effect_fraction <- if (S > 0L) effect_observed / S else rep(NA_real_, q)
  positive_ess <- pair_ess[pair_counts > 0 & is.finite(pair_ess)]
  spectral <- object$moment_diagnostics$transformed %||% list(
    negative_count = 0L,
    negative_mass = 0,
    total_spectral_mass = sum(abs(object$evals %||% 0))
  )
  total_mass <- spectral$total_spectral_mass %||% 0
  negative_mass <- spectral$negative_mass %||% 0

  result <- list(
    n_subjects = S,
    n_effects = q,
    rank = object$rank,
    rank_requested = object$rank_requested %||% object$rank,
    solver = object$solver %||% "pooled",
    representation = object$representation %||% "block_biprojector",
    representation_reasons = object$representation_reasons %||% character(0),
    moment_estimator = object$moment_estimator %||%
      .dkge_moment_estimator_name(object),
    effect_weighting = object$effect_weight_spec$method %||% "none",
    missingness = object$missingness %||% "none",
    effect_scaling = object$effect_scaling %||% "pooled_design",
    effect_observed_subjects = stats::setNames(
      as.numeric(effect_observed), object$effects %||% seq_len(q)
    ),
    effect_coverage_fraction = stats::setNames(
      as.numeric(effect_fraction), object$effects %||% seq_len(q)
    ),
    effect_coverage = .dkge_range_summary(effect_fraction),
    zero_coverage_pairs = sum(pair_counts == 0),
    total_pairs = length(pair_counts),
    pair_ess = .dkge_range_summary(positive_ess),
    negative_spectral_count = as.integer(spectral$negative_count %||% 0L),
    negative_spectral_mass = as.numeric(negative_mass),
    negative_spectral_fraction = if (is.finite(total_mass) && total_mass > 0) {
      as.numeric(negative_mass / total_mass)
    } else {
      0
    },
    error_models = .dkge_error_model_summary(object)
  )
  class(result) <- "summary.dkge"
  result
}

.dkge_format_counts <- function(x) {
  paste(paste(names(x), as.integer(x), sep = "="), collapse = ", ")
}

#' @export
print.summary.dkge <- function(x, ..., digits = max(3L, getOption("digits") - 3L)) {
  cat("Design-Kernel Group Embedding\n")
  cat(sprintf("  subjects: %d | effects: %d | rank: %d/%d\n",
              x$n_subjects, x$n_effects, x$rank, x$rank_requested))
  cat(sprintf("  estimator: %s | effect weights: %s | missingness: %s\n",
              x$moment_estimator, x$effect_weighting, x$missingness))
  cat(sprintf("  solver: %s | representation: %s | effect scaling: %s\n",
              x$solver, x$representation, x$effect_scaling))
  cat(sprintf(
    "  effect coverage: %.*f-%.*f (median %.*f) | zero-coverage pairs: %d/%d\n",
    digits, x$effect_coverage[["min"]], digits, x$effect_coverage[["max"]],
    digits, x$effect_coverage[["median"]], x$zero_coverage_pairs,
    x$total_pairs
  ))
  cat(sprintf("  pair ESS: %.*f-%.*f (median %.*f)\n",
              digits, x$pair_ess[["min"]], digits, x$pair_ess[["max"]],
              digits, x$pair_ess[["median"]]))
  cat(sprintf("  negative spectral mass: %.*g (%d eigenvalues; fraction %.*g)\n",
              digits, x$negative_spectral_mass, x$negative_spectral_count,
              digits, x$negative_spectral_fraction))
  cat(sprintf("  trial estimator: %s\n",
              .dkge_format_counts(x$error_models$estimator)))
  cat(sprintf("  trial dependence: %s\n",
              .dkge_format_counts(x$error_models$trial_dependence)))
  if (length(x$representation_reasons)) {
    cat("  representation reason: ",
        paste(x$representation_reasons, collapse = "; "), "\n", sep = "")
  }
  invisible(x)
}

#' Print a DKGE fit
#'
#' @param x A fitted `dkge` object.
#' @param ... Passed to the summary printer.
#' @return `x`, invisibly.
#' @export
print.dkge <- function(x, ...) {
  print(summary(x), ...)
  invisible(x)
}
