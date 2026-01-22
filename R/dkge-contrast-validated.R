# dkge-contrast-validated.R
# Dual-path (observed vs completed) DKGE contrasts with coverage diagnostics.

#' Dual-path DKGE contrasts with coverage diagnostics
#'
#' Computes DKGE contrasts twice on the same fold structure: once with an
#' observed-only coverage policy (typically `missingness = "rescale"`) and once with a
#' completed/penalised policy (e.g., `missingness = "shrink"`). Coverage metadata
#' from the fit and folds is returned together with simple sensitivity summaries.
#'
#' @param fit Fitted [dkge()] object.
#' @param contrasts Contrast specification accepted by [dkge_contrast()].
#' @param folds Fold definition (integer `k`, `dkge_folds`, data frame, or list).
#' @param ridge Ridge added to held-out compressed matrices.
#' @param parallel Logical; whether to parallelise held-out projections.
#' @param verbose Logical; emit progress messages.
#' @param align Logical; align held-out bases across folds.
#' @param observed_missingness Coverage policy for the observed-only path.
#' @param completed_missingness Coverage policy for the completed path.
#' @param observed_args Optional list of arguments for the observed policy
#'   (e.g., `list(min_pairs = 2)` for masking).
#' @param completed_args Optional list of arguments for the completed policy.
#' @param ... Additional arguments forwarded to the underlying cross-fitting
#'   helper.
#'
#' @return A list with class `dkge_contrast_validated` containing:
#'   - `observed`, `completed`: outputs from the respective paths.
#'   - `summary`: data frame with weighted means and sensitivity indices.
#'   - `provenance`: coverage metadata (effect IDs, subject masks, per-fold
#'     pair-count matrices).
#' @examples
#' \donttest{
#' toy <- dkge_sim_toy(
#'   factors = list(cond = list(L = 3)),
#'   active_terms = "cond", S = 4, P = 15, snr = 5
#' )
#' fit <- dkge(toy$B_list, toy$X_list, kernel = toy$K, rank = 2)
#' q <- nrow(fit$U)
#' c_vec <- rep(0, q)
#' c_vec[2] <- 1
#' c_vec[3] <- -1
#' res <- dkge_contrast_validated(fit,
#'                                contrasts = list(cond = c_vec),
#'                                folds = 3)
#' res$summary
#' }
#' @export
dkge_contrast_validated <- function(fit,
                                    contrasts,
                                    folds = NULL,
                                    ridge = 0,
                                    parallel = FALSE,
                                    verbose = FALSE,
                                    align = FALSE,
                                    observed_missingness = c("rescale", "mask", "none"),
                                    completed_missingness = c("shrink", "rescale", "none"),
                                    observed_args = list(),
                                    completed_args = list(),
                                    ...) {
  stopifnot(inherits(fit, "dkge"))
  observed_missingness <- match.arg(observed_missingness)
  completed_missingness <- match.arg(completed_missingness)

  contrast_list <- .normalize_contrasts(contrasts, fit)
  folds_obj <- .dkge_validated_coerce_folds(fit, folds)

  observed_res <- .dkge_contrast_kfold(
    fit,
    contrast_list = contrast_list,
    folds = folds_obj,
    ridge = ridge,
    parallel = parallel,
    verbose = verbose,
    align = align,
    missingness = observed_missingness,
    miss_args = observed_args,
    ...
  )

  completed_res <- .dkge_contrast_kfold(
    fit,
    contrast_list = contrast_list,
    folds = folds_obj,
    ridge = ridge,
    parallel = parallel,
    verbose = verbose,
    align = align,
    missingness = completed_missingness,
    miss_args = completed_args,
    ...
  )

  subject_ids <- fit$subject_ids %||% paste0("subject", seq_len(length(fit$Btil)))
  subject_weights <- fit$weights %||% rep(1, length(subject_ids))
  subject_weights <- as.numeric(subject_weights)
  sw <- sum(subject_weights)
  if (!is.finite(sw) || sw <= 0) {
    subject_weights <- rep(1 / length(subject_weights), length(subject_weights))
  } else {
    subject_weights <- subject_weights / sw
  }

  obs_scores <- .dkge_validated_subject_means(observed_res$values, subject_ids)
  comp_scores <- .dkge_validated_subject_means(completed_res$values, subject_ids)

  summary_tbl <- .dkge_validated_summary(obs_scores, comp_scores, subject_weights, names(contrast_list))

  provenance <- .dkge_validated_provenance(fit, folds_obj, observed_res, completed_res)

  out <- list(
    observed = observed_res,
    completed = completed_res,
    summary = summary_tbl,
    provenance = provenance
  )
  class(out) <- "dkge_contrast_validated"
  out
}

.dkge_validated_coerce_folds <- function(fit, folds) {
  if (is.null(folds)) {
    stop("folds must be provided for dkge_contrast_validated().", call. = FALSE)
  }
  if (is.numeric(folds) && length(folds) == 1) {
    return(dkge_define_folds(fit, type = "subject", k = folds))
  }
  if (inherits(folds, "dkge_folds")) {
    return(folds)
  }
  folds_obj <- as_dkge_folds(folds, fit_or_data = fit)
  if (!inherits(folds_obj, "dkge_folds")) {
    stop("Unable to coerce `folds` into a dkge_folds object.", call. = FALSE)
  }
  folds_obj
}

.dkge_validated_subject_means <- function(values, subject_ids) {
  contrasts <- names(values)
  res <- matrix(NA_real_, nrow = length(subject_ids), ncol = length(contrasts),
                dimnames = list(subject_ids, contrasts))
  for (contrast_name in contrasts) {
    subj_vals <- values[[contrast_name]]
    for (id in names(subj_vals)) {
      if (!id %in% subject_ids) next
      v <- subj_vals[[id]]
      if (is.null(v)) next
      res[id, contrast_name] <- mean(as.numeric(v), na.rm = TRUE)
    }
  }
  res
}

.dkge_validated_weighted_mean <- function(x, w) {
  keep <- !is.na(x) & !is.na(w) & w > 0
  if (!any(keep)) return(NA_real_)
  w_sub <- w[keep]
  w_sub <- w_sub / sum(w_sub)
  sum(x[keep] * w_sub)
}

.dkge_validated_summary <- function(obs_scores, comp_scores, weights, contrast_names) {
  w_vec <- weights
  names(w_vec) <- rownames(obs_scores)
  summary_rows <- lapply(seq_along(contrast_names), function(i) {
    cname <- contrast_names[[i]]
    obs <- obs_scores[, cname, drop = TRUE]
    comp <- comp_scores[, cname, drop = TRUE]
    mu_obs <- .dkge_validated_weighted_mean(obs, w_vec)
    mu_comp <- .dkge_validated_weighted_mean(comp, w_vec)
    delta <- mu_comp - mu_obs
    sens <- abs(delta) / (abs(mu_obs) + 1e-12)
    data.frame(
      contrast = cname,
      estimate_observed = mu_obs,
      estimate_completed = mu_comp,
      delta = delta,
      sensitivity = sens,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, summary_rows)
}

.dkge_validated_provenance <- function(fit, folds_obj, obs_res, comp_res) {
  assignments <- folds_obj$assignments
  fold_meta <- vector("list", length(assignments))
  effect_ids <- fit$provenance$effect_ids %||% rownames(fit$K)
  obs_mask <- fit$provenance$obs_mask %||% NULL
  for (i in seq_along(assignments)) {
    test_idx <- assignments[[i]]
    train_idx <- setdiff(seq_len(length(fit$Btil)), test_idx)
    fold_meta[[i]] <- list(
      fold = i,
      train_idx = train_idx,
      test_idx = test_idx,
      pair_counts_observed = obs_res$metadata$pair_counts[[i]],
      pair_counts_completed = comp_res$metadata$pair_counts[[i]]
    )
  }
  list(
    effect_ids = effect_ids,
    obs_mask = obs_mask,
    folds = fold_meta
  )
}

#' @export
print.dkge_contrast_validated <- function(x, ...) {
  cat("DKGE Validated Contrasts\n")
  cat("-------------------------\n")
  print(x$summary)
  invisible(x)
}
