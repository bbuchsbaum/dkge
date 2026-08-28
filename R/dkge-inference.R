
# dkge-inference.R
# Unified inference module for DKGE contrasts with multiple testing procedures

.dkge_validate_permutation_count <- function(x, argument = "B",
                                              minimum = 100L) {
  valid <- is.numeric(x) && is.null(dim(x)) && length(x) == 1L &&
    !is.na(x) && is.finite(x) && x == floor(x) &&
    x >= minimum && x <= .Machine$integer.max
  if (!valid) {
    .dkge_abort(
      sprintf(
        "`%s` must be one finite whole number between %d and %d.",
        argument, minimum, .Machine$integer.max
      ),
      "dkge_inference_permutation_error"
    )
  }
  as.integer(x)
}

.dkge_validate_signflip_matrix <- function(Y) {
  Y <- as.matrix(Y)
  if (!is.numeric(Y) || length(dim(Y)) != 2L || nrow(Y) < 5L ||
      ncol(Y) < 1L || any(!is.finite(Y))) {
    .dkge_abort(
      paste0(
        "Sign-flip data must be a finite numeric matrix with at least five ",
        "subjects and one feature."
      ),
      "dkge_inference_data_error"
    )
  }
  subject_ids <- rownames(Y)
  if (!is.null(subject_ids) &&
      (length(subject_ids) != nrow(Y) || anyNA(subject_ids) ||
       any(!nzchar(subject_ids)) || anyDuplicated(subject_ids))) {
    .dkge_abort("Sign-flip subject identifiers must be unique and non-empty.",
                "dkge_inference_subject_error")
  }
  feature_ids <- colnames(Y)
  if (!is.null(feature_ids) &&
      (length(feature_ids) != ncol(Y) || anyNA(feature_ids) ||
       any(!nzchar(feature_ids)) || anyDuplicated(feature_ids))) {
    .dkge_abort("Sign-flip feature identifiers must be unique and non-empty.",
                "dkge_inference_feature_error")
  }
  Y
}

#' One-sample sign-flip max-T inference on transported subject maps
#'
#' Computes cluster-wise one-sample t-statistics across subjects on transported
#' values (SxQ matrix), and calibrates p-values by the max-|t| distribution under
#' random subject-wise sign flips (symmetric null). The helper conditions on the
#' supplied matrix and is exact only when that complete matrix is jointly
#' row-sign invariant. It does not establish that an upstream adaptive DKGE
#' estimator has this property.
#'
#' @param Y SxQ matrix of aligned subject values on one identified reference
#'   support (rows = subjects, columns = support locations).
#' @param B number of sign-flip permutations
#' @param center "mean" or "median" for the location statistic (t uses mean)
#' @param tail "two.sided" | "greater" | "less"
#' @return A list with fields: `stat` (Q-vector of observed t-statistics), `p`
#'   (Q-vector of max-T family-wise-error adjusted p-values), `p_unadj`
#'   (Q-vector of per-column unadjusted permutation p-values), `maxnull`
#'   (B-vector of permutation maximum statistics), and `flips` (S-by-B sign
#'   matrix). Statistic and p-value names follow `colnames(Y)` (or stable
#'   `feature*` defaults); flip rows follow `rownames(Y)` (or `subject*`
#'   defaults).
#' @export
dkge_signflip_maxT <- function(Y, B = 2000, center = c("mean","median"),
                               tail = c("two.sided","greater","less")) {
  center <- match.arg(center); tail <- match.arg(tail)
  if (!identical(center, "mean")) {
    warning(
      "Non-mean `center` modes are deprecated because max-T uses a studentized mean statistic.",
      call. = FALSE
    )
    .dkge_abort(
      "`dkge_signflip_maxT()` supports only `center = 'mean'`.",
      "dkge_inference_center_error"
    )
  }
  B <- .dkge_validate_permutation_count(B, "B", minimum = 100L)
  Y <- .dkge_validate_signflip_matrix(Y)
  S <- nrow(Y); Q <- ncol(Y)
  subject_ids <- rownames(Y) %||% paste0("subject", seq_len(S))
  if (length(subject_ids) != S || anyNA(subject_ids) ||
      any(!nzchar(subject_ids)) || anyDuplicated(subject_ids)) {
    .dkge_abort("Sign-flip subject identifiers must be unique and non-empty.",
                "dkge_inference_subject_error")
  }
  feature_ids <- colnames(Y) %||% paste0("feature", seq_len(Q))
  permutation_ids <- paste0("perm", seq_len(B))
  canonical_order <- order(subject_ids, method = "radix")
  canonical_ids <- subject_ids[canonical_order]
  Y_work <- Y[canonical_order, , drop = FALSE]

  # observed t per cluster
  mu  <- colMeans(Y_work)
  sdv <- apply(Y_work, 2, stats::sd)
  t_obs <- mu / (sdv / sqrt(S) + 1e-12)

  # generate random sign matrix (SxB)
  flips_work <- matrix(sample(c(-1,1), S*B, replace=TRUE), S, B)
  dimnames(flips_work) <- list(canonical_ids, permutation_ids)
  Yc <- Y_work

  # Observed statistic on the tested side.
  obs_side <- switch(tail,
                     two.sided = abs(t_obs),
                     greater   = t_obs,
                     less      = -t_obs)

  # permutation max stat, matched to the tested tail so one-sided tests use a
  # signed (not absolute) null; using max(abs(.)) for a one-sided test inflates
  # the null and needlessly costs power. Accumulate per-cluster exceedances too
  # so an *uncorrected* p-value is available alongside the max-T adjusted one.
  maxnull <- numeric(B)
  exceed_unadj <- numeric(length(t_obs))
  for (b in seq_len(B)) {
    Yb <- flips_work[,b] * Yc
    mu_b  <- colMeans(Yb)
    sd_b  <- apply(Yb, 2, stats::sd)
    t_b   <- mu_b / (sd_b / sqrt(S) + 1e-12)
    stat_b <- switch(tail,
                     two.sided = abs(t_b),
                     greater   = t_b,
                     less      = -t_b)
    maxnull[b] <- max(stat_b)
    exceed_unadj <- exceed_unadj + (stat_b >= obs_side)
  }

  # p-values: max-T (strong FWER control) and per-cluster uncorrected.
  p <- sapply(obs_side, function(x) (1 + sum(maxnull >= x)) / (B + 1))
  p_unadj <- (1 + exceed_unadj) / (B + 1)
  names(t_obs) <- feature_ids
  names(p) <- feature_ids
  names(p_unadj) <- feature_ids
  names(maxnull) <- permutation_ids
  flips <- flips_work[match(subject_ids, canonical_ids), , drop = FALSE]
  rownames(flips) <- subject_ids

  list(stat = t_obs, p = p, p_unadj = p_unadj, maxnull = maxnull, flips = flips)
}

.dkge_rank_truncated_estimator_eligibility <- function(method) {
  method <- match.arg(method, c("loso", "kfold", "analytic"))
  eligibility <- dkge_alignment_eligibility(
    feature_source = "geometry_only",
    estimator_source = "same_data_rank_truncated",
    source_verified = TRUE
  )
  estimator_detail <- if (identical(method, "analytic")) {
    paste0(
      "; analytic LOSO uses a first-order basis approximation unless its ",
      "diagnostic selects an exact-LOSO fallback, and neither path makes the ",
      "cohort-estimated latent span invariant to joint null actions"
    )
  } else {
    paste0(
      "; cross-fitting omits the tested subject from its own basis but does ",
      "not make the cohort-estimated latent span invariant to joint null ",
      "actions"
    )
  }
  eligibility$reason <- paste0(eligibility$reason, estimator_detail)
  eligibility
}

.dkge_validate_inference_correction <- function(inference, correction) {
  if (identical(inference, "parametric") && identical(correction, "maxT")) {
    .dkge_abort(
      paste0(
        "Parametric inference does not define a max-T null distribution. ",
        "Choose `correction = 'fdr'`, `'bonferroni'`, or `'none'`, or use ",
        "`inference = 'signflip'` for max-T."
      ),
      "dkge_inference_correction_error"
    )
  }
  invisible(NULL)
}

#' Unified inference for DKGE contrasts
#'
#' High-level interface for statistical inference on DKGE contrasts with
#' integrated cross-fitting and multiple testing correction.
#'
#' @param fit A `dkge` object from [dkge_fit()] or [dkge()]
#' @param contrasts Contrast specification (see [dkge_contrast()])
#' @param method Cross-fitting method: "loso", "kfold", or "analytic"
#' @param inference Inference type:
#'   - `"signflip"`: Sign-flip permutation test (default)
#'   - `"freedman-lane"`: Freedman-Lane permutation (requires adapters)
#'   - `"parametric"`: Parametric t-test (assumes normality)
#' @param correction Multiple testing correction:
#'   - `"maxT"`: Family-wise error rate via max-T (default)
#'   - `"fdr"`: False discovery rate (Benjamini-Hochberg)
#'   - `"bonferroni"`: Bonferroni correction
#'   - `"none"`: No correction
#' @param n_perm Number of permutations for non-parametric tests
#' @param alpha Significance level for corrections
#' @param transported Deprecated. Must be `FALSE`. Rendering/alignment is no
#'   longer learned inside this inference helper.
#' @param transport Deprecated. Must be `NULL`. Use the typed two-stage workflow
#'   [dkge_transport_contrasts_to_reference()] then [dkge_infer_aligned()].
#' @param allow_approximate_alignment Logical; permit an explicitly labelled
#'   `"approximate"` alignment or same-data rank-truncated LOSO, K-fold, or
#'   analytic estimator. Ineligible/descriptive states are always refused. The
#'   default is fail-closed.
#' @param ... Additional arguments passed to [dkge_contrast()] and inference functions
#'
#' @return An object of class `dkge_inference` containing:
#'   - `contrasts`: The contrast results from cross-fitting
#'   - `statistics`: Test statistics per cluster/voxel
#'   - `p_values`: Raw p-values
#'   - `p_adjusted`: Adjusted p-values based on correction method
#'   - `significant`: Logical indicators of significance
#'   - `method`: Cross-fitting method used
#'   - `inference`: Inference type used
#'   - `correction`: Correction method applied
#'   - `metadata`: Additional information about the analysis
#'
#' @details
#' This function integrates the cross-fitting machinery from [dkge_contrast()]
#' with various statistical inference procedures. It first computes contrast
#' values using the specified cross-fitting method, then applies the chosen
#' inference procedure to obtain p-values, and finally applies multiple
#' testing correction.
#'
#' This helper performs inference only in the contrast result's native support.
#' It never learns correspondence. When subject supports differ, first build
#' typed alignment features, transport with
#' [dkge_transport_contrasts_to_reference()], extract its `dkge_aligned_maps`
#' object, and pass that object to [dkge_infer_aligned()].
#'
#' The workflow is:
#' 1. Compute contrast values via cross-fitting (LOSO/K-fold/analytic)
#' 2. Apply inference procedure (sign-flip/Freedman-Lane/parametric)
#' 3. Apply multiple testing correction (maxT/FDR/Bonferroni)
#'
#' Max-T uses one subject-sign action and one maximum over the entire
#' prespecified contrast-by-location family. Its conditional FWER guarantee
#' requires joint row-sign invariance of the supplied subject matrix. Standard
#' rank-truncated LOSO/K-fold DKGE estimates retain same-data latent-span
#' dependence and are therefore labelled approximate unless the full estimator
#' is rebuilt under every null action. The analytic method is also approximate:
#' it uses a first-order approximation to the LOSO basis except where its
#' diagnostic selects an exact-LOSO fallback.
#'
#' @examples
#' # Simulate and fit
#' toy <- dkge_sim_toy(
#'   factors = list(A = list(L = 2), B = list(L = 3)),
#'   active_terms = c("A", "B"), S = 6, P = 15, snr = 5
#' )
#' fit <- dkge(toy$B_list, toy$X_list, kernel = toy$K, rank = 2)
#'
#' # LOSO with sign-flip and maxT correction (fast with few perms for example)
#' \donttest{
#' results <- dkge_infer(
#'   fit, c(1, rep(0, 4)), n_perm = 100,
#'   allow_approximate_alignment = TRUE
#' )
#' results
#' }
#'
#' @seealso [dkge_contrast()], [dkge_infer_aligned()],
#'   [dkge_transport_contrasts_to_reference()], [dkge_signflip_maxT()],
#'   [dkge_freedman_lane()]
#' @export
dkge_infer <- function(fit, contrasts,
                      method = c("loso", "kfold", "analytic"),
                      inference = c("signflip", "freedman-lane", "parametric"),
                      correction = c("maxT", "fdr", "bonferroni", "none"),
                      n_perm = 2000,
                      alpha = 0.05,
                      transported = FALSE,
                      transport = NULL,
                      allow_approximate_alignment = FALSE,
                      ...) {
  method <- match.arg(method)
  inference <- match.arg(inference)
  correction <- match.arg(correction)
  .dkge_validate_inference_correction(inference, correction)
  alpha <- .dkge_validate_probability(alpha, "alpha")
  if (identical(inference, "freedman-lane")) {
    .dkge_abort(
      paste0(
        "Freedman-Lane inference is not implemented in `dkge_infer()`; ",
        "use `dkge_freedman_lane()` with the required time-series adapters."
      ),
      "dkge_inference_compatibility_error"
    )
  }
  if (identical(inference, "signflip")) {
    n_perm <- .dkge_validate_permutation_count(
      n_perm, "n_perm", minimum = 100L
    )
  }

  if (isTRUE(transported) || !is.null(transport)) {
    .dkge_abort(
      paste0(
        "Transport inside `dkge_infer()` is retired because that facade cannot ",
        "establish typed functional-correspondence provenance. Use ",
        "`dkge_transport_contrasts_to_reference()` followed by ",
        "`dkge_infer_aligned()`."
      ),
      "dkge_alignment_ineligible_error"
    )
  }

  contrast_results <- dkge_contrast(fit, contrasts, method = method, ...)

  mapped_values <- NULL
  estimator_eligibility <- .dkge_rank_truncated_estimator_eligibility(method)
  .dkge_assert_alignment_eligible(
    estimator_eligibility,
    allow_approximate = allow_approximate_alignment
  )

  # Step 2: Apply inference procedure
  infer_results <- switch(inference,
    signflip = .infer_signflip(contrast_results, n_perm, correction, mapped_values),
    `freedman-lane` = .infer_freedman_lane(contrast_results, n_perm, correction, ...),
    parametric = .infer_parametric(contrast_results, correction, mapped_values)
  )

  # Step 3: Apply correction (if not already done by inference method)
  if (inference != "signflip" || correction != "maxT") {
    infer_results <- .apply_correction(infer_results, correction, alpha)
  }

  # Mark significant results
  infer_results$significant <- lapply(infer_results$p_adjusted, function(p) p <= alpha)
  infer_results$alpha <- alpha

  infer_results$metadata$estimator <- list(
    status = estimator_eligibility$status,
    reason = estimator_eligibility$reason,
    method = method,
    approximate_override = isTRUE(allow_approximate_alignment)
  )

  structure(infer_results, class = "dkge_inference")
}

#' Group inference on already aligned subject maps
#'
#' This is the typed inference boundary for functional alignment. It consumes
#' one row per subject on one identified support, never fits correspondence,
#' and refuses ineligible alignment objects by default. Fit-level MFA pooling
#' weights are not inherited; the aligned object records its group weighting.
#'
#' @param aligned_maps An operator-bound [dkge_aligned_maps()] object returned by
#'   [dkge_transport_contrasts_to_reference()] or [dkge_align_to_template()].
#' @param inference `"signflip"` or `"parametric"`.
#' @param correction Multiple-testing correction.
#' @param n_perm Number of sign flips.
#' @param alpha Significance level.
#' @param allow_approximate_alignment Permit an explicitly labelled approximate
#'   alignment. Descriptive/ineligible objects remain errors.
#' @details
#' The subject is the sampling and exchangeability unit. Locations and
#' contrasts within a subject row are not independent replicates. Sign-flip
#' inference requires joint row-sign invariance under the null: one sign is
#' applied to each subject simultaneously across every prespecified
#' contrast-by-location cell. With `correction = "maxT"`, each null draw takes
#' the maximum statistic over that complete joint family. Marginal permutation
#' p-values use the same shared sign matrix, and FDR or Bonferroni correction is
#' applied once after flattening that complete family (then split back into the
#' original contrast shapes). `correction = "none"` leaves those joint-draw
#' marginal p-values unadjusted.
#'
#' Correspondence is fixed before this boundary. Independent-feature and
#' same-data-residualized alignments can still depend on a same-cohort,
#' rank-truncated latent span. The frozen calibration court's full-pipeline
#' re-estimation exact null-action comparator passed, but the negative-control
#' promotion gate did not; these frozen-plan modes therefore remain
#' `"approximate"`, fail closed by default, and require an explicit override.
#'
#' Parametric inference performs cellwise one-sample t tests and assumes
#' independent subjects plus an adequate normal approximation; it does not
#' define a max-T randomization distribution. Current aligned group inference
#' supports equal subject weighting only. Fit-level MFA weights and explicit
#' descriptive/bootstrap weights are not imported into this group test.
#' @return A `dkge_inference` object carrying the aligned estimand and weighting.
#' @export
dkge_infer_aligned <- function(
    aligned_maps,
    inference = c("signflip", "parametric"),
    correction = c("maxT", "fdr", "bonferroni", "none"),
    n_perm = 2000,
    alpha = 0.05,
    allow_approximate_alignment = FALSE) {
  .dkge_validate_aligned_maps(aligned_maps)
  .dkge_assert_alignment_eligible(
    aligned_maps,
    allow_approximate = allow_approximate_alignment
  )
  if (!identical(aligned_maps$subject_weighting, "equal_subject")) {
    .dkge_abort(
      paste0(
        "Aligned-map inference currently supports only equal subject weighting; ",
        "the object records '", aligned_maps$subject_weighting, "'."
      ),
      "dkge_alignment_weighting_error"
    )
  }
  inference <- match.arg(inference)
  correction <- match.arg(correction)
  .dkge_validate_inference_correction(inference, correction)
  alpha <- .dkge_validate_probability(alpha, "alpha")
  if (identical(inference, "signflip")) {
    n_perm <- .dkge_validate_permutation_count(
      n_perm, "n_perm", minimum = 100L
    )
  }
  contrast_ids <- as.character(aligned_maps$contrast_ids)
  contrast_stub <- list(
    contrasts = stats::setNames(contrast_ids, contrast_ids),
    values = aligned_maps$values,
    method = "aligned_reference_support"
  )
  infer_results <- switch(
    inference,
    signflip = .infer_signflip(
      contrast_stub, n_perm, correction, aligned_maps$values
    ),
    parametric = .infer_parametric(
      contrast_stub, correction, aligned_maps$values
    )
  )
  if (inference != "signflip" || correction != "maxT") {
    infer_results <- .apply_correction(infer_results, correction, alpha)
  }
  infer_results$significant <- lapply(
    infer_results$p_adjusted, function(p) p <= alpha
  )
  infer_results$alpha <- alpha
  infer_results$aligned_maps <- aligned_maps
  infer_results$metadata$alignment <- list(
    status = aligned_maps$eligibility$status,
    feature_source = aligned_maps$feature_source,
    reason = aligned_maps$eligibility$reason,
    estimand = aligned_maps$estimand,
    subject_weighting = aligned_maps$subject_weighting,
    approximate_override = isTRUE(allow_approximate_alignment)
  )
  structure(infer_results, class = "dkge_inference")
}

#' Sign-flip inference helper
#' @keywords internal
#' @noRd
.infer_signflip <- function(contrast_results, n_perm, correction,
                            mapped_values = NULL,
                            tail = "two.sided", center = "mean") {
  n_perm <- .dkge_validate_permutation_count(
    n_perm, "n_perm", minimum = 100L
  )
  n_contrasts <- length(contrast_results$contrasts)

  stats <- vector("list", n_contrasts)
  p_values <- vector("list", n_contrasts)
  p_adjusted <- vector("list", n_contrasts)

  Ys <- lapply(seq_len(n_contrasts), function(i) {
    if (!is.null(mapped_values)) mapped_values[[i]] else
      as.matrix(contrast_results, contrast = i)
  })
  Ys <- lapply(Ys, .dkge_validate_signflip_matrix)
  if (!length(Ys) || any(vapply(Ys, nrow, integer(1)) != nrow(Ys[[1L]]))) {
    .dkge_abort("Every contrast family member must contain the same subjects.",
                "dkge_inference_subject_error")
  }
  row_ids <- lapply(Ys, rownames)
  named_rows <- !vapply(row_ids, is.null, logical(1))
  if (any(named_rows) && !all(named_rows)) {
    .dkge_abort("Contrast matrices mix named and unnamed subject rows.",
                "dkge_inference_subject_error")
  }
  if (all(named_rows)) {
    target_ids <- row_ids[[1L]]
    if (anyNA(target_ids) || any(!nzchar(target_ids)) ||
        anyDuplicated(target_ids) || any(vapply(row_ids, function(ids) {
          anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids) ||
            !setequal(ids, target_ids)
        }, logical(1)))) {
      .dkge_abort("Contrast matrices do not identify the same subject cohort.",
                  "dkge_inference_subject_error")
    }
    target_ids <- sort(target_ids, method = "radix")
    Ys <- lapply(Ys, function(Y) Y[target_ids, , drop = FALSE])
  }
  contrast_ids <- names(mapped_values) %||% names(contrast_results$values) %||%
    names(contrast_results$contrasts)
  if (is.null(contrast_ids) && is.character(contrast_results$contrasts) &&
      length(contrast_results$contrasts) == n_contrasts) {
    contrast_ids <- contrast_results$contrasts
  }
  contrast_ids <- as.character(contrast_ids %||%
                                 paste0("contrast", seq_len(n_contrasts)))
  names(stats) <- contrast_ids
  names(p_values) <- contrast_ids
  names(p_adjusted) <- contrast_ids
  dimensions <- vapply(Ys, ncol, integer(1))
  names(dimensions) <- contrast_ids

  if (correction == "maxT") {
    feature_ids <- Map(function(Y, id) {
      paste0(id, "::", colnames(Y) %||% paste0("location", seq_len(ncol(Y))))
    }, Ys, contrast_ids)
    family <- do.call(cbind, unname(Ys))
    colnames(family) <- unlist(feature_ids, use.names = FALSE)
    result <- dkge_signflip_maxT(
      family, B = n_perm, tail = tail, center = center
    )
    ends <- cumsum(dimensions)
    starts <- c(1L, head(ends, -1L) + 1L)
    for (i in seq_len(n_contrasts)) {
      idx <- starts[[i]]:ends[[i]]
      local_names <- colnames(Ys[[i]]) %||%
        paste0("location", seq_len(dimensions[[i]]))
      stats[[i]] <- stats::setNames(unname(result$stat[idx]), local_names)
      p_values[[i]] <- stats::setNames(
        unname((result$p_unadj %||% result$p)[idx]), local_names
      )
      p_adjusted[[i]] <- stats::setNames(unname(result$p[idx]), local_names)
    }
    return(list(
      contrasts = contrast_results,
      statistics = stats,
      p_values = p_values,
      p_adjusted = p_adjusted,
      method = contrast_results$method,
      inference = "signflip",
      correction = correction,
      metadata = list(
        n_perm = n_perm,
        tail = tail,
        center = center,
        family_scope = "all_contrasts_and_support_locations",
        shared_subject_signs = TRUE,
        family_dimensions = dimensions,
        adjustment_scope = "all_contrasts_and_support_locations",
        adjustment_method = "maxT",
        n_family_tests = sum(dimensions),
        sign_flips = result$flips,
        max_null = result$maxnull
      )
    ))
  }

  n_subjects <- nrow(Ys[[1L]])
  flips <- matrix(
    sample(c(-1, 1), n_subjects * n_perm, replace = TRUE),
    n_subjects, n_perm
  )
  dimnames(flips) <- list(
    rownames(Ys[[1L]]) %||% paste0("subject", seq_len(n_subjects)),
    paste0("perm", seq_len(n_perm))
  )

  for (i in seq_len(n_contrasts)) {
    Y <- Ys[[i]]

    mu <- colMeans(Y)
    se <- apply(Y, 2, stats::sd) / sqrt(pmax(n_subjects, 1))
    stats[[i]] <- mu / (se + 1e-12)
    local_names <- colnames(Y) %||%
      paste0("location", seq_len(ncol(Y)))
    names(stats[[i]]) <- local_names
    null_stats <- matrix(vapply(seq_len(n_perm), function(b) {
      s <- flips[, b]
      Y_flip <- s * Y
      mu_flip <- colMeans(Y_flip)
      se_flip <- apply(Y_flip, 2, stats::sd) / sqrt(pmax(n_subjects, 1))
      mu_flip / (se_flip + 1e-12)
    }, numeric(ncol(Y))), nrow = ncol(Y), ncol = n_perm)
    observed_side <- switch(
      tail,
      two.sided = abs(stats[[i]]),
      greater = stats[[i]],
      less = -stats[[i]]
    )
    null_side <- switch(
      tail,
      two.sided = abs(null_stats),
      greater = null_stats,
      less = -null_stats
    )

    p_values[[i]] <- vapply(seq_along(observed_side), function(j) {
      (1 + sum(null_side[j, ] >= observed_side[[j]])) / (n_perm + 1)
    }, numeric(1))
    names(p_values[[i]]) <- local_names

    p_adjusted[[i]] <- p_values[[i]]
  }

  list(
    contrasts = contrast_results,
    statistics = stats,
    p_values = p_values,
    p_adjusted = p_adjusted,
    method = contrast_results$method,
    inference = "signflip",
    correction = correction,
    metadata = list(
      n_perm = n_perm,
      tail = tail,
      center = center,
      family_scope = "all_contrasts_and_support_locations",
      shared_subject_signs = TRUE,
      family_dimensions = dimensions,
      adjustment_scope = "none",
      adjustment_method = "none",
      n_family_tests = sum(dimensions),
      sign_flips = flips
    )
  )
}

#' Parametric inference helper
#' @keywords internal
#' @noRd
.infer_parametric <- function(contrast_results, correction, mapped_values = NULL) {
  n_contrasts <- length(contrast_results$contrasts)

  stats <- vector("list", n_contrasts)
  p_values <- vector("list", n_contrasts)
  df_vec <- numeric(n_contrasts)

  for (i in seq_len(n_contrasts)) {
    Y <- if (!is.null(mapped_values)) {
      mapped_values[[i]]
    } else {
      as.matrix(contrast_results, contrast = i)
    }
    n_subjects <- nrow(Y)
    df_vec[i] <- n_subjects - 1
    mu <- colMeans(Y)
    se <- apply(Y, 2, stats::sd) / sqrt(n_subjects)
    t_stats <- mu / (se + 1e-12)

    stats[[i]] <- t_stats
    p_values[[i]] <- 2 * pt(-abs(t_stats), df_vec[i])
  }

  dimensions <- vapply(p_values, length, integer(1))
  contrast_ids <- names(mapped_values) %||% names(contrast_results$values) %||%
    names(contrast_results$contrasts)
  if (is.null(contrast_ids) && is.character(contrast_results$contrasts) &&
      length(contrast_results$contrasts) == n_contrasts) {
    contrast_ids <- contrast_results$contrasts
  }
  contrast_ids <- as.character(
    contrast_ids %||% paste0("contrast", seq_len(n_contrasts))
  )
  names(stats) <- contrast_ids
  names(p_values) <- contrast_ids
  names(df_vec) <- contrast_ids
  names(dimensions) <- contrast_ids

  list(
    contrasts = contrast_results,
    statistics = stats,
    p_values = p_values,
    p_adjusted = p_values,  # Will be corrected next
    method = contrast_results$method,
    inference = "parametric",
    correction = correction,
    metadata = list(
      df = df_vec,
      family_scope = "all_contrasts_and_support_locations",
      family_dimensions = dimensions,
      adjustment_scope = "none",
      adjustment_method = "none",
      n_family_tests = sum(dimensions)
    )
  )
}

#' Apply multiple testing correction
#' @keywords internal
#' @noRd
.apply_correction <- function(infer_results, correction, alpha) {
  if (correction == "none") {
    infer_results$metadata$adjustment_scope <- "none"
    infer_results$metadata$adjustment_method <- "none"
    return(infer_results)
  }
  if (identical(correction, "maxT")) {
    .dkge_abort(
      "Max-T correction must be computed from a joint sign-flip family.",
      "dkge_inference_correction_error"
    )
  }

  dimensions <- lengths(infer_results$p_values)
  family_p <- unlist(
    lapply(infer_results$p_values, as.numeric), use.names = FALSE
  )
  family_adjusted <- stats::p.adjust(family_p, method = correction)
  cursor <- 0L
  for (i in seq_along(dimensions)) {
    size <- dimensions[[i]]
    idx <- cursor + seq_len(size)
    infer_results$p_adjusted[[i]] <- stats::setNames(
      family_adjusted[idx], names(infer_results$p_values[[i]])
    )
    cursor <- cursor + size
  }
  infer_results$metadata$family_scope <-
    "all_contrasts_and_support_locations"
  infer_results$metadata$adjustment_scope <-
    "all_contrasts_and_support_locations"
  infer_results$metadata$adjustment_method <- correction
  infer_results$metadata$n_family_tests <- length(family_p)

  infer_results
}

#' Freedman-Lane inference helper (placeholder)
#' @keywords internal
#' @noRd
.infer_freedman_lane <- function(contrast_results, n_perm, correction, ...) {
  stop("Freedman-Lane inference requires time-series data and GLM adapters. ",
       "See dkge_freedman_lane() for the scaffold implementation.")
}

#' Print method for dkge_inference
#'
#' @param x A dkge_inference object
#' @param ... Additional arguments (unused)
#' @export
print.dkge_inference <- function(x, ...) {
  cat("DKGE Inference Results\n")
  cat("----------------------\n")
  cat(sprintf("Cross-fitting: %s\n", x$method))
  cat(sprintf("Inference: %s\n", x$inference))
  cat(sprintf("Correction: %s\n", x$correction))
  if (!is.null(x$metadata$alignment)) {
    cat(sprintf("Alignment: %s (%s)\n",
                x$metadata$alignment$status,
                x$metadata$alignment$feature_source))
    cat(sprintf("Group weighting: %s\n",
                x$metadata$alignment$subject_weighting))
  }

  n_contrasts <- length(x$statistics)
  cat(sprintf("Contrasts: %d\n", n_contrasts))

  if (!is.null(x$alpha)) {
    cat(sprintf("Alpha level: %g\n", x$alpha))
  }

  for (i in seq_len(min(5, n_contrasts))) {
    n_sig <- sum(x$significant[[i]])
    n_total <- length(x$significant[[i]])
    cat(sprintf("  %s: %d/%d significant\n",
               names(x$contrasts$contrasts)[i], n_sig, n_total))
  }

  if (n_contrasts > 5) {
	cat("  ...\n")
  }

  invisible(x)
}

#' Convert DKGE inference results to a tidy data frame
#'
#' @param x A `dkge_inference` object
#' @param row.names NULL or a character vector giving the row names
#' @param optional Logical; if TRUE, setting row names is optional
#' @param ... Additional arguments passed to [base::data.frame()], including
#'   `stringsAsFactors`
#' @return Data frame with columns `contrast`, `cluster`, `statistic`, `p_value`,
#'   `p_adjusted`, and `significant`
#' @export
as.data.frame.dkge_inference <- function(x, row.names = NULL, optional = FALSE, ...) {
  dots <- list(...)
  stringsAsFactors <- dots$stringsAsFactors %||% FALSE
  dots$stringsAsFactors <- NULL

  contrast_names <- names(x$statistics)
  if (is.null(contrast_names) || any(!nzchar(contrast_names))) {
    contrast_names <- names(x$contrasts$contrasts)
  }
  if (is.null(contrast_names) || length(contrast_names) != length(x$statistics)) {
    contrast_names <- paste0("contrast", seq_along(x$statistics))
  }

  alpha_val <- x$alpha %||% NA_real_
  method_val <- x$method %||% NA_character_
  inference_val <- x$inference %||% NA_character_
  correction_val <- x$correction %||% NA_character_

  rows <- vector("list", length(x$statistics))
  for (i in seq_along(x$statistics)) {
    stats <- x$statistics[[i]]
    if (is.null(stats)) {
      next
    }
    cluster_ids <- names(stats)
    if (is.null(cluster_ids) || any(!nzchar(cluster_ids))) {
      cluster_ids <- paste0("cluster", seq_along(stats))
    }
    p_vals <- x$p_values[[i]] %||% rep(NA_real_, length(stats))
    padj <- x$p_adjusted[[i]] %||% p_vals
    signif_vec <- x$significant[[i]] %||% rep(NA, length(stats))

    rows[[i]] <- do.call(data.frame, c(list(
      contrast = rep(contrast_names[[i]], length(stats)),
      component = cluster_ids,
      statistic = as.numeric(stats),
      p_value = as.numeric(p_vals),
      p_adjusted = as.numeric(padj),
      significant = as.logical(signif_vec)
    ), dots, list(stringsAsFactors = stringsAsFactors)))
  }

  rows <- Filter(Negate(is.null), rows)
  result <- if (length(rows)) {
    do.call(rbind, rows)
  } else {
    do.call(data.frame, c(list(
      contrast = character(0),
      component = character(0),
      statistic = numeric(0),
      p_value = numeric(0),
      p_adjusted = numeric(0),
      significant = logical(0)
    ), dots, list(stringsAsFactors = stringsAsFactors)))
  }

  result$alpha <- rep(alpha_val, nrow(result))
  result$method <- rep(method_val, nrow(result))
  result$inference <- rep(inference_val, nrow(result))
  result$correction <- rep(correction_val, nrow(result))

  if (!is.null(row.names)) {
    rownames(result) <- row.names
  } else {
    rownames(result) <- NULL
  }

  result
}

# ---- Freedman-Lane scaffolding (heavy; requires time-series & GLM adapter) ----

#' Freedman-Lane permutations for DKGE (scaffold)
#'
#' This function orchestrates Freedman-Lane permutations at the *time-series* level:
#' for each subject, fit the reduced model (without the effect of interest), permute residuals,
#' reconstruct surrogate data, refit the full GLM to get B* betas, then re-run DKGE LOSO
#' to obtain a group statistic (e.g., max-|t| over medoid clusters). It requires the caller
#' to provide three adapter functions (or rely on 'fmrireg'/'neuroim2'):
#'   - fit_glm(Y_s, X_s, X0_s) -> list(beta = qxP, beta0 = q0xP, resid = TxP)
#'   - resample_resid(resid_s) -> resid_s* (TxP)  (permute or phase-randomize per run)
#'   - transport_and_stat(B_list, X_list, K, c) -> scalar (e.g., max-|t|)
#'
#' @param Y_list list of neuroim2 BrainVectors (or TxP matrices) per subject
#' @param X_list list of Txq design matrices (full)
#' @param X0_list list of Txq0 reduced designs (null space of the contrast)
#' @param K design kernel (qxq)
#' @param c contrast vector (qx1)
#' @param B number of permutations
#' @param adapters list with functions: fit_glm, resample_resid, transport_and_stat
#' @param seed RNG seed for reproducibility
#' @return list with fields: stat_obs, stat_null (B-vector), p, details
#' @export
dkge_freedman_lane <- function(Y_list, X_list, X0_list, K, c, B = 500,
                               adapters, seed = 123L) {
  stopifnot(length(Y_list) == length(X_list), length(X0_list) == length(X_list))
  B <- .dkge_validate_permutation_count(B, "B", minimum = 1L)
  set.seed(seed)
  S <- length(Y_list)

  # ---- 1) Fit reduced and full models once, compute observed pipeline stat ----
  message("Fitting observed data (full GLM per subject)...")
  fit_full <- lapply(seq_len(S), function(s) adapters$fit_glm(Y_list[[s]], X_list[[s]], X0_list[[s]]))
  B_list <- lapply(fit_full, `[[`, "beta")
  # Build observed DKGE fit and LOSO contrasts; user supplies the stat function
  stat_obs <- adapters$transport_and_stat(B_list, X_list, K, c)

  # ---- 2) Freedman-Lane permutations ----
  stat_null <- numeric(B)
  message("Running Freedman-Lane permutations...")
  for (b in seq_len(B)) {
    if (b %% max(1, B %/% 10) == 0) message(sprintf("  perm %d / %d", b, B))
    Bperm <- vector("list", S)
    for (s in seq_len(S)) {
      fs <- fit_full[[s]]
      # Y*_s = X_s beta0_s + P resid_s
      res_star <- adapters$resample_resid(fs$resid)  # TxP
      Ystar    <- X_list[[s]] %*% rbind(fs$beta0, matrix(0, nrow = ncol(X_list[[s]]) - nrow(fs$beta0), ncol = ncol(fs$beta0))) + res_star
      # refit full GLM on Ystar to get B*
      fs2 <- adapters$fit_glm(Ystar, X_list[[s]], X0_list[[s]])
      Bperm[[s]] <- fs2$beta
    }
    stat_null[b] <- adapters$transport_and_stat(Bperm, X_list, K, c)
  }

  # ---- 3) p-value (upper-tail by default if stat = max-|t|) ----
  p <- (1 + sum(stat_null >= stat_obs)) / (B + 1)
  list(stat_obs = stat_obs, stat_null = stat_null, p = p,
       details = list(seed = seed))
}
