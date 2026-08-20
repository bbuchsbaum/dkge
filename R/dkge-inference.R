
# dkge-inference.R
# Unified inference module for DKGE contrasts with multiple testing procedures

#' One-sample sign-flip max-T inference on transported subject maps
#'
#' Computes cluster-wise one-sample t-statistics across subjects on transported
#' values (SxQ matrix), and calibrates p-values by the max-|t| distribution under
#' random subject-wise sign flips (symmetric null). This does not re-estimate DKGE,
#' leveraging LOSO independence of each subject's value.
#'
#' @param Y SxQ matrix of subject values on the medoid parcellation (rows=subjects, cols=clusters)
#' @param B number of sign-flip permutations
#' @param center Location statistic. The beta API supports only `"mean"`,
#'   matching the one-sample t statistic used for observed and randomized data.
#' @param tail "two.sided" | "greater" | "less"
#' @param flips Optional precomputed S-by-B matrix of -1/+1 signs. This is an
#'   advanced reproducibility hook used to make serial and parallel execution
#'   consume exactly the same randomization descriptors.
#' @return A list with fields: `stat` (Q-vector of observed t-statistics), `p`
#'   (Q-vector of max-T family-wise-error adjusted p-values), `p_unadj`
#'   (Q-vector of per-column unadjusted permutation p-values), `maxnull`
#'   (B-vector of permutation maximum statistics), and `flips` (S-by-B sign
#'   matrix). Statistic and p-value names follow `colnames(Y)` (or stable
#'   `feature*` defaults); flip rows follow `rownames(Y)` (or `subject*`
#'   defaults).
#' @export
dkge_signflip_maxT <- function(Y, B = 2000, center = "mean",
                               tail = c("two.sided","greater","less"),
                               flips = NULL) {
  if (!is.character(center) || length(center) != 1L ||
      !identical(center, "mean")) {
    .dkge_abort(
      "`center` must be \"mean\"; median-centered sign-flip inference is not implemented.",
      "dkge_inference_compatibility_error"
    )
  }
  B <- .dkge_validate_resample_B(B)
  tail <- match.arg(tail)
  Y <- as.matrix(Y); S <- nrow(Y); Q <- ncol(Y)
  stopifnot(S >= 5, B >= 100)
  subject_ids <- rownames(Y) %||% paste0("subject", seq_len(S))
  feature_ids <- colnames(Y) %||% paste0("feature", seq_len(Q))
  permutation_ids <- paste0("perm", seq_len(B))

  # observed t per cluster
  mu  <- colMeans(Y)
  sdv <- apply(Y, 2, stats::sd)
  t_obs <- mu / (sdv / sqrt(S) + 1e-12)

  # generate random sign matrix (SxB)
  if (is.null(flips)) {
    flips <- matrix(sample(c(-1,1), S*B, replace=TRUE), S, B)
  } else {
    flips <- as.matrix(flips)
    if (!identical(dim(flips), c(S, as.integer(B))) ||
        any(!flips %in% c(-1, 1))) {
      .dkge_abort(
        sprintf("`flips` must be a %d-by-%d matrix containing only -1 and 1.", S, B),
        "dkge_inference_compatibility_error"
      )
    }
  }
  dimnames(flips) <- list(subject_ids, permutation_ids)
  Yc <- Y

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
    Yb <- flips[,b] * Yc
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

  list(stat = t_obs, p = p, p_unadj = p_unadj, maxnull = maxnull, flips = flips)
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
#'   - `"parametric"`: Parametric t-test (assumes normality)
#' @param correction Multiple testing correction:
#'   - `"maxT"`: Family-wise error rate via max-T (default)
#'   - `"fdr"`: False discovery rate (Benjamini-Hochberg)
#'   - `"bonferroni"`: Bonferroni correction
#'   - `"none"`: No correction
#' @param n_perm Number of permutations for non-parametric tests
#' @param alpha Significance level for corrections
#' @param center Location statistic for sign-flip inference. Only `"mean"` is
#'   implemented in the beta API.
#' @param parallel Logical; compute target-level inference through the parallel
#'   apply backend. Randomization descriptors are generated serially first, so
#'   serial and parallel results are identical for the same caller RNG state.
#' @param transported Logical; retained for backwards compatibility. Deprecated.
#' @param transport Optional list describing how to map subject clusters to a
#'   shared reference before inference. Provide `centroids`, `medoid`, and an
#'   optional mapper specification created via [dkge_mapper_spec()]. Additional
#'   parameters (e.g. `epsilon`, `lambda_emb`) are forwarded when constructing
#'   the default Sinkhorn mapper. Transported inference also requires an
#'   explicit `provenance` from [dkge_transport_provenance()]. A
#'   `fully_recomputed` declaration must provide
#'   `randomization_recompute(signs, target_index, target_name,
#'   contrast_results, observed_values)`, returning the randomized transported
#'   subject-by-feature matrix after rebuilding every data-dependent step.
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
#' When subject cluster counts differ across participants, supply the
#' `transport` argument so that contrasts are first mapped to a shared
#' parcellation before stacking. The resulting mapped subject matrices are
#' returned in the `transport` field of the output for downstream inspection.
#'
#' The workflow is:
#' 1. Compute contrast values via cross-fitting (LOSO/K-fold/analytic)
#' 2. Apply an implemented inference procedure (sign-flip or parametric)
#' 3. Apply multiple testing correction (maxT/FDR/Bonferroni)
#'
#' For sign-flip inference, the max-T correction provides strong FWER control.
#' For parametric inference, FDR may be more appropriate for exploratory analyses.
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
#' results <- dkge_infer(fit, c(1, rep(0, 4)), n_perm = 100)
#' results
#' }
#'
#' @seealso [dkge_contrast()], [dkge_signflip_maxT()]
#' @export
dkge_infer <- function(fit, contrasts,
                      method = c("loso", "kfold", "analytic"),
                      inference = c("signflip", "parametric"),
                      correction = c("maxT", "fdr", "bonferroni", "none"),
                      n_perm = 2000,
                      alpha = 0.05,
                      center = "mean",
                      parallel = FALSE,
                      transported = FALSE,
                      transport = NULL,
                      ...) {
  method <- match.arg(method)
  if (missing(inference)) {
    inference <- "signflip"
  } else if (!is.character(inference) || length(inference) != 1L) {
    .dkge_abort("`inference` must name one inference method.",
                "dkge_inference_compatibility_error")
  }
  correction <- match.arg(correction)
  n_perm <- .dkge_validate_resample_B(n_perm)
  alpha <- .dkge_validate_probability(alpha, "alpha")
  transport_recompute <- transport$randomization_recompute %||% NULL
  transport_provenance <- if (is.null(transport)) {
    NULL
  } else {
    transport$provenance %||% .dkge_data_derived_loading_provenance()
  }

  .dkge_validate_inference_compatibility(
    inference = inference,
    correction = correction,
    has_transport = !is.null(transport),
    transport_provenance = transport_provenance,
    transport_recompute = transport_recompute,
    center = center,
    n_targets = NA_integer_,
    parallel = parallel
  )

  if (!is.null(transported) && transported) {
    warning("`transported` argument is deprecated; supply transported matrices directly to inference helpers if needed.",
            call. = FALSE)
  }

  contrast_results <- dkge_contrast(
    fit, contrasts, method = method, parallel = parallel, ...
  )
  .dkge_validate_inference_compatibility(
    inference = inference,
    correction = correction,
    has_transport = !is.null(transport),
    transport_provenance = transport_provenance,
    transport_recompute = transport_recompute,
    center = center,
    n_targets = length(contrast_results$contrasts),
    parallel = parallel
  )

  mapped_values <- NULL
  transport_results <- NULL
  if (!is.null(transport) && !is.null(transport$centroids)) {
    medoid <- transport$medoid %||% 1L
    mapper_spec <- transport$mapper %||% NULL
    method_arg <- transport$method %||% "sinkhorn"
    mapper_args <- transport[intersect(names(transport),
                                       c("epsilon", "max_iter", "tol",
                                         "lambda_emb", "lambda_spa",
                                         "sigma_mm", "lambda_size",
                                         "value_type", "warm_start"))]
    args <- list(
      fit = fit,
      contrast_obj = contrast_results,
      medoid = medoid,
      centroids = transport$centroids,
      loadings = transport$loadings,
      betas = transport$betas,
      sizes = transport$sizes,
      mapper = mapper_spec,
      method = method_arg,
      transport_cache = transport$transport_cache,
      provenance = transport_provenance
    )
    args <- c(args, mapper_args)
    args <- args[!vapply(args, is.null, logical(1))]
    transport_results <- do.call(dkge_transport_contrasts_to_medoid, args)
    mapped_values <- lapply(transport_results, `[[`, "subj_values")
  }

  # Step 2: Apply inference procedure
  infer_results <- switch(inference,
    signflip = .infer_signflip(
      contrast_results, n_perm, correction, mapped_values,
      center = center, parallel = parallel,
      transport_recompute = transport_recompute,
      transport_provenance = transport_provenance
    ),
    parametric = .infer_parametric(contrast_results, correction)
  )

  # Step 3: Apply correction (if not already done by inference method)
  if (inference != "signflip" || correction != "maxT") {
    infer_results <- .apply_correction(infer_results, correction, alpha)
  }

  # Mark significant results
  infer_results$significant <- lapply(infer_results$p_adjusted, function(p) p <= alpha)
  infer_results$alpha <- alpha

  if (!is.null(transport_results)) {
  infer_results$transport <- transport_results
  }

  structure(infer_results, class = "dkge_inference")
}

#' Validate the complete dkge_infer compatibility contract
#' @keywords internal
#' @noRd
.dkge_validate_inference_compatibility <- function(
    inference,
    correction,
    has_transport,
    transport_provenance = NULL,
    transport_recompute = NULL,
    center,
    n_targets,
    parallel) {
  abort <- function(message) {
    .dkge_abort(message, "dkge_inference_compatibility_error")
  }
  if (!is.character(inference) || length(inference) != 1L ||
      !inference %in% c("signflip", "parametric", "freedman-lane")) {
    abort("`inference` must be \"signflip\" or \"parametric\".")
  }
  if (identical(inference, "freedman-lane")) {
    abort("Freedman-Lane is not implemented by `dkge_infer()` and is not a supported beta choice.")
  }
  if (!is.character(correction) || length(correction) != 1L ||
      !correction %in% c("maxT", "fdr", "bonferroni", "none")) {
    abort("`correction` must be one of maxT, fdr, bonferroni, or none.")
  }
  if (!is.logical(has_transport) || length(has_transport) != 1L ||
      is.na(has_transport)) {
    abort("`has_transport` must be one TRUE/FALSE value.")
  }
  if (!is.logical(parallel) || length(parallel) != 1L || is.na(parallel)) {
    abort("`parallel` must be one TRUE/FALSE value.")
  }
  if (!is.character(center) || length(center) != 1L ||
      !identical(center, "mean")) {
    abort("`center` must be \"mean\"; no other sign-flip location statistic is implemented.")
  }
  if (identical(inference, "parametric") && identical(correction, "maxT")) {
    abort("Parametric inference does not implement maxT correction; use fdr, bonferroni, or none.")
  }
  if (identical(inference, "parametric") && has_transport) {
    abort("Parametric inference with transport is unsupported; use signflip or omit transport.")
  }
  if (length(n_targets) != 1L ||
      (!is.na(n_targets) &&
       (!is.numeric(n_targets) || !is.finite(n_targets) ||
        n_targets < 1 || n_targets != as.integer(n_targets)))) {
    abort("`n_targets` must be a positive integer when known.")
  }
  if (has_transport) {
    if (is.null(transport_provenance)) {
      .dkge_abort(
        "Transported inference requires explicit transport provenance; descriptive loading-derived transport is not inferentially valid by default.",
        "dkge_transport_inference_error"
      )
    }
    .dkge_validate_transport_inference_provenance(
      transport_provenance,
      randomization_recompute = transport_recompute
    )
  }
  invisible(TRUE)
}

#' Sign-flip inference helper
#' @keywords internal
#' @noRd
.infer_signflip <- function(contrast_results, n_perm, correction,
                            mapped_values = NULL, center = "mean",
                            parallel = FALSE,
                            transport_recompute = NULL,
                            transport_provenance = NULL) {
  n_contrasts <- length(contrast_results$contrasts)
  target_names <- names(contrast_results$contrasts)
  if (is.null(target_names) || any(!nzchar(target_names))) {
    target_names <- paste0("target", seq_len(n_contrasts))
  }
  Y_list <- lapply(seq_len(n_contrasts), function(i) {
    Y <- if (!is.null(mapped_values)) {
      mapped_values[[i]]
    } else {
      as.matrix(contrast_results, contrast = i)
    }
    Y <- as.matrix(Y)
    if (is.null(colnames(Y))) {
      colnames(Y) <- paste0("feature", seq_len(ncol(Y)))
    }
    Y
  })
  names(Y_list) <- target_names

  # Draw all randomization descriptors before selecting an execution backend.
  flip_list <- lapply(Y_list, function(Y) {
    matrix(sample(c(-1, 1), nrow(Y) * n_perm, replace = TRUE),
           nrow(Y), n_perm)
  })

  one_target <- function(i) {
    Y <- Y_list[[i]]
    n_subjects <- nrow(Y)

    randomized_values <- function(b) {
      signs <- flip_list[[i]][, b]
      if (is.null(transport_recompute)) {
        return(signs * Y)
      }
      Y_randomized <- transport_recompute(
        signs = signs,
        target_index = i,
        target_name = target_names[[i]],
        contrast_results = contrast_results,
        observed_values = Y
      )
      Y_randomized <- as.matrix(Y_randomized)
      if (!is.numeric(Y_randomized) || any(!is.finite(Y_randomized)) ||
          !identical(dim(Y_randomized), dim(Y))) {
        .dkge_abort(
          sprintf(
            "Transport recomputation for target '%s' must return a finite numeric %d-by-%d matrix.",
            target_names[[i]], nrow(Y), ncol(Y)
          ),
          "dkge_transport_inference_error"
        )
      }
      dimnames(Y_randomized) <- dimnames(Y)
      Y_randomized
    }

    if (correction == "maxT") {
      if (is.null(transport_recompute)) {
        result <- dkge_signflip_maxT(
          Y, B = n_perm, center = center, flips = flip_list[[i]]
        )
      } else {
        mu <- colMeans(Y)
        se <- apply(Y, 2, stats::sd) / sqrt(pmax(n_subjects, 1))
        statistic <- mu / (se + 1e-12)
        names(statistic) <- colnames(Y)
        null_stats <- vapply(seq_len(n_perm), function(b) {
          Yb <- randomized_values(b)
          colMeans(Yb) /
            (apply(Yb, 2, stats::sd) / sqrt(pmax(n_subjects, 1)) + 1e-12)
        }, numeric(ncol(Y)))
        if (is.null(dim(null_stats))) {
          null_stats <- matrix(null_stats, nrow = ncol(Y))
        }
        maxnull <- apply(abs(null_stats), 2, max)
        p_unadj <- vapply(seq_along(statistic), function(j) {
          (1 + sum(abs(null_stats[j, ]) >= abs(statistic[j]))) / (n_perm + 1)
        }, numeric(1))
        p_adjusted <- vapply(abs(statistic), function(value) {
          (1 + sum(maxnull >= value)) / (n_perm + 1)
        }, numeric(1))
        names(p_unadj) <- names(statistic)
        names(p_adjusted) <- names(statistic)
        result <- list(stat = statistic, p_unadj = p_unadj, p = p_adjusted)
      }
      return(list(
        statistic = result$stat,
        p_value = result$p_unadj %||% result$p,
        p_adjusted = result$p
      ))
    }

    mu <- colMeans(Y)
    se <- apply(Y, 2, stats::sd) / sqrt(pmax(n_subjects, 1))
    statistic <- mu / (se + 1e-12)
    names(statistic) <- colnames(Y)

    null_stats <- vapply(seq_len(n_perm), function(b) {
      Y_flip <- randomized_values(b)
      mu_flip <- colMeans(Y_flip)
      se_flip <- apply(Y_flip, 2, stats::sd) / sqrt(pmax(n_subjects, 1))
      mu_flip / (se_flip + 1e-12)
    }, numeric(ncol(Y)))
    if (is.null(dim(null_stats))) {
      null_stats <- matrix(null_stats, nrow = ncol(Y))
    }

    p_value <- vapply(seq_along(statistic), function(j) {
      (1 + sum(abs(null_stats[j, ]) >= abs(statistic[j]))) / (n_perm + 1)
    }, numeric(1))
    names(p_value) <- names(statistic)
    list(statistic = statistic, p_value = p_value, p_adjusted = p_value)
  }

  target_results <- .dkge_apply(
    seq_len(n_contrasts), one_target, parallel = parallel
  )
  stats <- setNames(lapply(target_results, `[[`, "statistic"), target_names)
  p_values <- setNames(lapply(target_results, `[[`, "p_value"), target_names)
  p_adjusted <- setNames(lapply(target_results, `[[`, "p_adjusted"), target_names)

  list(
    contrasts = contrast_results,
    statistics = stats,
    p_values = p_values,
    p_adjusted = p_adjusted,
    method = contrast_results$method,
    inference = "signflip",
    correction = correction,
    exactness = transport_provenance$exactness %||% "conditional_exact",
    metadata = list(
      n_perm = n_perm,
      center = center,
      parallel = parallel,
      transport_provenance = transport_provenance
    )
  )
}

#' Parametric inference helper
#' @keywords internal
#' @noRd
.infer_parametric <- function(contrast_results, correction) {
  n_contrasts <- length(contrast_results$contrasts)
  target_names <- names(contrast_results$contrasts)
  if (is.null(target_names) || any(!nzchar(target_names))) {
    target_names <- paste0("target", seq_len(n_contrasts))
  }

  stats <- vector("list", n_contrasts)
  p_values <- vector("list", n_contrasts)
  df_vec <- numeric(n_contrasts)

  for (i in seq_len(n_contrasts)) {
    Y <- as.matrix(contrast_results, contrast = i)
    n_subjects <- nrow(Y)
    df_vec[i] <- n_subjects - 1
    mu <- colMeans(Y)
    se <- apply(Y, 2, stats::sd) / sqrt(n_subjects)
    t_stats <- mu / (se + 1e-12)
    feature_names <- colnames(Y) %||% paste0("feature", seq_len(ncol(Y)))
    names(t_stats) <- feature_names

    stats[[i]] <- t_stats
    p_values[[i]] <- 2 * pt(-abs(t_stats), df_vec[i])
    names(p_values[[i]]) <- feature_names
  }
  names(stats) <- target_names
  names(p_values) <- target_names

  list(
    contrasts = contrast_results,
    statistics = stats,
    p_values = p_values,
    p_adjusted = p_values,  # Will be corrected next
    method = contrast_results$method,
    inference = "parametric",
    correction = correction,
    exactness = "asymptotic",
    metadata = list(df = df_vec)
  )
}

#' Apply multiple testing correction
#' @keywords internal
#' @noRd
.apply_correction <- function(infer_results, correction, alpha) {
  if (correction == "none") {
    return(infer_results)
  }

  n_contrasts <- length(infer_results$p_values)

  for (i in seq_len(n_contrasts)) {
    p <- infer_results$p_values[[i]]

    infer_results$p_adjusted[[i]] <- switch(correction,
      fdr = p.adjust(p, method = "fdr"),
      bonferroni = pmin(p * length(p), 1),
      none = p
    )
  }

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
#' @return `x`, invisibly.
#' @method print dkge_inference
#' @export
print.dkge_inference <- function(x, ...) {
  cat("DKGE Inference Results\n")
  cat("----------------------\n")
  cat(sprintf("Cross-fitting: %s\n", x$method))
  cat(sprintf("Inference: %s\n", x$inference))
  cat(sprintf("Correction: %s\n", x$correction))

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
#'   - resample_resid(resid_s) -> resid_s* (TxP)  [permute or phase-randomize per run]
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
