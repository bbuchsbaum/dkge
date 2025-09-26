
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
#' @param center "mean" or "median" for the location statistic (t uses mean)
#' @param tail "two.sided" | "greater" | "less"
#' @return list with fields: stat (Q-vector), p (Q-vector), maxnull (B-vector), flips (SxB signs)
#' @export
dkge_signflip_maxT <- function(Y, B = 2000, center = c("mean","median"),
                               tail = c("two.sided","greater","less")) {
  center <- match.arg(center); tail <- match.arg(tail)
  Y <- as.matrix(Y); S <- nrow(Y); Q <- ncol(Y)
  stopifnot(S >= 5, B >= 100)

  # observed t per cluster
  mu  <- colMeans(Y)
  sdv <- apply(Y, 2, stats::sd)
  t_obs <- mu / (sdv / sqrt(S) + 1e-12)

  # generate random sign matrix (SxB)
  flips <- matrix(sample(c(-1,1), S*B, replace=TRUE), S, B)
  # center data if using "median" just for display; t uses mean in any case
  Yc <- Y

  # permutation max stat
  maxnull <- numeric(B)
  for (b in seq_len(B)) {
    Yb <- flips[,b] * Yc
    mu_b  <- colMeans(Yb)
    sd_b  <- apply(Yb, 2, stats::sd)
    t_b   <- mu_b / (sd_b / sqrt(S) + 1e-12)
    maxnull[b] <- max(abs(t_b))
  }

  # p-values (max-T, strong FWER control)
  if (tail == "two.sided") {
    p <- sapply(abs(t_obs), function(x) (1 + sum(maxnull >= x)) / (B + 1))
  } else if (tail == "greater") {
    p <- sapply(t_obs, function(x) (1 + sum(maxnull >= x)) / (B + 1)) # conservative
  } else {
    p <- sapply(-t_obs, function(x) (1 + sum(maxnull >= x)) / (B + 1))
  }

  list(stat = t_obs, p = p, maxnull = maxnull, flips = flips)
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
#' @param transported Logical; retained for backwards compatibility. Deprecated.
#' @param transport Optional list describing how to map subject clusters to a
#'   shared reference before inference. Provide `centroids`, `medoid`, and an
#'   optional mapper specification created via [dkge_mapper_spec()]. Additional
#'   parameters (e.g. `epsilon`, `lambda_emb`) are forwarded when constructing
#'   the default Sinkhorn mapper.
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
#' 2. Apply inference procedure (sign-flip/Freedman-Lane/parametric)
#' 3. Apply multiple testing correction (maxT/FDR/Bonferroni)
#'
#' For sign-flip inference, the max-T correction provides strong FWER control.
#' For parametric inference, FDR may be more appropriate for exploratory analyses.
#'
#' @examples
#' \dontrun{
#' # Standard LOSO with sign-flip and maxT correction
#' results <- dkge_infer(fit, c(1, -1, 0, 0, 0))
#'
#' # K-fold with parametric test and FDR
#' results <- dkge_infer(fit, c(1, -1, 0, 0, 0),
#'                      method = "kfold", folds = 5,
#'                      inference = "parametric",
#'                      correction = "fdr")
#'
#' # Multiple contrasts with analytic approximation
#' contrasts <- list(
#'   main1 = c(1, -1, 0, 0, 0),
#'   main2 = c(0, 0, 1, -1, 0)
#' )
#' results <- dkge_infer(fit, contrasts, method = "analytic")
#'
#' # Extract significant clusters for first contrast
#' sig_clusters <- which(results$significant[[1]])
#' }
#'
#' @seealso [dkge_contrast()], [dkge_signflip_maxT()], [dkge_freedman_lane()]
#' @export
dkge_infer <- function(fit, contrasts,
                      method = c("loso", "kfold", "analytic"),
                      inference = c("signflip", "freedman-lane", "parametric"),
                      correction = c("maxT", "fdr", "bonferroni", "none"),
                      n_perm = 2000,
                      alpha = 0.05,
                      transported = FALSE,
                      transport = NULL,
                      ...) {
  method <- match.arg(method)
  inference <- match.arg(inference)
  correction <- match.arg(correction)

  if (!is.null(transported) && transported) {
    warning("`transported` argument is deprecated; supply transported matrices directly to inference helpers if needed.",
            call. = FALSE)
  }

  contrast_results <- dkge_contrast(fit, contrasts, method = method, ...)

  mapped_values <- NULL
  transport_results <- NULL
  if (!is.null(transport) && !is.null(transport$centroids)) {
    medoid <- transport$medoid %||% 1L
    mapper_spec <- transport$mapper %||% NULL
    method_arg <- transport$method %||% "sinkhorn"
    mapper_args <- transport[intersect(names(transport),
                                       c("epsilon", "max_iter", "tol",
                                         "lambda_emb", "lambda_spa",
                                         "sigma_mm", "lambda_size"))]
    args <- list(
      fit = fit,
      contrast_obj = contrast_results,
      medoid = medoid,
      centroids = transport$centroids,
      loadings = transport$loadings,
      betas = transport$betas,
      sizes = transport$sizes,
      mapper = mapper_spec,
      method = method_arg
    )
    args <- c(args, mapper_args)
    args <- args[!vapply(args, is.null, logical(1))]
    transport_results <- do.call(dkge_transport_contrasts_to_medoid, args)
    mapped_values <- lapply(transport_results, `[[`, "subj_values")
  }

  # Step 2: Apply inference procedure
  infer_results <- switch(inference,
    signflip = .infer_signflip(contrast_results, n_perm, correction, mapped_values),
    `freedman-lane` = .infer_freedman_lane(contrast_results, n_perm, correction, ...),
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

#' Sign-flip inference helper
#' @keywords internal
#' @noRd
.infer_signflip <- function(contrast_results, n_perm, correction, mapped_values = NULL) {
  n_contrasts <- length(contrast_results$contrasts)

  stats <- vector("list", n_contrasts)
  p_values <- vector("list", n_contrasts)
  p_adjusted <- vector("list", n_contrasts)

  for (i in seq_len(n_contrasts)) {
    Y <- if (!is.null(mapped_values)) {
      mapped_values[[i]]
    } else {
      as.matrix(contrast_results, contrast = i)
    }
    n_subjects <- nrow(Y)

    if (correction == "maxT") {
      result <- dkge_signflip_maxT(Y, B = n_perm)
      stats[[i]] <- result$stat
      p_values[[i]] <- result$p
      p_adjusted[[i]] <- result$p
      next
    }

    mu <- colMeans(Y)
    se <- apply(Y, 2, stats::sd) / sqrt(pmax(n_subjects, 1))
    stats[[i]] <- mu / (se + 1e-12)

    flips <- matrix(sample(c(-1, 1), n_subjects * n_perm, replace = TRUE),
                    n_subjects, n_perm)
    null_stats <- apply(flips, 2, function(s) {
      Y_flip <- s * Y
      mu_flip <- colMeans(Y_flip)
      se_flip <- apply(Y_flip, 2, stats::sd) / sqrt(pmax(n_subjects, 1))
      mu_flip / (se_flip + 1e-12)
    })

    p_values[[i]] <- vapply(seq_along(stats[[i]]), function(j) {
      (1 + sum(abs(null_stats[j, ]) >= abs(stats[[i]][j]))) / (n_perm + 1)
    }, numeric(1))

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
    metadata = list(n_perm = n_perm)
  )
}

#' Parametric inference helper
#' @keywords internal
#' @noRd
.infer_parametric <- function(contrast_results, correction) {
  n_contrasts <- length(contrast_results$contrasts)

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

    stats[[i]] <- t_stats
    p_values[[i]] <- 2 * pt(-abs(t_stats), df_vec[i])
  }

  list(
    contrasts = contrast_results,
    statistics = stats,
    p_values = p_values,
    p_adjusted = p_values,  # Will be corrected next
    method = contrast_results$method,
    inference = "parametric",
    correction = correction,
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
#' @param ... Additional arguments passed to [base::data.frame()]
#' @param stringsAsFactors Logical; forwarded to [base::data.frame()]
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
