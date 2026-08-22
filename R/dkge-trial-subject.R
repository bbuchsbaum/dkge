# dkge-trial-subject.R
# Trialwise beta input reduced to q-space sufficient statistics.

#' Build a DKGE subject from trialwise beta maps
#'
#' Fits the second-stage model `Y = X B + E` without constructing any
#' voxel-by-voxel covariance. By default the estimator is IID OLS. Supplying a
#' trial covariance or precision uses GLS and stores `(X' W X)^-1` as the
#' unit-spatial-variance covariance of the effect estimates. The returned
#' subject retains q-space sufficient statistics and the shared `T x q`
#' design, but never the trialwise `y`. Use [dkge_trial_subject_chunks()] when
#' `T x P` itself is too large to materialize.
#'
#' @param y Numeric `T x P` matrix of trialwise beta maps.
#' @param design Numeric `T x q` second-stage design matrix. This may be a
#'   one-hot cell design or a general full-rank basis design.
#' @param id Optional subject identifier.
#' @param omega Optional voxel/parcel weighting passed to [dkge_subject()].
#' @param effect_precision Optional direct q-vector of effect precisions, or
#'   `"split_half"` to derive bounded effect precision from the stored halves
#'   with [dkge_split_effect_precision()].
#' @param trial_covariance Optional symmetric positive-definite `T x T` relative
#'   covariance of trial errors. Mutually exclusive with `trial_precision`.
#' @param trial_precision Optional symmetric positive-definite `T x T` trial
#'   precision, or a function mapping a `T x k` matrix to its precision-weighted
#'   counterpart. A function is materialised only in trial space for validation;
#'   no `P x P` matrix is formed.
#' @param effect_noise_cov Optional directly supplied q-by-q covariance
#'   multiplier for the estimated effects. When supplied it overrides the
#'   covariance implied by the OLS or GLS information matrix.
#' @param residual_variance Optional externally estimated length-P spatial noise
#'   variances. The default is the (weighted for GLS) residual sum of squares
#'   divided by the residual degrees of freedom.
#' @param noise_trace Optional externally supplied unweighted spatial noise
#'   trace. By default this is `sum(residual_variance)` when estimable.
#'   A supplied value is used by analytic debiasing when no extra spatial
#'   weights are applied.
#' @param split Optional split-half sufficient statistic. `"within_cell"`
#'   deterministically balances trials within each one-hot cell; `"alternate"`
#'   alternates all trial rows; `"run"` assigns whole runs to halves; and
#'   `"explicit"` uses `split_labels`. Both half-designs must remain full rank.
#' @param split_labels Optional length-T vector containing exactly two explicit
#'   half labels. Supplying it selects `split = "explicit"`.
#' @param run_labels Optional length-T run labels. With `split = "run"`, whole
#'   runs are assigned deterministically to halves. With another split mode,
#'   these labels are used to audit whether the halves are run-disjoint.
#' @param split_independent Logical; explicitly declare a non-run-disjoint
#'   partition independent. Without run-disjoint labels or this declaration,
#'   split moments are recorded as exploratory and a warning is emitted.
#' @param tol Numerical tolerance used for rank and one-hot checks.
#'
#' @return A `dkge_subject` carrying second-stage sufficient statistics.
#' @details
#' The analytic noise correction assumes a common trial-dependence matrix
#' across spatial columns and a separable error covariance: the effect
#' covariance is `(X' W X)^-1`, scaled by each spatial column's residual
#' variance. Trialwise beta maps from overlapping or temporally correlated
#' first-level estimates generally violate IID errors; supply
#' `trial_covariance` or `trial_precision` in that case. Using the IID default
#' when trial errors are correlated mis-scales both the analytic subtraction
#' and any precision derived from it.
#'
#' `split = "within_cell"` and `split = "alternate"` only create deterministic
#' partitions; they do not establish independence. A split cross-moment is
#' labeled independent only for run-disjoint halves or when
#' `split_independent = TRUE` is explicitly justified by the caller.
#' @export
#' @examples
#' set.seed(1)
#' X <- model.matrix(~ 0 + factor(rep(1:3, each = 4)))
#' Y <- X %*% matrix(rnorm(3 * 5), 3, 5) + matrix(rnorm(12 * 5), 12, 5)
#' subject <- dkge_trial_subject(Y, X, id = "s1")
dkge_trial_subject <- function(y, design, id = NULL, omega = NULL,
                               effect_precision = NULL,
                               trial_covariance = NULL,
                               trial_precision = NULL,
                               effect_noise_cov = NULL,
                               residual_variance = NULL,
                               noise_trace = NULL,
                               split = c("none", "within_cell", "alternate",
                                         "run", "explicit"),
                               split_labels = NULL,
                               run_labels = NULL,
                               split_independent = FALSE,
                               tol = 1e-8) {
  split_was_missing <- missing(split)
  split <- match.arg(split)
  split <- .dkge_resolve_trial_split(
    split, split_was_missing, split_labels, run_labels
  )
  y <- as.matrix(y)
  design <- as.matrix(design)
  if (!is.numeric(y) || !is.numeric(design) || any(!is.finite(y)) ||
      any(!is.finite(design))) {
    stop("`y` and `design` must be finite numeric matrices.", call. = FALSE)
  }
  if (nrow(y) != nrow(design)) {
    stop("`y` and `design` must have the same number of trial rows.",
         call. = FALSE)
  }
  if (!nrow(y) || !ncol(y) || !ncol(design)) {
    stop("`y` and `design` must have positive dimensions.", call. = FALSE)
  }
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol <= 0) {
    stop("`tol` must be a positive finite scalar.", call. = FALSE)
  }
  if (!is.logical(split_independent) || length(split_independent) != 1L ||
      is.na(split_independent)) {
    stop("`split_independent` must be TRUE or FALSE.", call. = FALSE)
  }
  if (identical(split, "explicit") && is.null(split_labels)) {
    stop("split='explicit' requires `split_labels`.", call. = FALSE)
  }
  if (!is.null(run_labels) && length(run_labels) != nrow(design)) {
    stop("`run_labels` must have one entry per trial row.", call. = FALSE)
  }
  if (!is.null(run_labels) && anyNA(run_labels)) {
    stop("`run_labels` cannot contain missing values.", call. = FALSE)
  }

  split_precision_requested <- is.character(effect_precision)
  if (split_precision_requested &&
      !identical(effect_precision, "split_half")) {
    stop("Character `effect_precision` must be 'split_half'.", call. = FALSE)
  }
  if (split_precision_requested && identical(split, "none")) {
    stop("effect_precision='split_half' requires a split-half estimate.",
         call. = FALSE)
  }
  effect_precision_input <- if (split_precision_requested) NULL else effect_precision

  effects <- colnames(design) %||% paste0("effect_", seq_len(ncol(design)))
  clusters <- colnames(y) %||% paste0("cluster_", seq_len(ncol(y)))
  colnames(design) <- effects
  colnames(y) <- clusters

  q <- ncol(design)
  qrx <- qr(design, tol = tol)
  if (qrx$rank < q) {
    stop(sprintf(
      "`design` must be full column rank for trialwise reduction (rank %d < %d).",
      qrx$rank, q
    ), call. = FALSE)
  }
  weighting <- .dkge_trial_weighting(
    trial_covariance, trial_precision, nrow(design), tol
  )
  effect_information <- if (identical(weighting$type, "iid")) {
    crossprod(design)
  } else {
    crossprod(design, weighting$apply(design))
  }
  effect_information <- (effect_information + t(effect_information)) / 2
  info_chol <- tryCatch(chol(effect_information), error = function(e) NULL)
  if (is.null(info_chol)) {
    stop("The weighted trial design does not have positive-definite effect information.",
         call. = FALSE)
  }
  effect_score <- if (identical(weighting$type, "iid")) {
    crossprod(design, y)
  } else {
    crossprod(design, weighting$apply(y))
  }
  beta <- if (identical(weighting$type, "iid")) {
    qr.coef(qrx, y)
  } else {
    backsolve(info_chol, forwardsolve(t(info_chol), effect_score))
  }
  rownames(beta) <- effects
  colnames(beta) <- clusters
  residuals <- y - design %*% beta
  residual_df <- nrow(design) - qrx$rank
  if (!identical(weighting$type, "iid") && residual_df <= 0 &&
      is.null(residual_variance)) {
    stop("Covariance-aware uncertainty requires positive residual degrees of freedom.",
         call. = FALSE)
  }
  residual_sum_squares <- if (identical(weighting$type, "iid")) {
    colSums(residuals^2)
  } else {
    weighted_residuals <- weighting$apply(residuals)
    rss <- colSums(residuals * weighted_residuals)
    rss_scale <- pmax(1, colSums(abs(residuals * weighted_residuals)))
    if (any(rss < -tol * rss_scale)) {
      stop("Trial precision produced a negative residual quadratic form.",
           call. = FALSE)
    }
    pmax(rss, 0)
  }
  supplied_residual_variance <- !is.null(residual_variance)
  residual_variance <- if (supplied_residual_variance) {
    .validate_effect_vector(residual_variance, clusters, "residual_variance")
  } else if (residual_df > 0) {
    stats::setNames(residual_sum_squares / residual_df, clusters)
  } else {
    stats::setNames(rep(NA_real_, ncol(y)), clusters)
  }
  names(residual_variance) <- clusters

  supplied_effect_covariance <- !is.null(effect_noise_cov)
  effect_noise_cov <- if (supplied_effect_covariance) {
    .validate_effect_covariance(effect_noise_cov, effects)
  } else {
    chol2inv(info_chol)
  }
  dimnames(effect_noise_cov) <- list(effects, effects)

  supplied_noise_trace <- !is.null(noise_trace)
  if (supplied_noise_trace) {
    noise_trace <- as.numeric(noise_trace)
    if (length(noise_trace) != 1L || !is.finite(noise_trace) ||
        noise_trace < 0) {
      stop("`noise_trace` must be a finite non-negative scalar.",
           call. = FALSE)
    }
  } else if (all(is.finite(residual_variance))) {
    noise_trace <- sum(residual_variance)
  } else {
    noise_trace <- NULL
  }

  one_hot <- all(abs(design - round(design)) <= tol) &&
    all(round(design) %in% c(0, 1)) &&
    all(rowSums(round(design)) == 1)
  effect_n <- if (one_hot) {
    stats::setNames(as.numeric(colSums(round(design))), effects)
  } else {
    NULL
  }

  split_betas <- NULL
  split_provenance <- NULL
  if (!identical(split, "none")) {
    if (identical(split, "within_cell") && !one_hot) {
      stop("split='within_cell' requires a one-hot design.", call. = FALSE)
    }
    split_group <- if (identical(split, "within_cell")) {
      cell <- max.col(round(design), ties.method = "first")
      stats::ave(seq_along(cell), cell,
                 FUN = function(idx) ((seq_along(idx) - 1L) %% 2L) + 1L)
    } else if (identical(split, "alternate")) {
      ((seq_len(nrow(design)) - 1L) %% 2L) + 1L
    } else if (identical(split, "run")) {
      if (is.null(run_labels)) {
        stop("split='run' requires `run_labels`.", call. = FALSE)
      }
      run_levels <- unique(as.character(run_labels))
      if (length(run_levels) < 2L) {
        stop("split='run' requires at least two runs.", call. = FALSE)
      }
      ((match(as.character(run_labels), run_levels) - 1L) %% 2L) + 1L
    } else {
      labels <- as.character(split_labels)
      if (length(labels) != nrow(design) || anyNA(labels)) {
        stop("`split_labels` must contain one non-missing label per trial row.",
             call. = FALSE)
      }
      label_levels <- unique(labels)
      if (length(label_levels) != 2L) {
        stop("`split_labels` must contain exactly two distinct labels.",
             call. = FALSE)
      }
      match(labels, label_levels)
    }
    half_indices <- list(first = which(split_group == 1L),
                         second = which(split_group == 2L))
    half_ranks <- integer(2)
    effect_counts <- matrix(0, q, 2,
                            dimnames = list(effects, names(half_indices)))
    split_betas <- lapply(seq_along(half_indices), function(h) {
      idx <- half_indices[[h]]
      Xh <- design[idx, , drop = FALSE]
      Yh <- y[idx, , drop = FALSE]
      qrh <- qr(Xh, tol = tol)
      half_ranks[[h]] <<- qrh$rank
      effect_counts[, h] <<- colSums(abs(Xh) > tol)
      if (qrh$rank < q) {
        stop(sprintf(
          "The '%s' split does not leave both half-designs full rank; add trials or choose split='none'.",
          split
        ), call. = FALSE)
      }
      half_weighting <- .dkge_subset_trial_weighting(weighting, idx, tol)
      Bh <- if (identical(half_weighting$type, "iid")) {
        qr.coef(qrh, Yh)
      } else {
        info_h <- crossprod(Xh, half_weighting$apply(Xh))
        info_h <- (info_h + t(info_h)) / 2
        chol_h <- tryCatch(chol(info_h), error = function(e) NULL)
        if (is.null(chol_h)) {
          stop(sprintf(
            "The '%s' split leaves a singular weighted half-design.", split
          ), call. = FALSE)
        }
        score_h <- crossprod(Xh, half_weighting$apply(Yh))
        backsolve(chol_h, forwardsolve(t(chol_h), score_h))
      }
      rownames(Bh) <- effects
      colnames(Bh) <- clusters
      Bh
    })
    names(split_betas) <- names(half_indices)

    run_disjoint <- FALSE
    if (!is.null(run_labels)) {
      run_character <- as.character(run_labels)
      run_disjoint <- !length(intersect(
        unique(run_character[half_indices$first]),
        unique(run_character[half_indices$second])
      ))
    }
    independent <- isTRUE(run_disjoint) || isTRUE(split_independent)
    independence_basis <- if (isTRUE(run_disjoint)) {
      "run_disjoint"
    } else if (isTRUE(split_independent)) {
      "user_declared"
    } else {
      "exploratory_same_source"
    }
    if (independent && !identical(weighting$type, "iid")) {
      .dkge_assert_split_covariance_independence(
        weighting, half_indices, tol
      )
    }
    if (!independent) {
      warning(sprintf(
        paste0("split='%s' is exploratory because the halves are not verified ",
               "independent. Supply run-disjoint `run_labels` or set ",
               "`split_independent = TRUE` only when independence is justified."),
        split
      ), call. = FALSE)
    }
    split_provenance <- list(
      mode = split,
      split_labels = stats::setNames(as.integer(split_group),
                                     rownames(design) %||% seq_len(nrow(design))),
      run_labels = if (is.null(run_labels)) NULL else as.character(run_labels),
      half_indices = half_indices,
      half_ranks = stats::setNames(half_ranks, names(half_indices)),
      effect_counts = effect_counts,
      independent = independent,
      independence_basis = independence_basis,
      run_disjoint = run_disjoint,
      user_declared_independent = isTRUE(split_independent)
    )
  }

  subject <- dkge_subject(
    beta,
    design = design,
    id = id,
    omega = omega,
    effect_n = effect_n,
    effect_precision = effect_precision_input,
    effect_noise_cov = effect_noise_cov,
    residual_variance = residual_variance,
    residual_df = residual_df,
    noise_trace = noise_trace,
    split_betas = split_betas
  )
  subject$input_type <- "trialwise"
  subject$n_trials <- nrow(design)
  dimnames(effect_information) <- list(effects, effects)
  dimnames(effect_score) <- list(effects, clusters)
  names(residual_sum_squares) <- clusters
  subject$effect_information <- effect_information
  subject$effect_score <- effect_score
  subject$residual_sum_squares <- residual_sum_squares
  subject$noise_trace_scope <- if (supplied_noise_trace) {
    "supplied"
  } else {
    "unweighted"
  }
  subject$split_provenance <- split_provenance
  subject$error_model <- list(
    estimator = if (identical(weighting$type, "iid")) "ols" else "gls",
    trial_dependence = weighting$type,
    effect_covariance_source = if (supplied_effect_covariance) {
      "supplied"
    } else {
      "inverse_information"
    },
    residual_scale_source = if (supplied_residual_variance) {
      "supplied"
    } else {
      "weighted_residual_mean_square"
    },
    noise_trace_source = if (supplied_noise_trace) "supplied" else "residual_variance",
    residual_df = residual_df,
    assumptions = "shared trial dependence across spatial columns; independent spatial errors for analytic trace reduction"
  )
  if (split_precision_requested) {
    subject$effect_precision <- dkge_split_effect_precision(subject)
    subject$split_reliability <- attr(subject$effect_precision, "diagnostics")
  }
  subject
}

#' Build a trialwise DKGE subject from feature chunks
#'
#' Reduces successive `T x P_block` trialwise matrices without retaining a full
#' `T x P` response matrix. The shared trial design and error model are applied
#' to every block, and the resulting q-space coefficients and per-feature
#' residual statistics are concatenated into one [dkge_subject()] record. The
#' returned object retains one `T x q` design and one `q x P` beta matrix.
#'
#' @param chunks Either a list of numeric `T x P_block` matrices or a function
#'   called as `chunks(i)` that returns the next matrix and returns `NULL` after
#'   the last block. A function source is the memory-bounded path.
#' @inheritParams dkge_trial_subject
#' @param omega Optional final voxel/parcel weighting of length total `P`, or a
#'   total-`P` square matrix. It is applied only after all chunks are reduced.
#' @param max_chunks Safety limit for function sources that fail to terminate.
#'
#' @return A `dkge_subject` with the sufficient statistics needed by fitting.
#'   The full trialwise response and the reconstructible q-by-P effect score
#'   (`effect_information %*% beta`) are not stored.
#' @export
#' @examples
#' set.seed(1)
#' X <- model.matrix(~ 0 + factor(rep(1:3, each = 4)))
#' blocks <- list(matrix(rnorm(12 * 2), 12, 2),
#'                matrix(rnorm(12 * 3), 12, 3))
#' subject <- dkge_trial_subject_chunks(blocks, X)
dkge_trial_subject_chunks <- function(
    chunks,
    design,
    id = NULL,
    omega = NULL,
    effect_precision = NULL,
    trial_covariance = NULL,
    trial_precision = NULL,
    effect_noise_cov = NULL,
    split = c("none", "within_cell", "alternate", "run", "explicit"),
    split_labels = NULL,
    run_labels = NULL,
    split_independent = FALSE,
    tol = 1e-8,
    max_chunks = 100000L) {
  split_was_missing <- missing(split)
  split <- match.arg(split)
  split <- .dkge_resolve_trial_split(
    split, split_was_missing, split_labels, run_labels
  )
  design <- as.matrix(design)
  if (!is.list(chunks) && !is.function(chunks)) {
    stop("`chunks` must be a list or a function source.", call. = FALSE)
  }
  if (!is.numeric(max_chunks) || length(max_chunks) != 1L ||
      !is.finite(max_chunks) || max_chunks < 1) {
    stop("`max_chunks` must be a positive finite scalar.", call. = FALSE)
  }
  split_precision <- is.character(effect_precision) &&
    identical(effect_precision, "split_half")
  block_precision <- if (split_precision) NULL else effect_precision

  records <- list()
  reduce_chunk <- function(y, i) {
    if (is.list(y) && !is.null(y$y)) y <- y$y
    y <- as.matrix(y)
    if (!is.numeric(y) || nrow(y) != nrow(design) || !ncol(y)) {
      stop(sprintf("Chunk %d must be a numeric T x P_block matrix.", i),
           call. = FALSE)
    }
    if (is.null(colnames(y))) {
      prior_width <- sum(vapply(records, function(record) ncol(record$beta),
                                integer(1)))
      colnames(y) <- paste0("feature", prior_width + seq_len(ncol(y)))
    }
    record <- withCallingHandlers(
      dkge_trial_subject(
        y = y,
        design = design,
        id = id,
        effect_precision = block_precision,
        trial_covariance = trial_covariance,
        trial_precision = trial_precision,
        effect_noise_cov = effect_noise_cov,
        split = split,
        split_labels = split_labels,
        run_labels = run_labels,
        split_independent = split_independent,
        tol = tol
      ),
      warning = function(w) {
        if (grepl("beta matrix has reduced rank", conditionMessage(w),
                  fixed = TRUE)) {
          invokeRestart("muffleWarning")
        }
      }
    )
    # The score is exactly effect_information %*% beta. Retaining it would
    # duplicate the dominant q-by-P object and defeat the memory benefit of
    # chunking at imaging scale.
    record$effect_score <- NULL
    records[[length(records) + 1L]] <<- record
  }

  if (is.list(chunks)) {
    if (!length(chunks)) {
      stop("`chunks` produced no feature blocks.", call. = FALSE)
    }
    for (i in seq_along(chunks)) {
      if (is.null(chunks[[i]])) {
        stop(sprintf(
          "Chunk %d is NULL; list chunks must be numeric matrices.", i
        ), call. = FALSE)
      }
      reduce_chunk(chunks[[i]], i)
    }
  } else {
    for (i in seq_len(as.integer(max_chunks))) {
      y <- chunks(i)
      if (is.null(y)) break
      reduce_chunk(y, i)
    }
    if (!length(records)) {
      stop("`chunks` produced no feature blocks.", call. = FALSE)
    }
    if (length(records) == as.integer(max_chunks)) {
      stop("Chunk source reached `max_chunks` without returning NULL.",
           call. = FALSE)
    }
  }

  cluster_ids <- unlist(lapply(records, `[[`, "cluster_ids"), use.names = FALSE)
  if (anyDuplicated(cluster_ids)) {
    stop("Feature names must be unique across chunks.", call. = FALSE)
  }
  first <- records[[1L]]
  beta <- do.call(cbind, lapply(records, `[[`, "beta"))
  residual_variance <- unlist(
    lapply(records, `[[`, "residual_variance"), use.names = TRUE
  )
  residual_sum_squares <- unlist(
    lapply(records, `[[`, "residual_sum_squares"), use.names = TRUE
  )
  split_betas <- if (is.null(first$split_betas)) {
    NULL
  } else {
    lapply(seq_along(first$split_betas), function(half) {
      do.call(cbind, lapply(records, function(record) record$split_betas[[half]]))
    })
  }
  block_traces <- lapply(records, `[[`, "noise_trace")
  noise_trace <- if (all(vapply(block_traces, Negate(is.null), logical(1)))) {
    sum(unlist(block_traces, use.names = FALSE))
  } else if (all(is.finite(residual_variance))) {
    sum(residual_variance)
  } else {
    NULL
  }

  subject <- dkge_subject(
    beta,
    design = design,
    id = id,
    omega = omega,
    effect_n = first$effect_n,
    effect_precision = block_precision,
    effect_noise_cov = first$effect_noise_cov,
    residual_variance = residual_variance,
    residual_df = first$residual_df,
    noise_trace = noise_trace,
    split_betas = split_betas
  )
  subject$input_type <- "trialwise_chunks"
  subject$n_trials <- nrow(design)
  subject$effect_information <- first$effect_information
  subject$effect_score <- NULL
  subject$residual_sum_squares <- residual_sum_squares
  subject$noise_trace_scope <- if (
    all(vapply(records, function(record) {
      identical(record$noise_trace_scope, "supplied")
    }, logical(1)))
  ) {
    "supplied"
  } else {
    "unweighted"
  }
  subject$split_provenance <- first$split_provenance
  subject$error_model <- first$error_model
  subject$error_model$storage <- "feature_chunked"
  subject$error_model$n_chunks <- length(records)
  subject$error_model$effect_score_retained <- FALSE
  if (split_precision) {
    subject$effect_precision <- dkge_split_effect_precision(subject)
    subject$split_reliability <- attr(subject$effect_precision, "diagnostics")
  }
  subject
}

#' Resolve implicit split modes shared by dense and chunked constructors
#'
#' @keywords internal
#' @noRd
.dkge_resolve_trial_split <- function(split, split_was_missing,
                                      split_labels, run_labels) {
  if (!is.null(split_labels)) {
    if (!split_was_missing && !split %in% c("none", "explicit")) {
      stop("`split_labels` cannot be combined with another split mode.",
           call. = FALSE)
    }
    return("explicit")
  }
  if (!is.null(run_labels) && isTRUE(split_was_missing)) {
    return("run")
  }
  split
}

#' Validate and construct a trial-space precision action
#'
#' @keywords internal
#' @noRd
.dkge_trial_weighting <- function(trial_covariance, trial_precision,
                                  n_trials, tol) {
  if (!is.null(trial_covariance) && !is.null(trial_precision)) {
    stop("Supply at most one of `trial_covariance` and `trial_precision`.",
         call. = FALSE)
  }
  if (is.null(trial_covariance) && is.null(trial_precision)) {
    return(list(type = "iid", apply = function(z) z,
                covariance = NULL, precision = NULL))
  }

  validate_spd <- function(x, argument) {
    x <- as.matrix(x)
    if (!identical(dim(x), c(n_trials, n_trials))) {
      stop(sprintf("`%s` must be a T x T matrix.", argument), call. = FALSE)
    }
    if (!is.numeric(x) || any(!is.finite(x))) {
      stop(sprintf("`%s` must be a finite numeric matrix.", argument),
           call. = FALSE)
    }
    scale <- max(1, max(abs(x)))
    if (max(abs(x - t(x))) > tol * scale) {
      stop(sprintf("`%s` must be symmetric.", argument), call. = FALSE)
    }
    x <- (x + t(x)) / 2
    ch <- tryCatch(chol(x), error = function(e) NULL)
    if (is.null(ch)) {
      stop(sprintf("`%s` must be positive definite.", argument),
           call. = FALSE)
    }
    list(matrix = x, chol = ch)
  }

  if (!is.null(trial_covariance)) {
    checked <- validate_spd(trial_covariance, "trial_covariance")
    apply_precision <- function(z) {
      z <- as.matrix(z)
      backsolve(checked$chol, forwardsolve(t(checked$chol), z))
    }
    return(list(type = "covariance", apply = apply_precision,
                covariance = checked$matrix, precision = NULL))
  }

  precision_type <- "precision"
  if (is.function(trial_precision)) {
    candidate <- tryCatch(
      trial_precision(diag(n_trials)),
      error = function(e) {
        stop(sprintf("`trial_precision` operator failed: %s", conditionMessage(e)),
             call. = FALSE)
      }
    )
    checked <- validate_spd(candidate, "trial_precision operator result")
    precision_type <- "precision_operator"
  } else {
    checked <- validate_spd(trial_precision, "trial_precision")
  }
  W <- checked$matrix
  list(type = precision_type, apply = function(z) W %*% as.matrix(z),
       covariance = NULL, precision = W)
}

#' Restrict a trial-error model to one split half
#'
#' Precision inputs are converted through their implied covariance before
#' subsetting, so the half receives the marginal rather than conditional error
#' covariance when the original precision is not block diagonal.
#'
#' @keywords internal
#' @noRd
.dkge_subset_trial_weighting <- function(weighting, indices, tol) {
  if (identical(weighting$type, "iid")) {
    return(list(type = "iid", apply = function(z) z,
                covariance = NULL, precision = NULL))
  }
  covariance <- weighting$covariance
  if (is.null(covariance)) {
    precision_chol <- chol(weighting$precision)
    covariance <- chol2inv(precision_chol)
  }
  covariance_half <- covariance[indices, indices, drop = FALSE]
  .dkge_trial_weighting(covariance_half, NULL, length(indices), tol)
}

#' Verify covariance independence of declared split halves
#'
#' @keywords internal
#' @noRd
.dkge_assert_split_covariance_independence <- function(weighting, half_indices,
                                                        tol) {
  covariance <- weighting$covariance
  if (is.null(covariance)) {
    covariance <- chol2inv(chol(weighting$precision))
  }
  cross_block <- covariance[half_indices$first, half_indices$second,
                            drop = FALSE]
  scale <- max(1, max(abs(covariance)))
  if (length(cross_block) && max(abs(cross_block)) > tol * scale) {
    stop(paste0(
      "The supplied trial error model has non-zero covariance across halves ",
      "declared independent. Use an exploratory split or provide a block-independent model."
    ), call. = FALSE)
  }
  invisible(TRUE)
}

#' Export split-half reliability as effect precision
#'
#' Converts the across-spatial-unit agreement between two stored effect
#' estimates into a finite precision in `[0, 1]`. Negative or undefined
#' correlations map to zero. A deterministic count factor shrinks effects with
#' few observations in either half.
#'
#' @param subject A `dkge_subject` containing two `split_betas` matrices.
#' @param method Reliability transform. `"spearman_brown"` applies the
#'   Spearman-Brown correction to the non-negative half correlation;
#'   `"positive_r2"` squares the non-negative correlation.
#' @param min_features Minimum number of spatial units required to estimate a
#'   correlation.
#' @param count_prior Non-negative pseudo-count controlling low-count shrinkage.
#'   Zero disables count shrinkage.
#' @return A named numeric effect-precision vector in `[0, 1]`. Diagnostics are
#'   stored in the `"diagnostics"` attribute.
#' @export
dkge_split_effect_precision <- function(
    subject,
    method = c("spearman_brown", "positive_r2"),
    min_features = 3L,
    count_prior = 2) {
  if (!inherits(subject, "dkge_subject")) {
    stop("`subject` must inherit from 'dkge_subject'.", call. = FALSE)
  }
  method <- match.arg(method)
  min_features <- as.integer(min_features)
  if (length(min_features) != 1L || is.na(min_features) || min_features < 2L) {
    stop("`min_features` must be an integer of at least 2.", call. = FALSE)
  }
  if (!is.numeric(count_prior) || length(count_prior) != 1L ||
      !is.finite(count_prior) || count_prior < 0) {
    stop("`count_prior` must be a finite non-negative scalar.", call. = FALSE)
  }
  halves <- subject$split_betas
  if (is.null(halves) || !is.list(halves) || length(halves) != 2L) {
    stop("Split effect precision requires two stored split beta matrices.",
         call. = FALSE)
  }
  B1 <- as.matrix(halves[[1]])
  B2 <- as.matrix(halves[[2]])
  if (!identical(dim(B1), dim(B2))) {
    stop("Stored split beta matrices must have identical dimensions.",
         call. = FALSE)
  }
  q <- nrow(B1)
  effects <- rownames(B1) %||% paste0("effect_", seq_len(q))
  correlation <- numeric(q)
  for (j in seq_len(q)) {
    x <- as.numeric(B1[j, ])
    y <- as.numeric(B2[j, ])
    finite <- is.finite(x) & is.finite(y)
    if (sum(finite) < min_features ||
        stats::sd(x[finite]) <= sqrt(.Machine$double.eps) ||
        stats::sd(y[finite]) <= sqrt(.Machine$double.eps)) {
      correlation[[j]] <- 0
    } else {
      value <- stats::cor(x[finite], y[finite])
      correlation[[j]] <- if (is.finite(value)) value else 0
    }
  }
  positive <- pmax(correlation, 0)
  reliability <- if (identical(method, "spearman_brown")) {
    2 * positive / (1 + positive)
  } else {
    positive^2
  }

  counts <- subject$split_provenance$effect_counts %||% NULL
  min_count <- if (is.null(counts)) {
    rep(Inf, q)
  } else {
    if (!identical(dim(counts), c(q, 2L))) {
      stop("Split provenance has incompatible effect counts.", call. = FALSE)
    }
    apply(counts, 1L, min)
  }
  count_factor <- if (count_prior == 0) {
    rep(1, q)
  } else {
    min_count / (min_count + count_prior)
  }
  count_factor[!is.finite(count_factor)] <- 1
  precision <- pmin(1, pmax(0, reliability * count_factor))
  names(precision) <- effects
  attr(precision, "diagnostics") <- list(
    method = method,
    correlation = stats::setNames(correlation, effects),
    reliability = stats::setNames(reliability, effects),
    min_half_count = stats::setNames(min_count, effects),
    count_factor = stats::setNames(count_factor, effects),
    independent = isTRUE(subject$split_provenance$independent),
    independence_basis = subject$split_provenance$independence_basis %||% NA_character_
  )
  precision
}
