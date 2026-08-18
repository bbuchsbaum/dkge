# dkge-trial-subject.R
# Trialwise beta input reduced to q-space sufficient statistics.

#' Build a DKGE subject from trialwise beta maps
#'
#' Fits the second-stage model `Y = X B + E` without constructing any
#' voxel-by-voxel covariance. The returned subject contains the q-by-P effect
#' estimates, `(X'X)^-1`, per-column residual variances, and (for a one-hot
#' design) per-effect trial counts. These are sufficient for count weighting and
#' analytic removal of the finite-trial noise moment in [dkge_fit()].
#'
#' @param y Numeric `T x P` matrix of trialwise beta maps.
#' @param design Numeric `T x q` second-stage design matrix. This may be a
#'   one-hot cell design or a general full-rank basis design.
#' @param id Optional subject identifier.
#' @param omega Optional voxel/parcel weighting passed to [dkge_subject()].
#' @param effect_precision Optional direct q-vector of effect precisions.
#' @param split Optional split-half sufficient statistic. `"within_cell"`
#'   alternates trials within each one-hot cell; `"alternate"` alternates all
#'   trial rows and is valid only when both half-designs remain full rank.
#' @param tol Numerical tolerance used for rank and one-hot checks.
#'
#' @return A `dkge_subject` carrying second-stage sufficient statistics.
#' @export
#' @examples
#' set.seed(1)
#' X <- model.matrix(~ 0 + factor(rep(1:3, each = 4)))
#' Y <- X %*% matrix(rnorm(3 * 5), 3, 5) + matrix(rnorm(12 * 5), 12, 5)
#' subject <- dkge_trial_subject(Y, X, id = "s1")
dkge_trial_subject <- function(y, design, id = NULL, omega = NULL,
                               effect_precision = NULL,
                               split = c("none", "within_cell", "alternate"),
                               tol = 1e-8) {
  split <- match.arg(split)
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

  effects <- colnames(design) %||% paste0("effect_", seq_len(ncol(design)))
  clusters <- colnames(y) %||% paste0("cluster_", seq_len(ncol(y)))
  colnames(design) <- effects
  colnames(y) <- clusters

  qrx <- qr(design, tol = tol)
  q <- ncol(design)
  if (qrx$rank < q) {
    stop(sprintf(
      "`design` must be full column rank for trialwise reduction (rank %d < %d).",
      qrx$rank, q
    ), call. = FALSE)
  }
  beta <- qr.coef(qrx, y)
  rownames(beta) <- effects
  colnames(beta) <- clusters
  residuals <- qr.resid(qrx, y)
  residual_df <- nrow(design) - qrx$rank
  residual_variance <- if (residual_df > 0) {
    colSums(residuals^2) / residual_df
  } else {
    rep(NA_real_, ncol(y))
  }
  names(residual_variance) <- clusters

  # Reuse the QR already computed above: X = QR gives X'X = R'R, so
  # (X'X)^{-1} = chol2inv(R). `qr()` may pivot, so undo the permutation.
  unpivot <- order(qrx$pivot)
  effect_noise_cov <- chol2inv(qr.R(qrx))[unpivot, unpivot, drop = FALSE]
  effect_noise_cov <- (effect_noise_cov + t(effect_noise_cov)) / 2
  dimnames(effect_noise_cov) <- list(effects, effects)

  one_hot <- all(abs(design - round(design)) <= tol) &&
    all(round(design) %in% c(0, 1)) &&
    all(rowSums(round(design)) == 1)
  effect_n <- if (one_hot) {
    stats::setNames(as.numeric(colSums(round(design))), effects)
  } else {
    NULL
  }

  split_betas <- NULL
  if (!identical(split, "none")) {
    if (identical(split, "within_cell") && !one_hot) {
      stop("split='within_cell' requires a one-hot design.", call. = FALSE)
    }
    split_group <- if (identical(split, "within_cell")) {
      cell <- max.col(round(design), ties.method = "first")
      stats::ave(seq_along(cell), cell,
                 FUN = function(idx) seq_along(idx) %% 2L)
    } else {
      seq_len(nrow(design)) %% 2L
    }
    half_indices <- list(first = which(split_group == 1L),
                         second = which(split_group == 0L))
    split_betas <- lapply(half_indices, function(idx) {
      Xh <- design[idx, , drop = FALSE]
      Yh <- y[idx, , drop = FALSE]
      qrh <- qr(Xh, tol = tol)
      if (qrh$rank < q) {
        stop(sprintf(
          "The '%s' split does not leave both half-designs full rank; add trials or choose split='none'.",
          split
        ), call. = FALSE)
      }
      Bh <- qr.coef(qrh, Yh)
      rownames(Bh) <- effects
      colnames(Bh) <- clusters
      Bh
    })
  }

  subject <- dkge_subject(
    beta,
    design = design,
    id = id,
    omega = omega,
    effect_n = effect_n,
    effect_precision = effect_precision,
    effect_noise_cov = effect_noise_cov,
    residual_variance = residual_variance,
    residual_df = residual_df,
    split_betas = split_betas
  )
  subject$input_type <- "trialwise"
  subject$n_trials <- nrow(design)
  subject
}

#' Build a trialwise DKGE subject from feature chunks
#'
#' Reduces successive `T x P_block` trialwise matrices without retaining a full
#' `T x P` response matrix. Each block is reduced with [dkge_trial_subject()],
#' then the q-space coefficients and per-feature residual statistics are joined
#' into one [dkge_subject()] record. The returned object retains one `T x q`
#' design and one `q x P` beta matrix; it never forms a `P x P` matrix.
#'
#' @param chunks Either a list of numeric `T x P_block` matrices or a function
#'   called as `chunks(i)` that returns the next matrix and returns `NULL` after
#'   the last block. A function source is the memory-bounded path.
#' @inheritParams dkge_trial_subject
#' @param omega Optional final voxel/parcel weighting of length total `P`, or a
#'   total-`P` square matrix. It is applied after all chunks are reduced.
#' @param max_chunks Safety limit for function sources that fail to terminate.
#'
#' @return A `dkge_subject` carrying the same OLS sufficient statistics as
#'   [dkge_trial_subject()] (`beta`, `effect_noise_cov`, `residual_variance`,
#'   `residual_df`, `effect_n`, optional `split_betas` / `omega`). Extra
#'   fields: `input_type = "trialwise_chunks"`, `n_trials`,
#'   `effect_information` (`X'X`), `residual_sum_squares`, and `error_model`
#'   (`storage = "feature_chunked"`, `n_chunks`, ...). `noise_trace` is not
#'   stored; [dkge_fit()] recomputes it from `residual_variance` so spatial
#'   and voxel weights are honoured. The full trialwise response is not stored.
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
    split = c("none", "within_cell", "alternate"),
    tol = 1e-8,
    max_chunks = 100000L) {
  split <- match.arg(split)
  design <- as.matrix(design)
  if (!is.list(chunks) && !is.function(chunks)) {
    stop("`chunks` must be a list or a function source.", call. = FALSE)
  }
  if (!is.numeric(max_chunks) || length(max_chunks) != 1L ||
      !is.finite(max_chunks) || max_chunks < 1 ||
      max_chunks > .Machine$integer.max) {
    stop("`max_chunks` must be a positive finite integer-sized scalar.",
         call. = FALSE)
  }
  max_chunks <- as.integer(max_chunks)
  if (is.list(chunks) && length(chunks) > max_chunks) {
    stop("The chunk list exceeds `max_chunks`.", call. = FALSE)
  }
  get_chunk <- if (is.function(chunks)) {
    chunks
  } else {
    function(i) if (i <= length(chunks)) chunks[[i]] else NULL
  }

  records <- list()
  for (i in seq_len(max_chunks)) {
    y <- get_chunk(i)
    if (is.null(y)) break
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
        effect_precision = effect_precision,
        split = split,
        tol = tol
      ),
      warning = function(w) {
        if (grepl("beta matrix has reduced rank", conditionMessage(w),
                  fixed = TRUE)) {
          invokeRestart("muffleWarning")
        }
      }
    )
    records[[length(records) + 1L]] <- record
  }
  if (!length(records)) {
    stop("`chunks` produced no feature blocks.", call. = FALSE)
  }
  if (is.function(chunks) && length(records) == max_chunks) {
    stop("Chunk source reached `max_chunks` without returning NULL.",
         call. = FALSE)
  }

  cluster_ids <- unlist(lapply(records, `[[`, "cluster_ids"),
                        use.names = FALSE)
  if (anyDuplicated(cluster_ids)) {
    stop("Feature names must be unique across chunks.", call. = FALSE)
  }
  first <- records[[1L]]
  beta <- do.call(cbind, lapply(records, `[[`, "beta"))
  residual_variance <- unlist(
    lapply(records, `[[`, "residual_variance"), use.names = TRUE
  )
  split_betas <- if (is.null(first$split_betas)) {
    NULL
  } else {
    lapply(seq_along(first$split_betas), function(half) {
      do.call(cbind, lapply(records, function(record) record$split_betas[[half]]))
    })
  }
  names(split_betas) <- names(first$split_betas)

  subject <- dkge_subject(
    beta,
    design = design,
    id = id,
    omega = omega,
    effect_n = first$effect_n,
    effect_precision = effect_precision,
    effect_noise_cov = first$effect_noise_cov,
    residual_variance = residual_variance,
    residual_df = first$residual_df,
    split_betas = split_betas
  )
  subject$input_type <- "trialwise_chunks"
  subject$n_trials <- nrow(design)
  subject$effect_information <- crossprod(design)
  dimnames(subject$effect_information) <- list(subject$effects, subject$effects)
  subject$effect_score <- NULL
  subject$residual_sum_squares <- if (is.null(subject$residual_df)) {
    NULL
  } else {
    residual_variance * subject$residual_df
  }
  subject$error_model <- list(
    estimator = "ols",
    trial_dependence = "iid",
    storage = "feature_chunked",
    n_chunks = length(records),
    effect_score_retained = FALSE,
    residual_df = subject$residual_df
  )
  subject
}
