# dkge-between-infer.R
# Formula-aware inference for between-subject DKGE models.

#' Resampling tests for between-subject DKGE RRR terms
#'
#' Performs term-specific global tests for a fitted [dkge_between_rrr()] model.
#' Residual-space rotation fixes the nuisance-design column space and rotates
#' its orthogonal complement before comparing the reduced and full reduced-rank
#' residual sums of squares. Freedman-Lane residual permutation remains
#' available as a legacy approximate method.
#'
#' @param object A `dkge_between_rrr` object.
#' @param terms Character vector of terms or model-matrix columns to test. When
#'   `NULL`, all non-intercept formula terms are tested.
#' @param method Resampling method. `"rotation"` uses Haar rotations in the
#'   orthogonal complement of the reduced design; `"freedman_lane"` permutes
#'   reduced-model residual rows.
#' @param B Number of rotations or permutations.
#' @param blocks Optional exchangeability blocks of length `n_subjects`.
#' @param seed Optional random seed.
#' @param adjust P-value adjustment method passed to [stats::p.adjust()].
#' @param statistic Test statistic. Currently `"frob"` for reduced-minus-full
#'   weighted residual sum of squares.
#' @param scope `"global"` for the term-level test only, `"features"` for
#'   featurewise tests only, or `"both"`.
#' @param feature_adjust Multiplicity adjustment for featurewise p-values.
#' @param parallel Logical; evaluate resampling replicates with
#'   [future.apply::future_lapply()]. Randomization descriptors are always
#'   drawn serially up front, so results are identical to the serial path for a
#'   given `seed`. The serial path streams replicates and is preferred when
#'   `scope != "global"` with many features.
#'
#' @details
#' For `method = "rotation"`, let `Q0` span the reduced design and `Qe` span
#' its orthogonal complement. Each null response is equivalent to
#' `Q0 Q0' Y + Qe G Qe' Y`, where `G` is Haar-uniform on the residual-space
#' orthogonal group. This fixes every reduced-model mean and preserves
#' `crossprod(Y)` exactly. The test is finite-sample exact under a matrix-normal
#' null with independent homoscedastic subject rows and arbitrary feature
#' covariance. Version 1 supports unweighted fits and one exchangeability block
#' only; unsupported weights or blocks are rejected. With non-Gaussian or
#' heteroscedastic subject errors, rotation is approximate rather than exact.
#'
#' Both methods refit the reduced model at the same requested `rank` as the
#' full model.
#' When the tested term removes enough design columns that the reduced model's
#' estimable rank drops below that value, the reduced fit is clipped to its own
#' maximum rank while the full fit is not, and the reduced-minus-full residual
#' statistic then also absorbs that rank difference. This affects the observed
#' statistic and every permutation replicate alike, but the null is not exactly
#' the null of the rank-matched comparison. Moreover, Freedman-Lane substitutes
#' reduced-model residuals for unobservable errors and is not an exact finite-
#' sample procedure for a general reduced model.
#'
#' The legacy Freedman-Lane method's frozen package audit
#' (`null-calibration-v1`) found empirical
#' type-I error of 0.074, 0.076, and 0.077 at nominal 0.05 for three terms with
#' 20 subjects (1,000 null data sets per term and 399 permutations). Distortion
#' attenuated but was mixed at 40 subjects; at 80 subjects the three estimates
#' were 0.030, 0.048, and 0.052. These results establish a material small-sample
#' limitation, not merely uncertainty from the number of permutations. Treat
#' small-sample p-values close to alpha as exploratory, report the limitation,
#' and use an independently justified exchangeability design. Increasing `B`
#' improves Monte Carlo resolution but does not remove this size distortion.
#' The frozen plan, replicate p-values, and summary are shipped under
#' `data-raw/dkge-null-calibration-*` and `inst/extdata/dkge-null-calibration-*`.
#'
#' A separate predeclared audit of residual rotation found global-null sizes of
#' 0.049, 0.055, and 0.058 and nonzero-nuisance null sizes of 0.040, 0.042,
#' and 0.066 for the same three terms at 20 subjects. One nuisance-interaction
#' result missed the frozen promotion ceiling by one rejection, and rotation's
#' strong-interaction power was 0.430 versus 0.580 for Freedman-Lane. That raw
#' gap is mostly Freedman-Lane size inflation (matched-condition sizes 0.050
#' versus 0.086; size-matched Freedman-Lane power about 0.49). The predeclared
#' gates are not relaxed, so rotation remains opt-in. Its frozen plan, 8,100
#' results, and decision report are shipped as
#' `data-raw/dkge-between-rotation-*` and
#' `inst/extdata/dkge-between-rotation-*`.
#'
#' @return Object of class `dkge_between_permutation`.
#' @export
dkge_between_permute <- function(object,
                                 terms = NULL,
                                 method = c("freedman_lane", "rotation"),
                                 B = 999L,
                                 blocks = NULL,
                                 seed = NULL,
                                 adjust = "none",
                                 statistic = c("frob"),
                                 scope = c("global", "features", "both"),
                                 feature_adjust = c("none", "fdr", "maxT"),
                                 parallel = FALSE) {
  stopifnot(inherits(object, "dkge_between_rrr"))
  method <- match.arg(method)
  statistic <- match.arg(statistic)
  scope <- match.arg(scope)
  feature_adjust <- match.arg(feature_adjust)
  B <- .dkge_validate_resample_B(B)

  target <- object$target
  design <- object$design
  mask <- object$target_mask
  Y <- target$Y[, mask, drop = FALSE]
  X <- design$X[match(target$subject_ids, design$subject_ids), , drop = FALSE]

  terms <- .dkge_between_terms_to_test(design, terms)
  if (!length(terms)) {
    stop("No terms available to test.", call. = FALSE)
  }

  if (is.null(blocks)) {
    blocks <- rep(1L, nrow(Y))
  } else {
    blocks <- as.factor(blocks)
    if (length(blocks) != nrow(Y)) {
      stop("`blocks` must have length equal to the number of subjects.", call. = FALSE)
    }
    if (anyNA(blocks)) {
      stop("`blocks` cannot contain missing values.", call. = FALSE)
    }
  }

  if (identical(method, "rotation")) {
    .dkge_validate_rotation_scope(object, blocks)
  }

  seed_state <- .dkge_seed_enter(seed)
  on.exit(.dkge_seed_exit(seed_state), add = TRUE)

  # Draw once, before any term evaluation (including parallel
  # future.seed advances). Every term then sees the same descriptors, so
  # a term's null does not depend on which other terms were requested.
  resample_descriptors <- .dkge_between_resample_descriptors(method, B, blocks)

  tests <- lapply(terms, function(term) {
    .dkge_between_test_term(object = object,
                            term = term,
                            method = method,
                            X = X,
                            Y = Y,
                            B = B,
                            blocks = blocks,
                            statistic = statistic,
                            scope = scope,
                            feature_adjust = feature_adjust,
                            parallel = parallel,
                            resample_descriptors = resample_descriptors)
  })
  names(tests) <- terms

  summary <- if (!identical(scope, "features")) {
    data.frame(
      term = terms,
      statistic = vapply(tests, `[[`, numeric(1), "statistic"),
      p = vapply(tests, `[[`, numeric(1), "p"),
      B = B,
      method = method,
      statistic_name = statistic,
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      term = terms,
      statistic = NA_real_,
      p = NA_real_,
      B = B,
      method = method,
      statistic_name = statistic,
      stringsAsFactors = FALSE
    )
  }
  summary$p_adjusted <- if (!identical(scope, "features")) {
    stats::p.adjust(summary$p, method = adjust)
  } else {
    NA_real_
  }

  feature_tests <- if (identical(scope, "global")) {
    NULL
  } else {
    stats::setNames(lapply(tests, `[[`, "feature_test"), terms)
  }
  term_maps <- stats::setNames(lapply(tests, `[[`, "term_map"), terms)

  out <- list(
    summary = summary,
    tests = tests,
    feature_tests = feature_tests,
    term_maps = term_maps,
    method = method,
    statistic = statistic,
    B = B,
    adjust = adjust,
    scope = scope,
    feature_adjust = feature_adjust,
    terms = terms,
    resampling_assumptions = .dkge_between_resampling_assumptions(method),
    call = match.call()
  )
  class(out) <- c("dkge_between_permutation", "list")
  out
}

.dkge_between_terms_to_test <- function(design, terms) {
  if (!is.null(terms)) return(as.character(terms))
  candidates <- setdiff(design$term_names, "(Intercept)")
  if (!length(candidates)) {
    candidates <- setdiff(colnames(design$X), "(Intercept)")
  }
  candidates
}

.dkge_between_term_columns <- function(design, term) {
  if (!is.null(design$term_columns[[term]])) {
    return(design$term_columns[[term]])
  }
  idx <- match(term, colnames(design$X))
  if (is.na(idx)) {
    stop("Unknown term or model-matrix column: ", term, call. = FALSE)
  }
  idx
}

.dkge_between_reduced_X <- function(design, X, term) {
  drop_cols <- .dkge_between_term_columns(design, term)
  keep <- setdiff(seq_len(ncol(X)), drop_cols)
  if (!length(keep)) {
    stop("Testing term '", term, "' would remove all design columns.", call. = FALSE)
  }
  X[, keep, drop = FALSE]
}

.dkge_between_resample_descriptors <- function(method, B, blocks) {
  if (identical(method, "freedman_lane")) {
    lapply(seq_len(B), function(b) .dkge_permute_within_blocks(blocks))
  } else {
    as.list(sample.int(.Machine$integer.max, B, replace = TRUE))
  }
}

.dkge_between_test_term <- function(object,
                                    term,
                                    method,
                                    X,
                                    Y,
                                    B,
                                    blocks,
                                    statistic,
                                    scope,
                                    feature_adjust,
                                    parallel = FALSE,
                                    resample_descriptors = NULL) {
  Xred <- .dkge_between_reduced_X(object$design, X, term)
  # The reduced model is refit under the *same* rank tolerance the full model
  # was accepted with. Falling back to the `.dkge_rrr_fit_core()` default would
  # reject a near-collinear design that `dkge_between_rrr(tol = ...)` fit
  # happily, so `tol` is threaded through every refit and QR setup below.
  tol <- object$tol %||% 1e-8
  red_obs <- .dkge_rrr_fit_core(X = Xred,
                                Y = Y,
                                rank = object$rank,
                                subject_weights = object$subject_weights,
                                feature_weights = object$feature_weights,
                                coef_ids = colnames(Xred),
                                feature_ids = colnames(Y),
                                subject_ids = rownames(Y),
                                tol = tol,
                                warn_rank = FALSE)
  stat_obs <- .dkge_between_stat(object, red_obs,
                                 subject_weights = object$subject_weights,
                                 feature_weights = object$feature_weights,
                                 statistic = statistic)
  term_map_obs <- .dkge_between_term_map_from_coef(object$coef,
                                                   design = object$design,
                                                   term = term,
                                                   drop = FALSE)
  feature_stat_obs <- .dkge_between_feature_stat(term_map_obs)
  feature_ge <- integer(length(feature_stat_obs))
  null_max <- numeric(B)

  null <- numeric(B)
  full_setup <- .dkge_between_rrr_fast_setup(
    X,
    subject_weights = object$subject_weights,
    feature_weights = object$feature_weights,
    tol = tol
  )
  red_setup <- .dkge_between_rrr_fast_setup(
    Xred,
    subject_weights = object$subject_weights,
    feature_weights = object$feature_weights,
    tol = tol
  )
  term_rows <- .dkge_between_term_rows(object$design, object$coef, term)
  global_only <- identical(scope, "global")
  # Coefficients are only needed for the featurewise statistic; under
  # scope = "global" the per-permutation backsolve is pure overhead.
  perm_coef_rows <- if (global_only) NULL else term_rows

  if (is.null(resample_descriptors)) {
    resample_descriptors <- .dkge_between_resample_descriptors(method, B, blocks)
  }
  if (length(resample_descriptors) != B) {
    stop("`resample_descriptors` must have length B.", call. = FALSE)
  }

  if (identical(method, "freedman_lane")) {
    resample_setup <- .dkge_between_permutation_setup(
      full_setup, red_setup, red_obs
    )
    make_input <- function(descriptor) {
      .dkge_between_permutation_eval_input(resample_setup, descriptor)
    }
  } else {
    resample_setup <- .dkge_between_rotation_setup(
      full_setup, red_setup, Y, tol = tol
    )
    make_input <- function(descriptor) {
      rotation <- .dkge_seeded_haar_rotation(
        descriptor, resample_setup$residual_dimension
      )
      .dkge_between_rotation_eval_input(resample_setup, rotation)
    }
  }

  eval_one <- function(descriptor) {
    eval_input <- make_input(descriptor)
    full_b <- .dkge_between_rrr_fast_eval_from_C(
      full_setup,
      eval_input$C_full,
      yss = eval_input$yss,
      rank = object$rank,
      coef_rows = perm_coef_rows,
      coef_ids = colnames(X),
      feature_ids = colnames(Y)
    )
    red_b <- .dkge_between_rrr_fast_eval_from_C(
      red_setup,
      eval_input$C_red,
      yss = eval_input$yss,
      rank = object$rank,
      coef_ids = colnames(Xred),
      feature_ids = colnames(Y)
    )
    list(stat = .dkge_between_fast_stat(full_b, red_b, statistic = statistic),
         feature_stat = if (global_only) NULL else .dkge_between_feature_stat(full_b$coef))
  }

  # Only materialise every replicate when running in parallel; the serial path
  # streams so that featurewise nulls stay O(features) rather than O(B x features).
  resample_results <- if (isTRUE(parallel)) {
    .dkge_apply(resample_descriptors, eval_one, parallel = TRUE,
                future.seed = TRUE)
  } else {
    NULL
  }

  feature_tol <- sqrt(.Machine$double.eps) * pmax(1, abs(feature_stat_obs))
  for (b in seq_len(B)) {
    res <- if (is.null(resample_results)) {
      eval_one(resample_descriptors[[b]])
    } else {
      resample_results[[b]]
    }
    null[b] <- res$stat
    if (!global_only) {
      feature_stat_b <- res$feature_stat
      feature_ge <- feature_ge + as.integer(feature_stat_b >= feature_stat_obs - feature_tol)
      null_max[b] <- max(feature_stat_b)
    }
  }
  stat_tol <- sqrt(.Machine$double.eps) * max(1, abs(stat_obs))
  p <- (1 + sum(null >= stat_obs - stat_tol)) / (B + 1)
  feature_test <- if (identical(scope, "global")) {
    NULL
  } else {
    p_raw <- (1 + feature_ge) / (B + 1)
    p_adj <- switch(
      feature_adjust,
      none = p_raw,
      fdr = stats::p.adjust(p_raw, method = "fdr"),
      maxT = vapply(feature_stat_obs,
                    function(s) {
                      tol <- sqrt(.Machine$double.eps) * max(1, abs(s))
                      (1 + sum(null_max >= s - tol)) / (B + 1)
                    },
                    numeric(1))
    )
    list(
      term = term,
      term_map = term_map_obs,
      statistic = feature_stat_obs,
      p = p_raw,
      p_adjusted = p_adj,
      adjust = feature_adjust,
      null_max = null_max,
      feature_ids = colnames(Y)
    )
  }

  list(term = term,
       statistic = stat_obs,
       p = p,
       null = null,
       term_map = term_map_obs,
       feature_test = feature_test)
}

#' Validate the supported residual-rotation scope
#'
#' @keywords internal
#' @noRd
.dkge_validate_rotation_scope <- function(object, blocks) {
  if (!is.null(object$subject_weights) || !is.null(object$feature_weights)) {
    stop(paste0(
      "method='rotation' currently requires an unweighted between-subject fit; ",
      "use method='freedman_lane' explicitly or refit with weights='none'."
    ), call. = FALSE)
  }
  if (nlevels(as.factor(blocks)) != 1L) {
    stop(paste0(
      "method='rotation' currently supports one exchangeability block only; ",
      "block-restricted rotations are not yet implemented."
    ), call. = FALSE)
  }
  invisible(TRUE)
}

#' Declare method-specific resampling assumptions
#'
#' @keywords internal
#' @noRd
.dkge_between_resampling_assumptions <- function(method) {
  if (identical(method, "rotation")) {
    list(
      exact_under = paste(
        "matrix-normal null with independent homoscedastic subject rows",
        "and arbitrary feature covariance"
      ),
      weights = "none",
      exchangeability_blocks = "single",
      non_gaussian = "approximate"
    )
  } else {
    list(
      exact_under = "exchangeable reduced-model errors with known nuisance effects",
      nuisance_residuals = "approximate in finite samples"
    )
  }
}

#' Construct compressed residual-rotation sufficient statistics
#'
#' @keywords internal
#' @noRd
.dkge_between_rotation_setup <- function(full_setup, red_setup, Y, tol = 1e-8) {
  Y <- as.matrix(Y)
  n <- nrow(Y)
  qr0 <- qr(red_setup$X, tol = tol)
  if (qr0$rank != ncol(red_setup$X)) {
    stop("Reduced design is rank deficient for residual-space rotation.",
         call. = FALSE)
  }
  residual_dimension <- n - qr0$rank
  if (residual_dimension < 2L) {
    stop(paste0(
      "method='rotation' requires at least two residual dimensions after ",
      "removing the tested term."
    ), call. = FALSE)
  }

  Q0 <- qr.Q(qr0, complete = FALSE)
  # Apply the compact Householder representation to only the trailing identity
  # columns.  qr.Q(qr0, complete = TRUE) would first materialise an n x n
  # subject matrix merely to retain these n x residual_dimension columns.
  complement_seed <- matrix(0, n, residual_dimension)
  complement_seed[cbind(qr0$rank + seq_len(residual_dimension),
                        seq_len(residual_dimension))] <- 1
  Qe <- qr.qy(qr0, complement_seed)
  fitted_coordinates <- crossprod(Q0, Y)
  residual_coordinates <- crossprod(Qe, Y)
  yss <- sum(fitted_coordinates * fitted_coordinates) +
    sum(residual_coordinates * residual_coordinates)
  yss_direct <- sum(Y * Y)
  scale <- max(1, yss_direct)
  if (abs(yss - yss_direct) > 100 * tol * scale) {
    stop("Residual-space rotation decomposition failed to preserve response energy.",
         call. = FALSE)
  }

  list(
    Q0 = Q0,
    Qe = Qe,
    fitted_coordinates = fitted_coordinates,
    residual_coordinates = residual_coordinates,
    full_fitted_map = crossprod(full_setup$Q, Q0),
    full_residual_map = crossprod(full_setup$Q, Qe),
    red_fitted_map = crossprod(red_setup$Q, Q0),
    red_residual_map = crossprod(red_setup$Q, Qe),
    yss = yss_direct,
    residual_dimension = residual_dimension,
    n = n,
    feature_ids = colnames(Y)
  )
}

#' Evaluate one residual-space rotation in compressed coordinates
#'
#' @keywords internal
#' @noRd
.dkge_between_rotation_eval_input <- function(rotation_setup, rotation) {
  rotation <- as.matrix(rotation)
  d <- rotation_setup$residual_dimension
  if (!identical(dim(rotation), c(d, d)) || !is.numeric(rotation) ||
      any(!is.finite(rotation))) {
    stop("`rotation` must be a finite residual-dimension square matrix.",
         call. = FALSE)
  }
  orthogonal_error <- max(abs(crossprod(rotation) - diag(d)))
  if (orthogonal_error > 1e-10) {
    stop("`rotation` must be orthogonal.", call. = FALSE)
  }
  rotated_residual <- rotation %*% rotation_setup$residual_coordinates
  list(
    C_full = rotation_setup$full_fitted_map %*%
      rotation_setup$fitted_coordinates +
      rotation_setup$full_residual_map %*% rotated_residual,
    C_red = rotation_setup$red_fitted_map %*%
      rotation_setup$fitted_coordinates +
      rotation_setup$red_residual_map %*% rotated_residual,
    yss = rotation_setup$yss
  )
}

#' Draw a Haar-uniform orthogonal matrix
#'
#' @keywords internal
#' @noRd
.dkge_haar_orthogonal <- function(d) {
  d <- as.integer(d)
  if (length(d) != 1L || is.na(d) || d < 1L) {
    stop("`d` must be a positive integer.", call. = FALSE)
  }
  decomposition <- qr(matrix(stats::rnorm(d * d), d, d))
  Q <- qr.Q(decomposition)
  signs <- sign(diag(qr.R(decomposition)))
  signs[signs == 0] <- 1
  sweep(Q, 2L, signs, "*")
}

#' Draw a rotation from an isolated deterministic RNG scope
#'
#' @keywords internal
#' @noRd
.dkge_seeded_haar_rotation <- function(seed, d) {
  # future.apply workers normally use L'Ecuyer-CMRG whereas the serial caller
  # normally uses Mersenne-Twister.  Merely calling set.seed(seed) therefore
  # does not identify the same draw on the two backends.  Pin the complete RNG
  # algorithm for this isolated draw, then restore both the caller's kind and
  # seed exactly.
  old_kind <- RNGkind()
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv,
                         inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    } else {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }
  }, add = TRUE)
  RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion",
          sample.kind = "Rejection")
  set.seed(as.integer(seed))
  .dkge_haar_orthogonal(d)
}

.dkge_between_permutation_setup <- function(full_setup, red_setup, red_obs) {
  fitted <- as.matrix(red_obs$fitted)
  residuals <- as.matrix(red_obs$residuals)
  if (nrow(fitted) != nrow(residuals) || ncol(fitted) != ncol(residuals)) {
    stop("Reduced fitted values and residuals must have matching dimensions.",
         call. = FALSE)
  }

  n <- nrow(fitted)
  subject_sqrt <- full_setup$subject_sqrt %||% rep(1, n)
  feature_sqrt <- full_setup$feature_sqrt

  fitted_fw <- fitted
  residuals_fw <- residuals
  if (!is.null(feature_sqrt)) {
    if (length(feature_sqrt) != ncol(fitted_fw)) {
      stop("Feature weights must match target columns.", call. = FALSE)
    }
    fitted_fw <- sweep(fitted_fw, 2L, feature_sqrt, "*")
    residuals_fw <- sweep(residuals_fw, 2L, feature_sqrt, "*")
  }

  subject_sq <- subject_sqrt^2
  fitted_w <- fitted_fw * subject_sqrt
  Q_all <- cbind(full_setup$Q, red_setup$Q)
  Q_resid <- Q_all * subject_sqrt

  list(
    residuals_fw = residuals_fw,
    C_fit = crossprod(Q_all, fitted_w),
    Q_resid = Q_resid,
    fitted_ss = sum(fitted_w * fitted_w),
    residual_row_ss = rowSums(residuals_fw * residuals_fw),
    fitted_residual_cross = fitted_fw %*% t(residuals_fw),
    subject_sq = subject_sq,
    n = n,
    n_full = ncol(full_setup$Q),
    n_red = ncol(red_setup$Q)
  )
}

.dkge_between_permutation_eval_input <- function(perm_setup, idx) {
  idx <- as.integer(idx)
  if (length(idx) != perm_setup$n || anyNA(idx)) {
    stop("Permutation index has the wrong length.", call. = FALSE)
  }
  source_to_dest <- match(seq_len(perm_setup$n), idx)
  if (anyNA(source_to_dest)) {
    stop("Permutation index must contain each subject exactly once.",
         call. = FALSE)
  }

  C_all <- perm_setup$C_fit +
    crossprod(perm_setup$Q_resid[source_to_dest, , drop = FALSE],
              perm_setup$residuals_fw)

  source <- seq_len(perm_setup$n)
  resid_ss <- sum(perm_setup$subject_sq[source_to_dest] *
                    perm_setup$residual_row_ss)
  cross_ss <- sum(perm_setup$subject_sq[source_to_dest] *
                    perm_setup$fitted_residual_cross[cbind(source_to_dest, source)])
  yss <- perm_setup$fitted_ss + resid_ss + 2 * cross_ss

  n_full <- perm_setup$n_full
  n_red <- perm_setup$n_red
  list(
    C_full = C_all[seq_len(n_full), , drop = FALSE],
    C_red = C_all[n_full + seq_len(n_red), , drop = FALSE],
    yss = yss
  )
}

.dkge_between_rrr_fast_setup <- function(X,
                                         subject_weights = NULL,
                                         feature_weights = NULL,
                                         tol = 1e-8) {
  X <- as.matrix(X)
  Xw <- X
  subject_sqrt <- NULL
  feature_sqrt <- NULL

  if (!is.null(subject_weights)) {
    subject_weights <- as.numeric(subject_weights)
    if (length(subject_weights) != nrow(X) ||
        any(!is.finite(subject_weights)) ||
        any(subject_weights < 0)) {
      stop("Subject weights must be non-negative, finite, and match input rows.",
           call. = FALSE)
    }
    subject_sqrt <- sqrt(pmax(subject_weights, 0))
    Xw <- Xw * subject_sqrt
  }
  if (!is.null(feature_weights)) {
    feature_weights <- as.numeric(feature_weights)
    if (any(!is.finite(feature_weights)) || any(feature_weights <= 0)) {
      stop("Feature weights must be positive and finite.", call. = FALSE)
    }
    feature_sqrt <- sqrt(feature_weights)
  }

  qrX <- qr(Xw, tol = tol)
  if (qrX$rank < ncol(Xw)) {
    stop("Weighted between-subject model matrix is rank deficient.", call. = FALSE)
  }

  list(
    X = X,
    Q = qr.Q(qrX),
    R = qr.R(qrX),
    subject_sqrt = subject_sqrt,
    feature_sqrt = feature_sqrt,
    tol = tol
  )
}

.dkge_between_rrr_fast_eval_from_C <- function(setup,
                                               C,
                                               yss,
                                               rank,
                                               coef_rows = NULL,
                                               coef_ids = colnames(setup$X),
                                               feature_ids = colnames(C)) {
  max_rank <- min(nrow(C), ncol(C))
  rank <- min(as.integer(rank), max_rank)

  gram <- tcrossprod(C)
  gram <- (gram + t(gram)) / 2
  eg <- eigen(gram, symmetric = TRUE)
  values <- pmax(eg$values, 0)
  keep <- seq_len(rank)
  U <- eg$vectors[, keep, drop = FALSE]
  explained <- sum(values[keep])
  sse <- max(0, yss - explained)

  coef <- NULL
  if (!is.null(coef_rows)) {
    projector <- U %*% t(U)
    coef_weighted <- backsolve(setup$R, projector)[coef_rows, , drop = FALSE] %*% C
    if (!is.null(setup$feature_sqrt)) {
      coef <- sweep(coef_weighted, 2L, setup$feature_sqrt, "/")
    } else {
      coef <- coef_weighted
    }
    rownames(coef) <- coef_ids[coef_rows]
    colnames(coef) <- feature_ids
  }

  list(explained = explained, sse = sse, coef = coef)
}

.dkge_between_fast_stat <- function(full, reduced, statistic = "frob") {
  if (!identical(statistic, "frob")) {
    stop("Unsupported statistic: ", statistic, call. = FALSE)
  }
  max(0, reduced$sse - full$sse)
}

.dkge_between_stat <- function(full,
                               reduced,
                               subject_weights = NULL,
                               feature_weights = NULL,
                               statistic = "frob") {
  if (!identical(statistic, "frob")) {
    stop("Unsupported statistic: ", statistic, call. = FALSE)
  }
  sse_red <- .dkge_between_weighted_sse(reduced$residuals, subject_weights, feature_weights)
  sse_full <- .dkge_between_weighted_sse(full$residuals, subject_weights, feature_weights)
  max(0, sse_red - sse_full)
}

.dkge_between_weighted_sse <- function(R, subject_weights = NULL, feature_weights = NULL) {
  R <- as.matrix(R)
  if (!is.null(subject_weights)) {
    R <- R * sqrt(pmax(as.numeric(subject_weights), 0))
  }
  if (!is.null(feature_weights)) {
    R <- sweep(R, 2L, sqrt(pmax(as.numeric(feature_weights), 0)), "*")
  }
  sum(R * R)
}

.dkge_between_feature_stat <- function(term_map) {
  term_map <- as.matrix(term_map)
  sqrt(colSums(term_map * term_map))
}

.dkge_permute_within_blocks <- function(blocks) {
  blocks <- as.factor(blocks)
  idx_all <- seq_along(blocks)
  out <- idx_all
  for (lvl in levels(blocks)) {
    idx <- which(blocks == lvl)
    out[idx] <- if (length(idx) <= 1L) idx else sample(idx, length(idx), replace = FALSE)
  }
  out
}

#' @export
print.dkge_between_permutation <- function(x, ...) {
  cat("<dkge_between_permutation>", "\n", sep = "")
  cat("  method       :", x$method, "\n")
  cat("  resamples    :", x$B, "\n")
  cat("  scope        :", x$scope, "\n")
  print(x$summary, row.names = FALSE)
  if (!is.null(x$feature_tests)) {
    cat("  feature tests:", paste(names(x$feature_tests), collapse = ", "), "\n")
  }
  invisible(x)
}
