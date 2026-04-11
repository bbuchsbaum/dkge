# dkge-between-infer.R
# Formula-aware inference for between-subject DKGE models.

#' Permutation tests for between-subject DKGE RRR terms
#'
#' Performs term-specific global tests for a fitted [dkge_between_rrr()] model.
#' The default uses Freedman-Lane residual permutations from a reduced model
#' that excludes the tested term, then compares the reduced and full
#' reduced-rank residual sums of squares.
#'
#' @param object A `dkge_between_rrr` object.
#' @param terms Character vector of terms or model-matrix columns to test. When
#'   `NULL`, all non-intercept formula terms are tested.
#' @param method Permutation method. Currently `"freedman_lane"`.
#' @param B Number of permutations.
#' @param blocks Optional exchangeability blocks of length `n_subjects`.
#' @param seed Optional random seed.
#' @param adjust P-value adjustment method passed to [stats::p.adjust()].
#' @param statistic Test statistic. Currently `"frob"` for reduced-minus-full
#'   weighted residual sum of squares.
#'
#' @return Object of class `dkge_between_permutation`.
#' @export
dkge_between_permute <- function(object,
                                 terms = NULL,
                                 method = c("freedman_lane"),
                                 B = 999L,
                                 blocks = NULL,
                                 seed = NULL,
                                 adjust = "none",
                                 statistic = c("frob"),
                                 scope = c("global", "features", "both"),
                                 feature_adjust = c("none", "fdr", "maxT")) {
  stopifnot(inherits(object, "dkge_between_rrr"))
  method <- match.arg(method)
  statistic <- match.arg(statistic)
  scope <- match.arg(scope)
  feature_adjust <- match.arg(feature_adjust)
  B <- as.integer(B)
  if (length(B) != 1L || is.na(B) || B < 1L) {
    stop("`B` must be a positive integer.", call. = FALSE)
  }

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
  }

  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    } else {
      NULL
    }
    on.exit({
      if (is.null(old_seed)) {
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      } else {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
  }

  tests <- lapply(terms, function(term) {
    .dkge_between_test_term(object = object,
                            term = term,
                            X = X,
                            Y = Y,
                            B = B,
                            blocks = blocks,
                            statistic = statistic,
                            scope = scope,
                            feature_adjust = feature_adjust)
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

.dkge_between_test_term <- function(object,
                                    term,
                                    X,
                                    Y,
                                    B,
                                    blocks,
                                    statistic,
                                    scope,
                                    feature_adjust) {
  Xred <- .dkge_between_reduced_X(object$design, X, term)
  full_obs <- .dkge_rrr_fit_core(X = X,
                                 Y = Y,
                                 rank = object$rank,
                                 subject_weights = object$subject_weights,
                                 feature_weights = object$feature_weights,
                                 coef_ids = colnames(X),
                                 feature_ids = colnames(Y),
                                 subject_ids = rownames(Y),
                                 warn_rank = FALSE)
  red_obs <- .dkge_rrr_fit_core(X = Xred,
                                Y = Y,
                                rank = object$rank,
                                subject_weights = object$subject_weights,
                                feature_weights = object$feature_weights,
                                coef_ids = colnames(Xred),
                                feature_ids = colnames(Y),
                                subject_ids = rownames(Y),
                                warn_rank = FALSE)
  stat_obs <- .dkge_between_stat(full_obs, red_obs,
                                 subject_weights = object$subject_weights,
                                 feature_weights = object$feature_weights,
                                 statistic = statistic)
  term_map_obs <- .dkge_between_term_map_from_coef(full_obs$coef,
                                                   design = object$design,
                                                   term = term,
                                                   drop = FALSE)
  feature_stat_obs <- .dkge_between_feature_stat(term_map_obs)
  feature_ge <- integer(length(feature_stat_obs))
  null_max <- numeric(B)

  null <- numeric(B)
  for (b in seq_len(B)) {
    idx <- .dkge_permute_within_blocks(blocks)
    Yb <- red_obs$fitted + red_obs$residuals[idx, , drop = FALSE]
    full_b <- .dkge_rrr_fit_core(X = X,
                                 Y = Yb,
                                 rank = object$rank,
                                 subject_weights = object$subject_weights,
                                 feature_weights = object$feature_weights,
                                 coef_ids = colnames(X),
                                 feature_ids = colnames(Y),
                                 subject_ids = rownames(Y),
                                 warn_rank = FALSE)
    red_b <- .dkge_rrr_fit_core(X = Xred,
                                Y = Yb,
                                rank = object$rank,
                                subject_weights = object$subject_weights,
                                feature_weights = object$feature_weights,
                                coef_ids = colnames(Xred),
                                feature_ids = colnames(Y),
                                subject_ids = rownames(Y),
                                warn_rank = FALSE)
    null[b] <- .dkge_between_stat(full_b, red_b,
                                  subject_weights = object$subject_weights,
                                  feature_weights = object$feature_weights,
                                  statistic = statistic)
    if (!identical(scope, "global")) {
      feature_stat_b <- .dkge_between_feature_stat(
        .dkge_between_term_map_from_coef(full_b$coef,
                                         design = object$design,
                                         term = term,
                                         drop = FALSE)
      )
      feature_ge <- feature_ge + as.integer(feature_stat_b >= feature_stat_obs)
      null_max[b] <- max(feature_stat_b)
    }
  }
  p <- (1 + sum(null >= stat_obs)) / (B + 1)
  feature_test <- if (identical(scope, "global")) {
    NULL
  } else {
    p_raw <- (1 + feature_ge) / (B + 1)
    p_adj <- switch(
      feature_adjust,
      none = p_raw,
      fdr = stats::p.adjust(p_raw, method = "fdr"),
      maxT = vapply(feature_stat_obs,
                    function(s) (1 + sum(null_max >= s)) / (B + 1),
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
  cat("  permutations :", x$B, "\n")
  cat("  scope        :", x$scope, "\n")
  print(x$summary, row.names = FALSE)
  if (!is.null(x$feature_tests)) {
    cat("  feature tests:", paste(names(x$feature_tests), collapse = ", "), "\n")
  }
  invisible(x)
}
