# dkge-transport.R
# Transport DKGE maps to a common reference parcellation via pluggable mappers.

.dkge_pairwise_sqdist <- function(A, B) {
  A <- as.matrix(A)
  B <- as.matrix(B)
  cpp_fun <- get0("pairwise_sqdist_cpp", mode = "function")
  if (!is.null(cpp_fun)) {
    return(cpp_fun(A, B))
  }

  n <- nrow(A); m <- nrow(B)
  An <- rowSums(A * A)
  Bn <- rowSums(B * B)
  D <- matrix(An, nrow = n, ncol = m)
  D <- D + matrix(Bn, nrow = n, ncol = m, byrow = TRUE)
  D <- D - 2 * tcrossprod(A, B)
  D[D < 0] <- 0
  D
}

.dkge_cost_matrix <- function(Aemb_s, Aemb_ref, X_s = NULL, X_ref = NULL,
                              lambda_emb = 1, lambda_spa = 0.5, sigma_mm = 15,
                              sizes_s = NULL, sizes_ref = NULL, lambda_size = 0) {
  C_emb <- .dkge_pairwise_sqdist(Aemb_s, Aemb_ref)
  C_spa <- 0
  if (!is.null(X_s) && !is.null(X_ref)) {
    C_spa <- .dkge_pairwise_sqdist(X_s / sigma_mm, X_ref / sigma_mm)
  }
  C <- lambda_emb * C_emb + lambda_spa * C_spa
  if (!is.null(sizes_s) && !is.null(sizes_ref) && lambda_size > 0) {
    la <- log(pmax(sizes_s, 1e-8))
    lb <- log(pmax(sizes_ref, 1e-8))
    C <- C + lambda_size * outer(la, lb, function(x, y) (x - y)^2)
  }
  C
}

.dkge_logsumexp <- function(x) {
  m <- max(x)
  if (!is.finite(m)) {
    return(m)
  }
  m + log(sum(exp(x - m)))
}

.dkge_sinkhorn_cache <- new.env(parent = emptyenv())
assign(".order", character(0), envir = .dkge_sinkhorn_cache)

.dkge_sinkhorn_cache_fetch <- function(key) {
  state <- .dkge_sinkhorn_cache[[key]]
  if (!is.null(state)) {
    order <- get(".order", envir = .dkge_sinkhorn_cache, inherits = FALSE)
    order <- c(setdiff(order, key), key)
    assign(".order", order, envir = .dkge_sinkhorn_cache)
  }
  state
}

.dkge_sinkhorn_cache_store <- function(key, state, max_entries = 64) {
  assign(key, state, envir = .dkge_sinkhorn_cache)
  order <- get(".order", envir = .dkge_sinkhorn_cache, inherits = FALSE)
  order <- c(setdiff(order, key), key)
  if (length(order) > max_entries) {
    drop <- order[seq_len(length(order) - max_entries)]
    if (length(drop)) {
      rm(list = drop, envir = .dkge_sinkhorn_cache)
    }
    order <- order[(length(order) - max_entries + 1):length(order)]
  }
  assign(".order", order, envir = .dkge_sinkhorn_cache)
}

.dkge_sinkhorn_plan <- function(C, mu, nu, epsilon = 0.05,
                                max_iter = 5000L, tol = 1e-4,
                                warm_start = TRUE,
                                return_diagnostics = FALSE) {
  if (!is.matrix(C) || !is.numeric(C) || any(!is.finite(C))) {
    stop("`C` must be a finite numeric cost matrix.", call. = FALSE)
  }
  if (length(mu) != nrow(C) || length(nu) != ncol(C)) {
    stop("`mu` and `nu` must match the rows and columns of `C`.", call. = FALSE)
  }
  mu <- as.numeric(mu)
  nu <- as.numeric(nu)
  if (any(!is.finite(mu)) || any(!is.finite(nu)) ||
      any(mu <= 0) || any(nu <= 0)) {
    stop("`mu` and `nu` must contain finite, strictly positive masses.", call. = FALSE)
  }
  if (length(epsilon) != 1L || !is.finite(epsilon) || epsilon <= 0 ||
      length(max_iter) != 1L || !is.finite(max_iter) || max_iter < 1 ||
      length(tol) != 1L || !is.finite(tol) || tol <= 0) {
    stop("`epsilon`, `max_iter`, and `tol` must be finite and positive.", call. = FALSE)
  }
  max_iter <- as.integer(max_iter)
  warm_start <- isTRUE(warm_start)
  if (abs(sum(mu) - sum(nu)) > 1e-6) {
    stop("mu and nu must sum to the same total mass")
  }

  sinkhorn_fun <- get0("sinkhorn_plan_cpp", mode = "function")
  if (is.null(sinkhorn_fun)) {
    stop("sinkhorn_plan_cpp() is unavailable; reinstall dkge with compiled code or install a binary build.",
         call. = FALSE)
  }

  # Duals are specific to every entry and its position, not merely to summary
  # moments of the cost and masses. A full serialized digest prevents different
  # problems with the same mean/SD from sharing a warm start.
  key <- digest::digest(list(C = C, mu = mu, nu = nu, epsilon = epsilon),
                        algo = "xxhash64", serialize = TRUE)
  state <- if (warm_start) .dkge_sinkhorn_cache_fetch(key) else NULL
  cache_hit <- !is.null(state)

  # An already-converged cached plan is the deterministic answer for an equal
  # or looser requested tolerance. Return it exactly instead of perturbing it
  # with another scaling cycle.
  if (cache_hit && !is.null(state$plan) &&
      is.finite(state$marginal_error) && state$marginal_error <= tol) {
    diagnostics <- list(
      converged = TRUE,
      iterations = 0L,
      marginal_error = state$marginal_error,
      row_marginal_error = state$row_marginal_error,
      column_marginal_error = state$column_marginal_error,
      cache_hit = TRUE,
      warm_started = TRUE,
      tolerance = tol
    )
    if (return_diagnostics) {
      return(list(plan = state$plan, diagnostics = diagnostics,
                  log_u = state$log_u, log_v = state$log_v))
    }
    return(state$plan)
  }
  log_u_init <- if (!is.null(state)) state$log_u else NULL
  log_v_init <- if (!is.null(state)) state$log_v else NULL

  res <- sinkhorn_fun(C, mu, nu, epsilon, as.integer(max_iter), tol,
                      log_u_init = log_u_init,
                      log_v_init = log_v_init,
                      keep_duals = TRUE)
  plan <- res$plan
  row_err <- max(abs(rowSums(plan) - mu))
  col_err <- max(abs(colSums(plan) - nu))
  marg_err <- max(row_err, col_err)
  converged <- is.finite(marg_err) && marg_err <= tol
  iterations <- as.integer(res$iterations %||% max_iter)
  diagnostics <- list(
    converged = converged,
    iterations = iterations,
    marginal_error = marg_err,
    row_marginal_error = row_err,
    column_marginal_error = col_err,
    cache_hit = cache_hit,
    warm_started = cache_hit,
    tolerance = tol
  )
  if (!converged) {
    warning(sprintf(
      "Sinkhorn did not converge in %d iterations (marginal error %.2e > tol %.2e); increase max_iter or epsilon.",
      iterations, marg_err, tol
    ), call. = FALSE)
  }
  # Failed iterates are deliberately not cached: a later solve must not inherit
  # a state that never satisfied the advertised numerical contract.
  if (warm_start && converged && !is.null(res$log_u) && !is.null(res$log_v)) {
    .dkge_sinkhorn_cache_store(key, list(log_u = res$log_u,
                                         log_v = res$log_v,
                                         plan = plan,
                                         marginal_error = marg_err,
                                         row_marginal_error = row_err,
                                         column_marginal_error = col_err,
                                         iterations = iterations))
  }
  if (return_diagnostics) {
    return(list(plan = plan, diagnostics = diagnostics,
                log_u = res$log_u, log_v = res$log_v))
  }
  plan
}

# Convert a joint coupling into the linear operator appropriate for the value
# being transported. Intensive fields are target-conditional expectations;
# extensive values distribute each source total across targets.
.dkge_transport_operator <- function(plan, mu, nu,
                                     value_type = c("intensive", "extensive")) {
  value_type <- match.arg(value_type)
  if (!is.matrix(plan) || length(mu) != nrow(plan) || length(nu) != ncol(plan)) {
    stop("`plan`, `mu`, and `nu` have incompatible dimensions.", call. = FALSE)
  }
  if (value_type == "intensive") {
    target_mass <- colSums(plan)
    if (any(!is.finite(target_mass)) || any(target_mass <= 0)) {
      stop("The transport plan has an empty or invalid target marginal.", call. = FALSE)
    }
    return(sweep(plan, 2L, target_mass, "/"))
  }
  source_mass <- rowSums(plan)
  if (any(!is.finite(source_mass)) || any(source_mass <= 0)) {
    stop("The transport plan has an empty or invalid source marginal.", call. = FALSE)
  }
  sweep(plan, 1L, source_mass, "/")
}

#' Clear cached dual variables for Sinkhorn warm-starts
#'
#' Releases memory held by the internal Sinkhorn cache. Call this after large
#' transport batches or before serialising objects.
#'
#' @return Logical `TRUE` invisibly.
#' @export
#' @examples
#' dkge_clear_sinkhorn_cache()
dkge_clear_sinkhorn_cache <- function() {
  keys <- setdiff(ls(.dkge_sinkhorn_cache, all.names = TRUE), ".order")
  if (length(keys)) {
    rm(list = keys, envir = .dkge_sinkhorn_cache)
  }
  assign(".order", character(0), envir = .dkge_sinkhorn_cache)
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# Helper for mapping lists --------------------------------------------------
# ---------------------------------------------------------------------------

.dkge_prepare_mapping_inputs <- function(loadings, centroids, sizes) {
  S <- length(loadings)
  stopifnot(length(centroids) == S)
  if (is.null(sizes)) {
    sizes <- lapply(loadings, function(A) rep(1, nrow(A)))
  } else {
    stopifnot(length(sizes) == S)
    sizes <- Map(function(sz, A) if (is.null(sz)) rep(1, nrow(A)) else sz,
                 sizes, loadings)
  }
  feat <- lapply(loadings, function(A) {
    A <- as.matrix(A)
    norms <- pmax(sqrt(rowSums(A^2)), 1e-8)
    A / norms
  })
  list(features = feat, sizes = sizes)
}

#' Resolve fold-safe functional features from contrast receipts
#'
#' Every subject's loading is expressed in its own held-out DKGE gauge. Before
#' row distances can be compared, map those gauges to the reference subject's
#' held-out basis with K-Procrustes. The scalar contrast values themselves are
#' not rotated and remain invariant to this coordinate choice.
#'
#' @keywords internal
#' @noRd
.dkge_fold_receipt_loadings <- function(fit, contrast_obj,
                                        reference_subject) {
  receipts <- contrast_obj$metadata$alignment_receipts %||% NULL
  if (!inherits(receipts, "dkge_alignment_receipts") || !length(receipts)) {
    .dkge_abort(
      paste0(
        "Inferential functional transport requires fold-safe alignment receipts. ",
        "Recompute contrasts with current LOSO/K-fold DKGE code, or select ",
        "`alignment_mode = 'descriptive'` explicitly."
      ),
      "dkge_alignment_receipt_error"
    )
  }
  S <- length(fit$Btil)
  if (length(receipts) != S) {
    .dkge_abort(
      "Alignment receipt count does not match the fitted subject cohort.",
      "dkge_alignment_receipt_error"
    )
  }
  reference_subject <- as.integer(reference_subject)
  if (length(reference_subject) != 1L || !is.finite(reference_subject) ||
      reference_subject < 1L || reference_subject > S) {
    .dkge_abort("`medoid` is not a valid reference-subject index.",
                "dkge_alignment_receipt_error")
  }
  expected_ids <- as.character(fit$subject_ids %||% seq_len(S))
  actual_ids <- unname(vapply(receipts, `[[`, character(1), "subject_id"))
  if (!identical(actual_ids, expected_ids)) {
    .dkge_abort(
      "Alignment receipt subject identity/order does not match the fit.",
      "dkge_alignment_receipt_error"
    )
  }
  ineligible <- !vapply(receipts, function(x) isTRUE(x$inference$eligible),
                        logical(1))
  if (any(ineligible)) {
    details <- paste0(actual_ids[ineligible], " (",
                      vapply(receipts[ineligible], function(x) {
                        x$inference$reason %||% "unknown"
                      }, character(1)), ")")
    .dkge_abort(
      paste0(
        "Fold receipts are not eligible for inferential functional transport: ",
        paste(details, collapse = ", "), "."
      ),
      "dkge_alignment_ineligible_error"
    )
  }

  reference_basis <- receipts[[reference_subject]]$basis
  loadings <- vector("list", S)
  gauges <- vector("list", S)
  names(loadings) <- names(gauges) <- expected_ids
  for (s in seq_len(S)) {
    receipt <- receipts[[s]]
    if (!identical(receipt$basis_hash, .dkge_object_hash(receipt$basis)) ||
        !identical(receipt$loadings_hash, .dkge_object_hash(receipt$loadings))) {
      .dkge_abort(
        sprintf("Alignment receipt for subject '%s' was mutated after fitting.",
                expected_ids[[s]]),
        "dkge_alignment_receipt_error"
      )
    }
    if (nrow(receipt$loadings) != receipt$preprocessing$n_clusters ||
        !identical(receipt$preprocessing$cluster_order_hash,
                   .dkge_object_hash(receipt$preprocessing$cluster_order))) {
      .dkge_abort(
        sprintf("Alignment receipt cluster provenance is invalid for subject '%s'.",
                expected_ids[[s]]),
        "dkge_alignment_receipt_error"
      )
    }
    pr <- dkge_procrustes_K(reference_basis, receipt$basis, fit$K,
                            allow_reflection = TRUE)
    aligned <- receipt$loadings %*% pr$R
    loadings[[s]] <- aligned
    delta <- pr$U_aligned - reference_basis
    gauges[[s]] <- list(
      reference_subject = reference_subject,
      rotation = pr$R,
      principal_cosines = pr$cosines,
      k_procrustes_residual = sqrt(max(
        0, sum(delta * (fit$K %*% delta))
      )) / sqrt(max(1, ncol(reference_basis))),
      source_basis_hash = receipt$basis_hash,
      aligned_loading_hash = .dkge_object_hash(aligned)
    )
  }
  list(
    loadings = loadings,
    receipts = receipts,
    gauges = gauges,
    source = "fold_receipts_k_procrustes",
    feature_source = "descriptive_adaptive",
    estimator_source = "same_data_rank_truncated",
    recompute_under_null = FALSE,
    inferentially_eligible = FALSE,
    eligibility_status = "ineligible",
    contract = paste0(
      attr(receipts, "inferential_contract"), " ",
      "The frozen functional-alignment court found static fold-loading ",
      "transport materially inflated; it is descriptive unless every ",
      "beta-dependent quantity is recomputed under the null."
    )
  )
}

.dkge_operator_diagnostics <- function(operator, plan, source_feat, target_feat,
                                       self_map = FALSE) {
  operator <- as.matrix(operator)
  absolute <- abs(operator)
  target_mass <- colSums(absolute)
  probabilities <- sweep(
    absolute, 2L, pmax(target_mass, .Machine$double.xmin), "/"
  )
  target_entropy <- apply(probabilities, 2L, function(p) {
    p <- p[p > 0]
    if (!length(p)) return(0)
    -sum(p * log(p))
  })
  effective_points <- exp(target_entropy)
  diagonal_mass <- rep(NA_real_, ncol(operator))
  self_mass <- NA_real_
  if (isTRUE(self_map) && nrow(operator) == ncol(operator)) {
    diagonal_mass <- diag(probabilities)
    if (!is.null(plan)) {
      plan_matrix <- as.matrix(plan)
      self_mass <- sum(diag(plan_matrix)) /
        max(sum(plan_matrix), .Machine$double.xmin)
    } else {
      self_mass <- mean(diagonal_mass)
    }
  }
  mapped_constant <- as.numeric(
    .dkge_apply_operator(operator, rep(1, nrow(operator)))
  )
  constant_error <- max(abs(mapped_constant - 1))
  constant_ratio <- sqrt(mean(mapped_constant^2))
  mapped_features <- .dkge_apply_operator(operator, source_feat)
  source_rms <- sqrt(mean(as.matrix(source_feat)^2))
  target_rms <- sqrt(mean(as.matrix(target_feat)^2))
  mapped_rms <- sqrt(mean(as.matrix(mapped_features)^2))
  feature_reconstruction <- sqrt(mean(
    (as.matrix(mapped_features) - as.matrix(target_feat))^2
  )) / max(target_rms, .Machine$double.xmin)
  list(
    point_spread = list(
      target_effective_points = as.numeric(effective_points),
      mean_effective_points = mean(effective_points),
      median_effective_points = stats::median(effective_points),
      max_effective_points = max(effective_points),
      target_diagonal_mass = diagonal_mass,
      mean_diagonal_mass = if (all(is.na(diagonal_mass))) {
        NA_real_
      } else {
        mean(diagonal_mass, na.rm = TRUE)
      },
      self_mass = self_mass
    ),
    amplitude_preservation = list(
      constant_rms_ratio = constant_ratio,
      constant_max_abs_error = constant_error,
      source_feature_rms = source_rms,
      mapped_feature_rms = mapped_rms,
      target_feature_rms = target_rms,
      mapped_to_target_rms_ratio = mapped_rms /
        max(target_rms, .Machine$double.xmin),
      feature_reconstruction_nrmse = feature_reconstruction,
      operator_l2_gain = norm(operator, "F") /
        sqrt(max(1, ncol(operator)))
    ),
    self_map = isTRUE(self_map)
  )
}

.dkge_epsilon_calibration <- function(mapper_spec) {
  calibration <- mapper_spec$params$epsilon_calibration %||% NULL
  if (is.null(calibration)) return(NULL)
  if (!identical(mapper_spec$strategy, "sinkhorn")) {
    .dkge_abort("Point-spread epsilon calibration is available only for Sinkhorn.",
                "dkge_alignment_calibration_error")
  }
  if (!is.list(calibration)) {
    .dkge_abort("`epsilon_calibration` must be a list.",
                "dkge_alignment_calibration_error")
  }
  target <- calibration$target_effective_points
  grid <- sort(unique(as.numeric(calibration$epsilon_grid)))
  data_hash <- calibration$calibration_data_hash
  if (!is.numeric(target) || length(target) != 1L || !is.finite(target) ||
      target < 1 || !length(grid) || any(!is.finite(grid)) || any(grid <= 0) ||
      !is.character(data_hash) || length(data_hash) != 1L ||
      is.na(data_hash) || !nzchar(data_hash)) {
    .dkge_abort(
      paste0(
        "Epsilon calibration requires target_effective_points >= 1, a ",
        "positive epsilon_grid, and a non-empty calibration_data_hash."
      ),
      "dkge_alignment_calibration_error"
    )
  }
  list(
    target_effective_points = as.numeric(target),
    epsilon_grid = grid,
    calibration_data_hash = data_hash,
    source = calibration$source %||% "separate_operator_calibration"
  )
}

.dkge_mapping_numerical_status <- function(mapping, mapper_spec = NULL) {
  strategy <- mapping$strategy %||% mapping$type %||%
    mapper_spec$strategy %||% NA_character_
  if (!identical(strategy, "sinkhorn")) {
    return(list(
      required = FALSE, passed = TRUE, converged = NA,
      marginal_error = NA_real_, tolerance = NA_real_, failures = character()
    ))
  }
  diagnostics <- mapping$fit_info$diagnostics %||%
    mapping$stats$diagnostics %||% list()
  tolerance <- diagnostics$tolerance %||% mapping$fit_info$tol %||%
    mapping$tol %||% mapper_spec$params$tol %||% 1e-4
  marginal_error <- diagnostics$marginal_error %||% NA_real_
  finite_operator <- !is.null(mapping$operator) &&
    all(is.finite(as.matrix(mapping$operator)))
  finite_plan <- !is.null(mapping$plan) &&
    all(is.finite(as.matrix(mapping$plan)))
  converged <- isTRUE(diagnostics$converged)
  marginal_ok <- is.numeric(marginal_error) && length(marginal_error) == 1L &&
    is.finite(marginal_error) && marginal_error <= tolerance
  failures <- character()
  if (!finite_operator) failures <- c(failures, "non-finite operator")
  if (!finite_plan) failures <- c(failures, "non-finite plan")
  if (!converged) failures <- c(failures, "solver did not converge")
  if (!marginal_ok) failures <- c(failures, "marginal error exceeds tolerance")
  list(
    required = TRUE,
    passed = finite_operator && finite_plan && converged && marginal_ok,
    converged = converged,
    marginal_error = as.numeric(marginal_error),
    tolerance = as.numeric(tolerance),
    finite_operator = finite_operator,
    finite_plan = finite_plan,
    failures = unique(failures)
  )
}

.dkge_assert_mapping_numerically_valid <- function(
    mapping,
    mapper_spec = NULL,
    context = "Sinkhorn mapping",
    class = "dkge_alignment_numerical_error") {
  status <- .dkge_mapping_numerical_status(mapping, mapper_spec)
  if (!isTRUE(status$passed)) {
    .dkge_abort(
      sprintf("%s is numerically invalid: %s.",
              context, paste(status$failures, collapse = "; ")),
      class
    )
  }
  invisible(status)
}

.dkge_fit_mapper_policy <- function(mapper_spec, source_feat, target_feat,
                                    source_weights, target_weights,
                                    source_xyz, target_xyz,
                                    self_map = FALSE) {
  calibration <- .dkge_epsilon_calibration(mapper_spec)
  fit_once <- function(spec) {
    fit_mapper(
      spec,
      source_feat = source_feat,
      target_feat = target_feat,
      source_weights = source_weights,
      target_weights = target_weights,
      source_xyz = source_xyz,
      target_xyz = target_xyz
    )
  }
  if (is.null(calibration)) return(fit_once(mapper_spec))

  candidates <- lapply(calibration$epsilon_grid, function(epsilon) {
    candidate_spec <- mapper_spec
    candidate_spec$params$epsilon <- epsilon
    candidate_spec$params$epsilon_calibration <- NULL
    mapping <- tryCatch(
      fit_once(candidate_spec),
      dkge_alignment_numerical_error = function(e) {
        .dkge_abort(
          sprintf("Epsilon-calibration candidate %.8g is numerically invalid: %s",
                  epsilon, conditionMessage(e)),
          "dkge_alignment_calibration_error"
        )
      }
    )
    .dkge_assert_mapping_numerically_valid(
      mapping, candidate_spec,
      context = sprintf("Epsilon-calibration candidate %.8g", epsilon),
      class = "dkge_alignment_calibration_error"
    )
    diagnostics <- .dkge_operator_diagnostics(
      mapping$operator, mapping$plan %||% NULL,
      source_feat, target_feat, self_map = self_map
    )
    list(mapping = mapping,
         effective_points = diagnostics$point_spread$mean_effective_points)
  })
  spread <- vapply(candidates, `[[`, numeric(1), "effective_points")
  score <- abs(log(pmax(spread, .Machine$double.xmin)) -
                 log(calibration$target_effective_points))
  best_score <- min(score)
  tied <- which(abs(score - best_score) <= 1e-12)
  best <- tied[[which.min(calibration$epsilon_grid[tied])]]
  mapping <- candidates[[best]]$mapping
  diagnostics <- mapping$fit_info$diagnostics %||% list()
  diagnostics$epsilon_calibration <- list(
    enabled = TRUE,
    target_effective_points = calibration$target_effective_points,
    epsilon_grid = calibration$epsilon_grid,
    achieved_effective_points = spread,
    objective = score,
    selected_epsilon = calibration$epsilon_grid[[best]],
    calibration_data_hash = calibration$calibration_data_hash,
    source = calibration$source
  )
  mapping$fit_info$epsilon <- calibration$epsilon_grid[[best]]
  mapping$fit_info$diagnostics <- diagnostics
  mapping
}

.dkge_run_mapper <- function(mapper_spec, value_list, feature_list, feature_ref,
                             centroids_list, centroid_ref,
                             sizes_list, size_ref,
                             reference_subject,
                             operators = NULL,
                             plans = NULL,
                             solver_diagnostics = NULL) {
  S <- length(feature_list)
  Q <- nrow(feature_ref)
  compute_values <- !is.null(value_list)
  if (compute_values) {
    stopifnot(length(value_list) == S)
  }

  mapped <- if (compute_values) matrix(NA_real_, S, Q) else NULL
  operators_out <- vector("list", S)
  plans_out <- vector("list", S)
  diagnostics_out <- vector("list", S)
  value_type <- mapper_spec$params$value_type %||% "intensive"

  for (s in seq_len(S)) {
    use_cached <- !is.null(operators) && length(operators) >= s &&
      !is.null(operators[[s]])

    if (use_cached) {
      operator <- operators[[s]]
      plan <- if (!is.null(plans) && length(plans) >= s) plans[[s]] else NULL
      diagnostics <- if (!is.null(solver_diagnostics) &&
                         length(solver_diagnostics) >= s) {
        solver_diagnostics[[s]]
      } else {
        list(reused_operator = TRUE)
      }
      diagnostics$reused_operator <- TRUE
    } else {
      map_fit <- .dkge_fit_mapper_policy(
        mapper_spec,
        source_feat = feature_list[[s]],
        target_feat = feature_ref,
        source_weights = sizes_list[[s]],
        target_weights = size_ref,
        source_xyz = centroids_list[[s]],
        target_xyz = centroid_ref,
        self_map = identical(s, as.integer(reference_subject))
      )
      operator <- map_fit$operator
      plan <- map_fit$plan %||% NULL
      diagnostics <- map_fit$fit_info$diagnostics %||% list()
      diagnostics$epsilon <- map_fit$fit_info$epsilon %||%
        mapper_spec$params$epsilon %||% NA_real_
    }

    operator_diagnostics <- .dkge_operator_diagnostics(
      operator, plan,
      feature_list[[s]], feature_ref,
      self_map = identical(s, as.integer(reference_subject))
    )
    diagnostics <- utils::modifyList(
      diagnostics %||% list(), operator_diagnostics, keep.null = TRUE
    )

    operators_out[[s]] <- operator
    plans_out[[s]] <- plan
    diagnostics_out[[s]] <- diagnostics

    if (compute_values) {
      vals <- value_list[[s]]
      mapped[s, ] <- as.numeric(.dkge_apply_operator(operator, vals))
    }
  }

  list(
    value = if (compute_values) apply(mapped, 2, stats::median, na.rm = TRUE) else NULL,
    subj_values = mapped,
    plans = plans_out,
    operators = operators_out,
    diagnostics = diagnostics_out,
    mapper = list(strategy = mapper_spec$strategy,
                  params = mapper_spec$params,
                  value_type = value_type),
    feature_ref = feature_ref,
    size_ref = size_ref
  )
}

.dkge_resolve_mapper_spec <- function(mapper, method, dots) {
  if (inherits(mapper, "dkge_mapper_spec")) {
    return(mapper)
  }

  if (is.character(mapper)) {
    if (identical(mapper, "sinkhorn_cpp")) {
      warning("`sinkhorn_cpp` is deprecated; `sinkhorn` already uses the compiled backend.",
              call. = FALSE)
      mapper <- "sinkhorn"
    }
    mapper <- dkge_mapper_spec(mapper)
  }

  if (inherits(mapper, "dkge_mapper_spec")) {
    if (length(dots)) {
      mapper$params[names(dots)] <- dots
    }
    return(mapper)
  }

  strategy <- method %||% "sinkhorn"
  if (identical(strategy, "sinkhorn_cpp")) {
    warning("`sinkhorn_cpp` is deprecated; `sinkhorn` already uses the compiled backend.",
            call. = FALSE)
    strategy <- "sinkhorn"
  }
  allowed <- c("epsilon", "max_iter", "tol", "lambda_emb", "lambda_spa",
               "sigma_mm", "lambda_size", "value_type", "warm_start")
  params <- dots[intersect(names(dots), allowed)]
  if (length(params)) {
    do.call(dkge_mapper_spec, c(list(type = strategy), params))
  } else {
    dkge_mapper_spec(strategy)
  }
}

#' Fingerprint every structural input to a fitted alignment
#'
#' @keywords internal
#' @noRd
.dkge_alignment_structural_receipt <- function(mapper_spec, feature_list,
                                               size_list, centroids,
                                               reference_subject,
                                               subject_ids = NULL,
                                               preprocessing = NULL,
                                               reference_support = NULL,
                                               functional_template = NULL,
                                               eligibility = NULL,
                                               reference_selection = NULL) {
  S <- length(feature_list)
  subject_ids <- as.character(subject_ids %||% names(feature_list) %||%
                                seq_len(S))
  components <- list(
    features = feature_list,
    centroids = centroids,
    masses = size_list,
    reference_subject = as.integer(reference_subject),
    mapper = mapper_spec,
    preprocessing = preprocessing,
    subject_order = subject_ids,
    value_semantics = mapper_spec$params$value_type %||% "intensive",
    reference_support = reference_support,
    functional_template = functional_template,
    eligibility = eligibility,
    reference_selection = reference_selection
  )
  component_hashes <- vapply(components, .dkge_object_hash, character(1))
  structure(
    list(
      schema_version = "1.0.0",
      component_hashes = component_hashes,
      structural_hash = .dkge_object_hash(component_hashes)
    ),
    class = c("dkge_alignment_structural_receipt", "list")
  )
}

#' Build immutable fitted alignment state from one mapper fit
#'
#' @keywords internal
#' @noRd
.dkge_new_fitted_alignment <- function(mapper_run, mapper_spec, feature_list,
                                       size_list, centroids,
                                       reference_subject,
                                       subject_ids = NULL,
                                       preprocessing = NULL,
                                       reference_selection = NULL) {
  S <- length(feature_list)
  subject_ids <- as.character(subject_ids %||% names(feature_list) %||%
                                seq_len(S))
  contract <- .dkge_alignment_objects_from_inputs(
    feature_list, size_list, centroids, reference_subject,
    preprocessing = preprocessing,
    reference_selection = reference_selection
  )
  contract$eligibility <- .dkge_alignment_solver_eligibility(
    contract$eligibility,
    mapper_spec = mapper_spec,
    diagnostics = mapper_run$diagnostics,
    operators = mapper_run$operators,
    plans = mapper_run$plans,
    subject_ids = subject_ids
  )
  receipt <- .dkge_alignment_structural_receipt(
    mapper_spec, feature_list, size_list, centroids, reference_subject,
    subject_ids = subject_ids, preprocessing = preprocessing,
    reference_support = contract$support,
    functional_template = contract$template,
    eligibility = contract$eligibility,
    reference_selection = reference_selection
  )
  solution_hashes <- vapply(
    list(
      operators = mapper_run$operators,
      plans = mapper_run$plans,
      diagnostics = mapper_run$diagnostics
    ),
    .dkge_object_hash,
    character(1)
  )
  out <- list(
      schema_version = "1.0.0",
      operators = mapper_run$operators,
      plans = mapper_run$plans,
      diagnostics = mapper_run$diagnostics,
      mapper_spec = mapper_spec,
      feature_list = feature_list,
      size_list = size_list,
      feature_ref = feature_list[[reference_subject]],
      size_ref = size_list[[reference_subject]],
      centroids = centroids,
      medoid = as.integer(reference_subject),
      reference_subject = as.integer(reference_subject),
      subject_ids = subject_ids,
      preprocessing = preprocessing,
      reference_support = contract$support,
      functional_template = contract$template,
      feature_source = contract$feature_source,
      estimator_source = contract$estimator_source,
      eligibility = contract$eligibility,
      reference_selection = reference_selection,
      reference_method = reference_selection$method %||% "explicit_legacy",
      reference_is_medoid = !is.null(reference_selection$medoid),
      structural_receipt = receipt,
      solution_hashes = solution_hashes
  )
  out$fitted_hash <- .dkge_object_hash(c(
    receipt$structural_hash,
    solution_hashes,
    eligibility = .dkge_object_hash(contract$eligibility)
  ))
  structure(out, class = c("dkge_fitted_alignment", "list"))
}

#' Validate a fitted alignment against live structural arguments
#'
#' @keywords internal
#' @noRd
.dkge_validate_fitted_alignment <- function(alignment, mapper_spec,
                                            feature_list, size_list, centroids,
                                            reference_subject,
                                            subject_ids = NULL,
                                            preprocessing = NULL,
                                            reference_selection = NULL) {
  .dkge_validate_fitted_alignment_object(alignment)
  required <- c("operators", "mapper_spec", "feature_list", "size_list",
                "feature_ref", "size_ref", "centroids", "medoid")
  absent <- setdiff(required, names(alignment))
  if (length(absent)) {
    .dkge_abort(
      sprintf("Fitted alignment is missing structural field(s): %s.",
              paste(absent, collapse = ", ")),
      "dkge_alignment_cache_mismatch"
    )
  }
  cached_reference <- alignment$reference_subject %||% alignment$medoid
  cached_subject_ids <- alignment$subject_ids %||%
    as.character(names(alignment$feature_list) %||%
                   seq_along(alignment$feature_list))
  cached_preprocessing <- alignment$preprocessing %||% NULL
  cached_now <- .dkge_alignment_structural_receipt(
    alignment$mapper_spec,
    alignment$feature_list,
    alignment$size_list,
    alignment$centroids,
    cached_reference,
    subject_ids = cached_subject_ids,
    preprocessing = cached_preprocessing,
    reference_support = alignment$reference_support,
    functional_template = alignment$functional_template,
    eligibility = alignment$eligibility,
    reference_selection = alignment$reference_selection
  )
  stored <- alignment$structural_receipt %||% cached_now
  if (!identical(stored$structural_hash, cached_now$structural_hash)) {
    bad <- names(stored$component_hashes)[
      stored$component_hashes != cached_now$component_hashes
    ]
    .dkge_abort(
      sprintf("Fitted alignment state was mutated after fitting: %s.",
              paste(bad, collapse = ", ")),
      "dkge_alignment_cache_mismatch"
    )
  }

  live <- .dkge_alignment_structural_receipt(
    mapper_spec, feature_list, size_list, centroids, reference_subject,
    subject_ids = subject_ids, preprocessing = preprocessing,
    reference_support = alignment$reference_support,
    functional_template = alignment$functional_template,
    eligibility = alignment$eligibility,
    reference_selection = reference_selection
  )
  mismatch <- names(stored$component_hashes)[
    stored$component_hashes != live$component_hashes
  ]
  if (length(mismatch)) {
    .dkge_abort(
      sprintf(
        "Fitted alignment cannot be reused because structural input(s) changed: %s.",
        paste(mismatch, collapse = ", ")
      ),
      "dkge_alignment_cache_mismatch"
    )
  }
  list(
    alignment = alignment,
    provenance = list(
      hit = TRUE,
      validation = "exact_structural_match",
      structural_hash = stored$structural_hash,
      legacy_upgraded = FALSE
    )
  )
}

.dkge_transport_to_medoid <- function(mapper_spec, value_list, loadings,
                                      centroids, sizes, medoid,
                                      transport_cache = NULL,
                                      subject_ids = NULL,
                                      preprocessing = NULL,
                                      reference_selection = NULL,
                                      alignment_features = NULL) {
  subject_ids <- as.character(subject_ids %||% names(loadings) %||%
                                paste0("subject", seq_along(loadings)))
  canonical <- .dkge_reference_centroids(centroids, subject_ids)
  centroids <- canonical$centroids
  subject_ids <- canonical$subject_ids
  sizes <- .dkge_reference_sizes(sizes, centroids, subject_ids)
  if (is.null(reference_selection) && !is.null(transport_cache)) {
    reference_selection <- transport_cache$reference_selection %||% NULL
  }
  if (is.null(reference_selection)) {
    reference_selection <- dkge_select_reference_subject(
      centroids = centroids,
      sizes = sizes,
      method = "explicit",
      reference_subject = medoid,
      subject_ids = subject_ids,
      provenance = list(source = "explicit_transport_argument")
    )
  } else {
    .dkge_validate_reference_selection(
      reference_selection,
      subject_ids = subject_ids,
      centroids = centroids,
      sizes = sizes,
      alignment_features = alignment_features,
      mapper_spec = mapper_spec
    )
    if (!is.null(reference_selection$alignment_features_hash) &&
        is.null(alignment_features)) {
      .dkge_abort(
        "Functional reference selection requires its typed alignment features.",
        "dkge_reference_selection_error"
      )
    }
    selected <- reference_selection$reference_subject
    if (!identical(as.integer(medoid), as.integer(selected))) {
      .dkge_abort(
        "Explicit reference index conflicts with the fitted reference selection.",
        "dkge_reference_selection_error"
      )
    }
  }
  medoid <- reference_selection$reference_subject
  prep <- .dkge_prepare_mapping_inputs(loadings, centroids, sizes)
  cache_provenance <- NULL
  if (is.null(transport_cache)) {
    feature_list <- prep$features
    size_list <- prep$sizes
    operators <- NULL
    feature_ref <- feature_list[[medoid]]
    size_ref <- size_list[[medoid]]
    centroids_used <- centroids
    medoid_used <- medoid
  } else {
    validated <- .dkge_validate_fitted_alignment(
      transport_cache, mapper_spec, prep$features, prep$sizes, centroids,
      medoid, subject_ids = subject_ids, preprocessing = preprocessing,
      reference_selection = reference_selection
    )
    alignment <- validated$alignment
    cache_provenance <- validated$provenance
    feature_list <- alignment$feature_list
    size_list <- alignment$size_list
    operators <- alignment$operators
    plans <- alignment$plans
    solver_diagnostics <- alignment$diagnostics
    feature_ref <- alignment$feature_ref
    size_ref <- alignment$size_ref
    centroids_used <- alignment$centroids
    medoid_used <- alignment$reference_subject %||% alignment$medoid
  }
  if (is.null(transport_cache)) {
    plans <- NULL
    solver_diagnostics <- NULL
  }

  centroid_ref <- centroids_used[[medoid_used]]
  res <- .dkge_run_mapper(mapper_spec, value_list, feature_list, feature_ref,
                          centroids_used, centroid_ref, size_list, size_ref,
                          medoid_used,
                          operators = operators,
                          plans = plans,
                          solver_diagnostics = solver_diagnostics)

  res$feature_list <- feature_list
  res$size_list <- size_list
  res$centroids <- centroids_used
  res$medoid <- medoid_used
  if (is.null(transport_cache)) {
    fitted <- .dkge_new_fitted_alignment(
      res, mapper_spec, feature_list, size_list, centroids_used, medoid_used,
      subject_ids = subject_ids, preprocessing = preprocessing,
      reference_selection = reference_selection
    )
    cache_provenance <- list(
      hit = FALSE,
      validation = "new_fit",
      structural_hash = fitted$structural_receipt$structural_hash,
      legacy_upgraded = FALSE
    )
  } else {
    fitted <- transport_cache
  }
  res$fitted_alignment <- fitted
  res$reference_selection <- reference_selection
  res$cache_provenance <- cache_provenance
  res
}

#' Prepare subject-to-medoid transport operators for reuse
#'
#' Computes and caches subject-to-medoid transport matrices so downstream
#' routines (e.g. bootstraps) can reuse a fixed consensus mapping without
#' re-solving the transport problem on every call.
#'
#' @param fit A `dkge` object.
#' @param centroids List of subject centroid matrices. Defaults to the centroids
#'   stored on `fit` or `fit$input`.
#' @param loadings Optional list of subject loadings (`P_s x r`). When omitted,
#'   they are recomputed from `fit$Btil` or the supplied `betas`.
#' @param betas Optional list of subject betas used to recompute loadings when
#'   `loadings` is `NULL`.
#' @param sizes Optional list of cluster masses.
#' @param mapper Mapper specification or shorthand passed to
#'   [dkge_mapper_spec()].
#' @param medoid Index (1-based) of the reference subject.
#' @param reference_selection Optional typed selection from
#'   [dkge_select_reference_subject()]. When supplied, its selected subject is
#'   authoritative and `medoid` is only a compatibility alias.
#' @param alignment_features Optional typed feature object used both for
#'   reference selection and mapper fitting.
#' @param preprocessing Optional immutable provenance for feature construction.
#'   Caller-authored provenance is accepted only when it exactly matches a typed
#'   `alignment_features` object. Loose loadings/betas are always recorded as
#'   descriptive and inferentially ineligible.
#' @param ... Additional mapper arguments such as `epsilon` or `lambda_spa`.
#'
#' @return A list containing cached application `operators`, joint transport
#'   `plans`, solver `diagnostics`, `mapper_spec`, `feature_list`, `size_list`,
#'   `feature_ref`, `size_ref`, `centroids`, and `medoid`.
#' @keywords internal
#' @export
dkge_prepare_transport <- function(fit,
                                   centroids = NULL,
                                   loadings = NULL,
                                   betas = NULL,
                                   sizes = NULL,
                                   mapper = "sinkhorn",
                                   medoid = 1L,
                                   reference_selection = NULL,
                                   alignment_features = NULL,
                                   preprocessing = NULL,
                                   ...) {
  stopifnot(inherits(fit, "dkge"))

  loose_source <- if (!is.null(loadings)) {
    "caller_supplied_loose_loadings"
  } else if (!is.null(betas)) {
    "caller_supplied_loose_betas"
  } else {
    "full_fit_loadings"
  }

  medoid_missing <- missing(medoid)
  if (!is.null(reference_selection)) {
    .dkge_validate_reference_selection(reference_selection)
    selected <- reference_selection$reference_subject
    if (!medoid_missing && !identical(as.integer(medoid), as.integer(selected))) {
      .dkge_abort("`medoid` conflicts with `reference_selection`.",
                  "dkge_reference_selection_error")
    }
    medoid <- selected
  }

  if (!is.null(alignment_features)) {
    .dkge_validate_alignment_features(alignment_features, fit = fit)
    if (!is.null(loadings) || !is.null(betas)) {
      .dkge_abort(
        "Typed `alignment_features` cannot be combined with loose loadings or betas.",
        "dkge_alignment_feature_error"
      )
    }
    loadings <- alignment_features$features
    feature_preprocessing <- .dkge_alignment_feature_preprocessing(
      alignment_features
    )
    if (!is.null(preprocessing) &&
        !identical(preprocessing, feature_preprocessing)) {
      .dkge_abort("Feature preprocessing provenance was supplied twice and differs.",
                  "dkge_alignment_feature_error")
    }
    preprocessing <- feature_preprocessing
  } else {
    if (!is.null(preprocessing)) {
      .dkge_abort(
        paste0(
          "Caller-authored `preprocessing` cannot establish inferential ",
          "provenance for loose loadings or betas; supply a validated typed ",
          "`alignment_features` object instead."
        ),
        "dkge_alignment_preprocessing_error"
      )
    }
    preprocessing <- list(
      schema_version = "1.0.0",
      source = loose_source,
      feature_source = "descriptive_adaptive",
      estimator_source = "descriptive",
      recompute_under_null = FALSE,
      inferentially_eligible = FALSE,
      eligibility_status = "ineligible",
      feature_provenance = list(
        source = loose_source,
        verified = FALSE,
        contract = paste0(
          "Loose loadings/betas are accepted only for descriptive transport; ",
          "strings and caller hashes are not independent-data evidence."
        )
      )
    )
  }

  if (is.null(loadings)) {
    if (!is.null(betas)) {
      loadings <- dkge_predict_loadings(fit, betas)
    } else if (!is.null(fit$Btil)) {
      loadings <- .dkge_fit_subject_loadings(fit)
    } else {
      stop("Provide betas or pre-computed loadings.")
    }
  }

  if (is.null(centroids)) {
    centroids <- fit$centroids %||% fit$input$centroids %||%
      stop("Centroids required for transport preparation; none found in fit or arguments.")
  }

  stopifnot(length(centroids) == length(loadings))

  subject_ids <- as.character(
    fit$subject_ids %||% names(loadings) %||%
      paste0("subject", seq_along(loadings))
  )
  canonical <- .dkge_reference_centroids(centroids, subject_ids)
  centroids <- canonical$centroids
  subject_ids <- canonical$subject_ids
  sizes <- .dkge_reference_sizes(sizes, centroids, subject_ids)

  mapper_spec <- .dkge_resolve_mapper_spec(mapper, method = NULL, dots = list(...))

  prep <- .dkge_prepare_mapping_inputs(loadings, centroids, sizes)
  feature_list <- prep$features
  size_list <- prep$sizes
  reference_selection <- reference_selection %||%
    dkge_select_reference_subject(
      centroids = centroids,
      sizes = size_list,
      method = "explicit",
      reference_subject = medoid,
      subject_ids = subject_ids,
      provenance = list(source = "explicit_prepare_argument")
    )
  .dkge_validate_reference_selection(
    reference_selection,
    subject_ids = subject_ids,
    centroids = centroids,
    sizes = size_list,
    alignment_features = alignment_features,
    mapper_spec = mapper_spec
  )
  if (!is.null(reference_selection$alignment_features_hash) &&
      is.null(alignment_features)) {
    .dkge_abort(
      "Functional reference selection requires its typed alignment features.",
      "dkge_reference_selection_error"
    )
  }
  feature_ref <- feature_list[[medoid]]
  size_ref <- size_list[[medoid]]
  centroid_ref <- centroids[[medoid]]

  mapper_run <- .dkge_run_mapper(mapper_spec,
                                 value_list = NULL,
                                 feature_list = feature_list,
                                 feature_ref = feature_ref,
                                 centroids_list = centroids,
                                 centroid_ref = centroid_ref,
                                 sizes_list = size_list,
                                 size_ref = size_ref,
                                 reference_subject = medoid,
                                 operators = NULL)

  .dkge_new_fitted_alignment(
    mapper_run = mapper_run,
    mapper_spec = mapper_spec,
    feature_list = feature_list,
    size_list = size_list,
    centroids = centroids,
    reference_subject = medoid,
    subject_ids = subject_ids,
    preprocessing = preprocessing,
    reference_selection = reference_selection
  )
}

# ---------------------------------------------------------------------------
# Public wrappers -----------------------------------------------------------
# ---------------------------------------------------------------------------

#' Low-level descriptive Sinkhorn transport to a reference subject
#'
#' This compatibility primitive fits correspondence from caller-supplied loose
#' features (`A_list`) and therefore returns an inferentially ineligible fitted
#' alignment. It is useful for descriptive transport diagnostics. Inferential
#' workflows should use typed features with [dkge_prepare_alignment()] or
#' [dkge_transport_contrasts_to_reference()].
#'
#' @param v_list List of subject-level cluster values (length P_s each).
#' @param A_list List of subject loadings (P_s x r).
#' @param centroids List of subject cluster centroids (each P_s x 3 matrix).
#' @param sizes Optional list of cluster masses (defaults to uniform weights).
#' @param reference_subject Integer index of the fixed reference subject
#'   (1-based). This argument does not select or certify a medoid.
#' @param lambda_emb,lambda_spa Cost weights for embedding and spatial terms.
#' @param sigma_mm Spatial rescaling (millimetres).
#' @param epsilon,max_iter,tol Sinkhorn parameters.
#' @param value_type Value semantics. `"intensive"` transports field values as
#'   target-conditional averages and preserves constants; `"extensive"`
#'   distributes source totals and preserves their sum.
#' @param warm_start Logical; reuse converged dual variables for an identical
#'   cost-and-mass problem.
#' @param transport_cache Optional fitted alignment returned by
#'   [dkge_prepare_transport()]. Cached operators are reused only after every
#'   structural input is fingerprint-validated.
#' @return List containing summary statistics, transported subject maps, and
#'   per-subject joint `plans`, application `operators`, and solver
#'   `diagnostics`.
#' @export
dkge_transport_to_reference_sinkhorn <- function(v_list, A_list, centroids,
                                                 sizes = NULL,
                                                 reference_subject,
                                                 lambda_emb = 1,
                                                 lambda_spa = 0.5,
                                                 sigma_mm = 15,
                                                 epsilon = 0.05,
                                                 max_iter = 5000L,
                                                 tol = 1e-4,
                                                 value_type = c("intensive", "extensive"),
                                                 warm_start = TRUE,
                                                 transport_cache = NULL) {
  sizes_supplied <- !missing(sizes)
  mapper_supplied <- !missing(lambda_emb) || !missing(lambda_spa) ||
    !missing(sigma_mm) || !missing(epsilon) || !missing(max_iter) ||
    !missing(tol) || !missing(value_type) || !missing(warm_start)
  value_type <- match.arg(value_type)
  if (!is.null(transport_cache) && !sizes_supplied &&
      inherits(transport_cache, "dkge_fitted_alignment")) {
    sizes <- transport_cache$size_list
  }
  if (!is.null(transport_cache) && !mapper_supplied) {
    mapper_spec <- transport_cache$mapper_spec
  } else {
    mapper_spec <- dkge_mapper_spec("sinkhorn",
                                    lambda_emb = lambda_emb,
                                    lambda_spa = lambda_spa,
                                    sigma_mm = sigma_mm,
                                    epsilon = epsilon,
                                    max_iter = max_iter,
                                    tol = tol,
                                    value_type = value_type,
                                    warm_start = warm_start)
  }
  result <- .dkge_transport_to_medoid(
    mapper_spec, v_list, A_list, centroids, sizes, reference_subject,
    transport_cache = transport_cache
  )
  if (!isTRUE(result$fitted_alignment$eligibility$solver_converged)) {
    .dkge_abort(
      paste0(
        "Sinkhorn transport did not satisfy its numerical contract: ",
        result$fitted_alignment$eligibility$reason
      ),
      "dkge_alignment_numerical_error"
    )
  }
  result
}

#' Transport cluster values to a medoid (deprecated)
#'
#' Compatibility alias for [dkge_transport_to_reference_sinkhorn()].
#'
#' @inheritParams dkge_transport_to_reference_sinkhorn
#' @param medoid Deprecated name for `reference_subject`.
#' @export
dkge_transport_to_medoid_sinkhorn <- function(v_list, A_list, centroids,
                                              sizes = NULL, medoid,
                                              lambda_emb = 1,
                                              lambda_spa = 0.5,
                                              sigma_mm = 15,
                                              epsilon = 0.05,
                                              max_iter = 5000L,
                                              tol = 1e-4,
                                              value_type = c("intensive", "extensive"),
                                              warm_start = TRUE,
                                              transport_cache = NULL) {
  sizes_supplied <- !missing(sizes)
  mapper_supplied <- !missing(lambda_emb) || !missing(lambda_spa) ||
    !missing(sigma_mm) || !missing(epsilon) || !missing(max_iter) ||
    !missing(tol) || !missing(value_type) || !missing(warm_start)
  .Deprecated(
    "dkge_transport_to_reference_sinkhorn",
    package = "dkge",
    msg = paste0(
      "`dkge_transport_to_medoid_sinkhorn()` is deprecated; use the ",
      "explicit fixed-reference name. Both are low-level descriptive APIs."
    )
  )
  if (!is.null(transport_cache) && !mapper_supplied) {
    if (sizes_supplied) {
      return(dkge_transport_to_reference_sinkhorn(
        v_list, A_list, centroids, sizes = sizes,
        reference_subject = medoid,
        transport_cache = transport_cache
      ))
    }
    return(dkge_transport_to_reference_sinkhorn(
      v_list, A_list, centroids,
      reference_subject = medoid,
      transport_cache = transport_cache
    ))
  }
  dkge_transport_to_reference_sinkhorn(
    v_list, A_list, centroids, sizes = sizes,
    reference_subject = medoid,
    lambda_emb = lambda_emb, lambda_spa = lambda_spa, sigma_mm = sigma_mm,
    epsilon = epsilon, max_iter = max_iter, tol = tol,
    value_type = value_type, warm_start = warm_start,
    transport_cache = transport_cache
  )
}

#' @rdname dkge_transport_to_reference_sinkhorn
#' @param medoid Deprecated name for `reference_subject`.
#' @param return_plans Logical; if TRUE, include transport plans in the output.
#' @description
#' `dkge_transport_to_medoid_sinkhorn_cpp()` is a deprecated compatibility
#' alias. The reference-oriented function already uses the compiled Sinkhorn
#' backend.
#' @export
dkge_transport_to_medoid_sinkhorn_cpp <- function(v_list, A_list, centroids, sizes = NULL,
                                                  medoid,
                                                  lambda_emb = 1, lambda_spa = 0.5, sigma_mm = 15,
                                                  epsilon = 0.05, max_iter = 5000L, tol = 1e-4,
                                                  value_type = c("intensive", "extensive"),
                                                  warm_start = TRUE,
                                                  return_plans = FALSE,
                                                  transport_cache = NULL) {
  sizes_supplied <- !missing(sizes)
  mapper_supplied <- !missing(lambda_emb) || !missing(lambda_spa) ||
    !missing(sigma_mm) || !missing(epsilon) || !missing(max_iter) ||
    !missing(tol) || !missing(value_type) || !missing(warm_start)
  .Deprecated("dkge_transport_to_reference_sinkhorn")
  if (!is.null(transport_cache) && !mapper_supplied) {
    res <- if (sizes_supplied) {
      dkge_transport_to_reference_sinkhorn(
        v_list, A_list, centroids, sizes = sizes,
        reference_subject = medoid,
        transport_cache = transport_cache
      )
    } else {
      dkge_transport_to_reference_sinkhorn(
        v_list, A_list, centroids,
        reference_subject = medoid,
        transport_cache = transport_cache
      )
    }
  } else {
    value_type <- match.arg(value_type)
    res <- dkge_transport_to_reference_sinkhorn(
      v_list, A_list, centroids, sizes = sizes,
      reference_subject = medoid,
      lambda_emb = lambda_emb, lambda_spa = lambda_spa,
      sigma_mm = sigma_mm, epsilon = epsilon, max_iter = max_iter, tol = tol,
      value_type = value_type, warm_start = warm_start,
      transport_cache = transport_cache
    )
  }
  if (!return_plans) {
    res$plans <- NULL
  }
  res
}

#' Transport component loadings for legacy descriptive display
#'
#' This compatibility helper learns correspondence from full-fit component
#' loadings and is therefore descriptive/ineligible. It cannot supply an
#' inferential alignment. New workflows should build typed alignment features
#' and call [dkge_prepare_alignment()].
#'
#' @param fit A `dkge` object used to compute the loadings.
#' @param medoid Integer index of the reference subject (1-based).
#' @param centroids List of subject cluster centroids (each P_s x 3 matrix).
#' @param loadings Optional list of subject loadings (P_s x r). When omitted,
#'   they are recomputed from `betas`.
#' @param betas Optional list of subject betas used to recompute loadings when
#'   `loadings` is `NULL`.
#' @param sizes Optional list of cluster masses (defaults to uniform weights).
#' @param mapper Optional mapper specification created by [dkge_mapper_spec()].
#'   When `NULL`, defaults to Sinkhorn with the supplied parameters.
#' @param method Mapper strategy (`"sinkhorn"`, `"ridge"`, or `"ols"`). The
#'   legacy `"sinkhorn_cpp"` name is a deprecated alias for `"sinkhorn"`.
#' @param transport_cache Optional fitted alignment from
#'   [dkge_prepare_transport()]. Reuse requires exact structural fingerprints.
#' @param ... Additional parameters passed when building the default mapper
#'   specification (e.g. `epsilon`, `lambda_emb`).
#' @return List with `group` (medoid cluster vectors per component),
#'   `subjects` (per-subject transported values), and `cache` (transport cache
#'   reused for future calls).
#' @export
dkge_transport_loadings_to_medoid <- function(fit, medoid, centroids,
                                               loadings = NULL,
                                               betas = NULL,
                                               sizes = NULL,
                                               mapper = NULL,
                                               method = c("sinkhorn", "ridge", "ols", "sinkhorn_cpp"),
                                               transport_cache = NULL,
                                               ...) {
  .Deprecated(
    "dkge_prepare_alignment",
    package = "dkge",
    msg = paste0(
      "`dkge_transport_loadings_to_medoid()` is deprecated and descriptive ",
      "only; use typed alignment features with `dkge_prepare_alignment()`."
    )
  )
  stopifnot(inherits(fit, "dkge"))
  if (is.null(loadings)) {
    if (!is.null(betas)) {
      loadings <- dkge_predict_loadings(fit, betas)
    } else if (!is.null(fit$Btil)) {
      loadings <- .dkge_fit_subject_loadings(fit)
    } else {
      stop("Provide betas or pre-computed loadings.")
    }
  }
  S <- length(loadings)
  stopifnot(length(centroids) == S)

  dots <- list(...)
  mapper_supplied <- !is.null(mapper) || !missing(method) || length(dots) > 0L
  if (!is.null(transport_cache) && !mapper_supplied) {
    mapper_spec <- transport_cache$mapper_spec
  } else {
    mapper_spec <- .dkge_resolve_mapper_spec(mapper, method = method[1], dots = dots)
  }

  if (is.null(sizes)) {
    sizes <- lapply(loadings, function(A) rep(1, nrow(A)))
  } else {
    sizes <- Map(function(sz, A) if (is.null(sz)) rep(1, nrow(A)) else sz, sizes, loadings)
  }
  rank <- ncol(loadings[[1]])

  subj_vals <- vector("list", rank)
  group_vals <- vector("list", rank)
  cache_local <- transport_cache
  for (j in seq_len(rank)) {
    v_list <- lapply(loadings, function(A) A[, j])
    tr <- .dkge_transport_to_medoid(mapper_spec, v_list, loadings, centroids,
                                    sizes = sizes, medoid = medoid,
                                    transport_cache = cache_local,
                                    subject_ids = fit$subject_ids)
    subj_vals[[j]] <- tr$subj_values
    group_vals[[j]] <- tr$value
    if (is.null(cache_local)) {
      cache_local <- tr$fitted_alignment
    }
  }
  list(
    group = group_vals,
    subjects = subj_vals,
    cache = cache_local,
    eligibility = cache_local$eligibility,
    metadata = list(status = "descriptive", inferential = FALSE)
  )
}


#' Transport subject contrasts to a medoid parcellation
#'
#' @param fit A `dkge` object used to compute the contrasts.
#' @param contrast_obj A `dkge_contrasts` result.
#' @param medoid Integer index of the reference subject (1-based).
#' @param centroids List of subject cluster centroids (each P_s x 3 matrix).
#' @param loadings Optional list of subject loadings (P_s x r).
#' @param betas Optional list of subject betas used to recompute loadings.
#' @param sizes Optional list of cluster masses.
#' @param mapper Optional mapper specification created by [dkge_mapper_spec()].
#' @param method Mapper strategy (`"sinkhorn"`, `"ridge"`, or `"ols"`). The
#'   legacy `"sinkhorn_cpp"` name is a deprecated alias for `"sinkhorn"`.
#' @param transport_cache Optional fitted alignment from
#'   [dkge_prepare_transport()]. Reuse requires exact structural fingerprints.
#' @param reference_selection Optional typed selection from
#'   [dkge_select_reference_subject()]. When supplied, its selected reference
#'   is authoritative and `medoid` is only a compatibility alias.
#' @param alignment_features Optional typed object from
#'   [dkge_alignment_features()]. Independent and contrast-orthogonal modes
#'   require this object so feature construction and provenance cannot drift
#'   apart.
#' @param alignment_mode Functional-feature provenance. The default
#'   `"fold_safe"` requires typed held-out receipts stored on `contrast_obj` and
#'   K-Procrustes-aligns their low-r gauges; the frozen court classifies this
#'   legacy path as descriptive/ineligible. `"independent"` consumes
#'   provenance-declared over-ranked features from separately identified data.
#'   `"contrast_orthogonal"` consumes
#'   covariance-residualized kernel-image features and is explicitly
#'   approximate. `"descriptive"` permits transductive/full-fit loadings.
#' @param ... Additional parameters passed when building the default mapper
#'   specification.
#' @return Named list of transport results (one per contrast) with an attached
#'   `cache` element for reuse.
#' @keywords internal
#' @noRd
.dkge_transport_contrasts_to_reference_core <- function(fit, contrast_obj,
                                               medoid = NULL,
                                               centroids = NULL,
                                               loadings = NULL,
                                               betas = NULL,
                                               sizes = NULL,
                                               mapper = NULL,
                                               method = c("sinkhorn", "ridge", "ols", "sinkhorn_cpp"),
                                               transport_cache = NULL,
                                               reference_selection = NULL,
                                               alignment_features = NULL,
                                               alignment_mode = c(
                                                 "fold_safe", "independent",
                                                 "contrast_orthogonal", "descriptive"
                                               ),
                                               .method_missing = missing(method),
                                               .alignment_mode_missing =
                                                 missing(alignment_mode),
                                               ...) {
  stopifnot(inherits(fit, "dkge"), inherits(contrast_obj, "dkge_contrasts"))
  medoid_missing <- missing(medoid) || is.null(medoid)
  if (!is.null(reference_selection)) {
    .dkge_validate_reference_selection(reference_selection)
    selected <- reference_selection$reference_subject
    if (!medoid_missing && !identical(as.integer(medoid),
                                      as.integer(selected))) {
      .dkge_abort("`medoid` conflicts with `reference_selection`.",
                  "dkge_reference_selection_error")
    }
    medoid <- selected
  } else if (medoid_missing) {
    .dkge_abort(
      paste0(
        "Supply a typed `reference_selection`, or pass the legacy `medoid` ",
        "index explicitly."
      ),
      "dkge_reference_selection_error"
    )
  }
  alignment_mode_missing <- .alignment_mode_missing
  alignment_mode <- match.arg(alignment_mode)
  alignment_provenance <- NULL
  if (!is.null(alignment_features)) {
    .dkge_validate_alignment_features(
      alignment_features, fit = fit, contrast_obj = contrast_obj
    )
    expected_mode <- switch(
      alignment_features$feature_source,
      independent = "independent",
      same_data_residualized = "contrast_orthogonal",
      .dkge_abort(
        "The typed alignment-feature source is unsupported by this transport path.",
        "dkge_alignment_feature_error"
      )
    )
    if (alignment_mode_missing) {
      alignment_mode <- expected_mode
    } else if (!identical(alignment_mode, expected_mode)) {
      .dkge_abort(
        sprintf(
          "`alignment_mode = '%s'` does not match feature source '%s'.",
          alignment_mode, alignment_features$feature_source
        ),
        "dkge_alignment_feature_error"
      )
    }
    if (!is.null(loadings) || !is.null(betas)) {
      .dkge_abort(
        "Typed `alignment_features` cannot be combined with loose loadings or betas.",
        "dkge_alignment_feature_error"
      )
    }
    loadings <- alignment_features$features
    alignment_provenance <- .dkge_alignment_feature_preprocessing(
      alignment_features
    )
  } else if (alignment_mode %in% c("independent", "contrast_orthogonal")) {
    .dkge_abort(
      sprintf(
        "`alignment_mode = '%s'` requires a typed `alignment_features` object.",
        alignment_mode
      ),
      "dkge_alignment_feature_error"
    )
  } else if (identical(alignment_mode, "fold_safe")) {
    if (!is.null(loadings) || !is.null(betas)) {
      .dkge_abort(
        paste0(
          "`alignment_mode = 'fold_safe'` consumes immutable loadings from ",
          "the contrast receipts; do not replace them with `loadings` or `betas`."
        ),
        "dkge_alignment_receipt_error"
      )
    }
    alignment_provenance <- .dkge_fold_receipt_loadings(
      fit, contrast_obj, reference_subject = medoid
    )
    loadings <- alignment_provenance$loadings
  } else if (is.null(loadings)) {
    if (!is.null(betas)) {
      loadings <- dkge_predict_loadings(fit, betas)
    } else if (!is.null(fit$Btil)) {
      loadings <- .dkge_fit_subject_loadings(fit)
    } else {
      stop("Provide betas or pre-computed loadings for descriptive transport.")
    }
    alignment_provenance <- list(
      source = if (!is.null(betas)) "predicted_full_fit" else "full_fit_loadings",
      feature_source = "descriptive_adaptive",
      estimator_source = "descriptive",
      recompute_under_null = FALSE,
      inferentially_eligible = FALSE,
      contract = "Descriptive/transductive functional alignment; not group-inference eligible."
    )
  } else {
    alignment_provenance <- list(
      source = "explicit_descriptive_loadings",
      feature_source = "descriptive_adaptive",
      estimator_source = "descriptive",
      recompute_under_null = FALSE,
      inferentially_eligible = FALSE,
      contract = "Caller-supplied descriptive functional alignment; no inferential receipt."
    )
  }
  S <- length(loadings)
  if (is.null(centroids)) {
    centroids <- fit$centroids %||% fit$input$centroids %||%
      stop("Centroids required for transport; none found in fit or arguments.")
  }
  stopifnot(length(centroids) == S)

  dots <- list(...)
  mapper_supplied <- !is.null(mapper) || !.method_missing || length(dots) > 0L
  if (!is.null(transport_cache) && !mapper_supplied) {
    mapper_spec <- transport_cache$mapper_spec
  } else {
    mapper_spec <- .dkge_resolve_mapper_spec(mapper, method = method[1], dots = dots)
  }

  if (is.null(sizes)) {
    sizes <- lapply(loadings, function(A) rep(1, nrow(A)))
  } else {
    sizes <- Map(function(sz, A) if (is.null(sz)) rep(1, nrow(A)) else sz, sizes, loadings)
  }

  cache_local <- transport_cache
  out <- vector("list", length(contrast_obj$values))
  for (i in seq_along(contrast_obj$values)) {
    tr <- .dkge_transport_to_medoid(mapper_spec,
                                    contrast_obj$values[[i]],
                                    loadings,
                                    centroids,
                                    sizes = sizes,
                                    medoid = medoid,
                                    transport_cache = cache_local,
                                    subject_ids = fit$subject_ids,
                                    preprocessing = alignment_provenance,
                                    reference_selection = reference_selection,
                                    alignment_features = alignment_features)
    out[[i]] <- tr
    if (is.null(cache_local)) {
      cache_local <- tr$fitted_alignment
    }
  }
  names(out) <- names(contrast_obj$values)
  aligned <- .dkge_apply_fitted_alignment(
    contrast_obj$values,
    fitted_alignment = cache_local,
    contrast_obj = contrast_obj,
    contrast_ids = names(out),
    estimand = list(
      population = "analysis subjects represented on one identified reference support",
      contrast_estimator = contrast_obj$method %||% "cross_fitted",
      alignment_status = cache_local$eligibility$status,
      aggregation = "equal_subject",
      value_semantics = cache_local$mapper_spec$params$value_type %||% "intensive"
    ),
    application_context = "dkge_transport_contrasts_to_reference"
  )
  attr(out, "cache") <- cache_local
  attr(out, "fitted_alignment") <- cache_local
  attr(out, "aligned_maps") <- aligned
  attr(out, "alignment_provenance") <- alignment_provenance
  attr(out, "alignment_features_hash") <-
    alignment_features$structural_hash %||% NULL
  attr(out, "alignment_mode") <- alignment_mode
  attr(out, "reference_selection") <- cache_local$reference_selection
  out
}
