#!/usr/bin/env Rscript

if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("The simulation runner requires the suggested package 'devtools'.")
}
devtools::load_all(quiet = TRUE)

`%or%` <- function(x, y) if (is.null(x)) y else x

rep_count <- suppressWarnings(as.integer(Sys.getenv("DKGE_SIM_REPS", "20")))
if (!is.finite(rep_count) || rep_count < 1L || rep_count > 20L) {
  stop("DKGE_SIM_REPS must be an integer from 1 through 20.")
}
replicate_seeds <- 73101:73120
replicate_seeds <- replicate_seeds[seq_len(rep_count)]

scenarios <- data.frame(
  scenario = seq_len(10),
  signal = c("null", "null", "ordinal", "ordinal", "nonlinear",
             "interaction", "low_rank", "ordinal", "nonlinear", "low_rank"),
  noise = c("iid", "heteroscedastic", "iid", "heteroscedastic", "ar",
            "ar", "unequal_precision", "unequal_precision",
            "heteroscedastic", "ar"),
  missingness = c("mcar", "response", "mcar", "response", "mcar",
                  "response", "mcar", "response", "response", "mcar"),
  stringsAsFactors = FALSE
)

grid <- dkge_effect_grid(list(
  A = list(L = 3, type = "nominal", levels = paste0("A", 1:3)),
  B = list(L = 5, type = "nominal", levels = paste0("B", 1:5)),
  response = list(L = 4, type = "ordinal", levels = as.character(1:4), l = 1)
))
kernel <- design_kernel(
  grid,
  terms = list("A", "B", "response", c("A", "B"),
               c("A", "response"), c("B", "response"),
               c("A", "B", "response")),
  basis = "cell", normalize = "max_diag"
)$K
cells <- grid$cells
effect_ids <- grid$cell_labels
q <- length(effect_ids)
P <- 12L
S <- 18L

center_scale <- function(x) {
  x <- as.numeric(x)
  x <- x - mean(x)
  x / sqrt(sum(x^2))
}
response_linear <- center_scale(as.numeric(cells$response))
response_curve <- center_scale((as.numeric(cells$response) - 2.5)^2)
a_score <- center_scale(as.numeric(factor(cells$A, levels = paste0("A", 1:3))))
b_score <- center_scale(as.numeric(factor(cells$B, levels = paste0("B", 1:5))))
interaction_score <- center_scale(response_linear * (a_score + 0.7 * b_score))

spatial <- qr.Q(qr(cbind(
  sin(seq(0, 2 * pi, length.out = P)),
  cos(seq(0, 3 * pi, length.out = P))
)))[, 1:2, drop = FALSE]

signal_basis <- function(signal) {
  switch(
    signal,
    null = matrix(0, q, 1),
    ordinal = matrix(response_linear, q, 1),
    nonlinear = matrix(response_curve, q, 1),
    interaction = matrix(interaction_score, q, 1),
    low_rank = cbind(response_linear, center_scale(response_curve + a_score))
  )
}

block_correlation <- function(run, rho) {
  Tn <- length(run)
  out <- diag(Tn)
  if (rho == 0) return(out)
  for (label in unique(run)) {
    idx <- which(run == label)
    out[idx, idx] <- outer(seq_along(idx), seq_along(idx),
                           function(i, j) rho^abs(i - j))
  }
  out
}

make_counts <- function(noise, missingness) {
  counts <- matrix(0L, S, q, dimnames = list(paste0("s", 1:S), effect_ids))
  response_index <- as.numeric(cells$response)
  for (s in seq_len(S)) {
    exposure <- exp(rnorm(1, sd = if (noise == "unequal_precision") 0.65 else 0.35))
    lambda <- 1.8 * exposure * c(1.3, 1.1, 0.9, 0.65)[response_index]
    drop_probability <- if (missingness == "mcar") {
      rep(0.14, q)
    } else {
      0.04 + 0.08 * (response_index - 1)
    }
    present <- runif(q) > drop_probability
    counts[s, present] <- 4L + rpois(sum(present), lambda[present])
  }
  absent <- which(colSums(counts) == 0)
  if (length(absent)) counts[1L, absent] <- 4L
  counts
}

make_dataset <- function(scenario, seed) {
  set.seed(seed + 1000L * scenario$scenario)
  phi <- signal_basis(scenario$signal)
  true_rank <- if (scenario$signal == "low_rank") 2L else 1L
  signal_amplitude <- if (scenario$signal == "null") 0 else 1.25
  counts <- make_counts(scenario$noise, scenario$missingness)
  ols_subjects <- vector("list", S)
  gls_subjects <- vector("list", S)
  truth <- vector("list", S)

  for (s in seq_len(S)) {
    local <- which(counts[s, ] > 0)
    local_ids <- effect_ids[local]
    coefficient <- 1 + rnorm(ncol(phi), sd = 0.12)
    B_true <- signal_amplitude * phi %*% diag(coefficient, ncol(phi)) %*%
      t(spatial[, seq_len(ncol(phi)), drop = FALSE])
    rownames(B_true) <- effect_ids
    colnames(B_true) <- paste0("v", seq_len(P))
    truth[[s]] <- B_true

    trial_cell <- unlist(Map(function(cell, n) rep(cell, n),
                             local, counts[s, local]), use.names = FALSE)
    trial_run <- unlist(lapply(counts[s, local], function(n) {
      rep(seq_len(4), length.out = n)
    }), use.names = FALSE)
    order_index <- order(trial_run, runif(length(trial_run)))
    trial_cell <- trial_cell[order_index]
    trial_run <- trial_run[order_index]
    local_column <- match(trial_cell, local)
    X <- matrix(0, length(trial_cell), length(local),
                dimnames = list(NULL, local_ids))
    X[cbind(seq_along(local_column), local_column)] <- 1

    response_trial <- as.numeric(cells$response[trial_cell])
    subject_scale <- if (scenario$noise == "unequal_precision") {
      exp(seq(log(0.55), log(2.2), length.out = S)[s])
    } else {
      1
    }
    trial_sd <- switch(
      scenario$noise,
      iid = rep(1, nrow(X)),
      ar = rep(1, nrow(X)),
      heteroscedastic = subject_scale * (0.65 + 0.28 * response_trial),
      unequal_precision = subject_scale * (0.75 + 0.18 * response_trial)
    )
    feature_variance <- if (scenario$noise == "unequal_precision") {
      seq(0.5, 2, length.out = P)
    } else {
      rep(1, P)
    }
    rho <- if (scenario$noise == "ar") 0.60 else 0
    correlation <- block_correlation(trial_run, rho)
    covariance <- correlation * tcrossprod(trial_sd)
    error <- t(chol(covariance)) %*%
      matrix(rnorm(nrow(X) * P), nrow(X), P)
    error <- sweep(error, 2, sqrt(feature_variance), "*")
    Y <- X %*% B_true[local, , drop = FALSE] + error
    colnames(Y) <- colnames(B_true)

    ols <- suppressWarnings(dkge_trial_subject(
      Y, X, id = paste0("s", s), split = "run",
      run_labels = paste0("run", trial_run),
      effect_precision = "split_half"
    ))
    gls <- suppressWarnings(dkge_trial_subject(
      Y, X, id = paste0("s", s), trial_covariance = covariance
    ))
    precision <- 1 / pmax(
      diag(gls$effect_noise_cov) * mean(gls$residual_variance), 1e-10
    )
    names(precision) <- local_ids
    ols$known_precision <- precision
    ols_subjects[[s]] <- ols
    gls_subjects[[s]] <- gls
  }
  list(
    ols = ols_subjects, gls = gls_subjects, truth = truth,
    counts = counts, true_rank = true_rank
  )
}

fit_methods <- function(dataset) {
  explicit <- lapply(dataset$ols, function(subject) {
    subject$effect_precision <- subject$known_precision
    subject
  })
  common <- list(K = kernel, rank = dataset$true_rank, w_method = "none",
                 missingness = "rescale")
  call_fit <- function(subjects, ...) {
    suppressWarnings(do.call(
      dkge_fit, c(list(data = dkge_data(subjects)), common, list(...))
    ))
  }
  list(
    legacy = call_fit(dataset$ols),
    count = call_fit(dataset$ols,
                     effect_weights = dkge_effect_weights("count")),
    explicit_precision = call_fit(
      explicit, effect_weights = dkge_effect_weights("precision")
    ),
    random_effects = call_fit(
      dataset$ols,
      effect_weights = dkge_effect_weights(
        "random_effects", within = "count", max_ratio = 10
      )
    ),
    iid_analytic = call_fit(dataset$ols, debias = "analytic"),
    covariance_analytic = call_fit(dataset$gls, debias = "analytic"),
    split_half = call_fit(
      dataset$ols, effect_weights = dkge_effect_weights("precision"),
      debias = "split_half"
    )
  )
}

truth_pool <- function(fit, truth) {
  masks <- getFromNamespace(".dkge_obs_masks_from_provenance", "dkge")(
    fit$provenance, fit$subject_ids, q
  )
  active <- !identical(fit$effect_weight_spec$method %or% "none", "none")
  numerator <- matrix(0, q, q)
  denominator <- matrix(0, q, q)
  for (s in seq_len(S)) {
    mask <- as.numeric(masks[[s]])
    moment <- tcrossprod(truth[[s]])
    if (active) {
      p <- pmax(as.numeric(fit$effect_precision[[s]]), 0)
      h <- fit$weights[[s]] * tcrossprod(sqrt(p)) * tcrossprod(mask)
    } else {
      h <- fit$weights[[s]] * tcrossprod(mask)
    }
    numerator <- numerator + h * moment
    denominator <- denominator + h
  }
  out <- matrix(0, q, q)
  valid <- denominator > 0
  out[valid] <- numerator[valid] / denominator[valid]
  if (active) out <- sum(fit$weights) * out
  (out + t(out)) / 2
}

principal_angle <- function(fit, target_moment, rank) {
  target_chat <- fit$Khalf %*% t(fit$R) %*% target_moment %*%
    fit$R %*% fit$Khalf
  target_eigen <- eigen((target_chat + t(target_chat)) / 2, symmetric = TRUE)
  target_basis <- target_eigen$vectors[, seq_len(rank), drop = FALSE]
  estimated <- fit$Khalf %*% fit$U
  if (ncol(estimated) < rank) return(90)
  singular <- svd(crossprod(target_basis, estimated[, seq_len(rank), drop = FALSE]),
                  nu = 0, nv = 0)$d
  max(acos(pmin(1, pmax(0, singular)))) * 180 / pi
}

subject_statistic <- function(moment, contrast) {
  as.numeric(crossprod(contrast, moment %*% contrast))
}

method_metrics <- function(fit, dataset, method, scenario, seed) {
  target <- truth_pool(fit, dataset$truth)
  nonnull <- scenario$signal != "null"
  angle <- if (nonnull) principal_angle(fit, target, dataset$true_rank) else NA_real_
  leading <- if (ncol(fit$U)) {
    fit$Khalf %*% fit$U[, 1, drop = FALSE]
  } else {
    matrix(0, q, 1)
  }
  inverse_count <- colMeans(ifelse(dataset$counts > 0, 1 / dataset$counts, NA),
                            na.rm = TRUE)
  inverse_correlation <- if (sd(inverse_count) > 0) {
    abs(cor(as.numeric(leading)^2, inverse_count))
  } else {
    NA_real_
  }

  contrast <- response_linear / sqrt(sum(response_linear^2))
  target_scale <- max(norm(target, "F"), 1e-8)
  contrast_error <- if (nonnull) {
    subject_statistic(fit$effect_moment - target, contrast) / target_scale
  } else {
    NA_real_
  }

  true_subject <- vapply(dataset$truth, function(B) {
    subject_statistic(tcrossprod(B), contrast)
  }, numeric(1))
  estimated_subject <- vapply(fit$effect_moments, subject_statistic,
                              numeric(1), contrast = contrast)
  set.seed(seed + 700000L + scenario$scenario)
  boot <- replicate(399, mean(sample(estimated_subject, S, replace = TRUE)))
  interval <- stats::quantile(boot, c(0.025, 0.975), names = FALSE)
  coverage <- interval[[1]] <= mean(true_subject) && interval[[2]] >= mean(true_subject)

  linear_subject <- vapply(seq_len(S), function(s) {
    B <- fit$Braw[[s]]
    as.numeric(crossprod(contrast, B) %*% spatial[, 1])
  }, numeric(1))
  set.seed(seed + 900000L + scenario$scenario)
  null_p <- dkge_signflip_maxT(matrix(linear_subject, ncol = 1), B = 999)$p[[1]]

  precision <- fit$effect_precision
  influence <- vapply(seq_len(q), function(effect) {
    weight <- vapply(seq_len(S), function(s) {
      fit$weights[[s]] * pmax(as.numeric(precision[[s]])[[effect]], 0)
    }, numeric(1))
    if (sum(weight) > 0) max(weight) / sum(weight) else NA_real_
  }, numeric(1))

  data.frame(
    scenario = scenario$scenario,
    seed = seed,
    signal = scenario$signal,
    noise = scenario$noise,
    missingness = scenario$missingness,
    method = method,
    subspace_angle_deg = angle,
    inverse_count_correlation = inverse_correlation,
    contrast_error = contrast_error,
    interval_covered = as.numeric(coverage),
    null_reject = as.numeric(scenario$signal == "null" && null_p <= 0.05),
    null_p = null_p,
    high_count_influence = mean(influence, na.rm = TRUE),
    empty_cell_fraction = mean(dataset$counts == 0),
    stringsAsFactors = FALSE
  )
}

replicate_rows <- vector("list", nrow(scenarios) * length(replicate_seeds))
cursor <- 0L
for (scenario_index in seq_len(nrow(scenarios))) {
  scenario <- scenarios[scenario_index, , drop = FALSE]
  for (seed in replicate_seeds) {
    cursor <- cursor + 1L
    message(sprintf("scenario %d/10 seed %d", scenario$scenario, seed))
    dataset <- make_dataset(scenario, seed)
    fits <- fit_methods(dataset)
    replicate_rows[[cursor]] <- do.call(rbind, Map(
      function(fit, method) method_metrics(fit, dataset, method, scenario, seed),
      fits, names(fits)
    ))
  }
}
replicates <- do.call(rbind, replicate_rows)

split_groups <- split(replicates, interaction(replicates$scenario,
                                               replicates$method, drop = TRUE))
summary_rows <- lapply(split_groups, function(x) {
  data.frame(
    row_type = "scenario_method",
    scenario = x$scenario[[1]], signal = x$signal[[1]], noise = x$noise[[1]],
    missingness = x$missingness[[1]], method = x$method[[1]],
    n_replicates = nrow(x),
    median_subspace_angle_deg = stats::median(x$subspace_angle_deg, na.rm = TRUE),
    median_inverse_count_correlation = stats::median(x$inverse_count_correlation,
                                                     na.rm = TRUE),
    contrast_bias = mean(x$contrast_error, na.rm = TRUE),
    contrast_rmse = sqrt(mean(x$contrast_error^2, na.rm = TRUE)),
    interval_coverage = mean(x$interval_covered),
    false_positive_rate = if (x$signal[[1]] == "null") mean(x$null_reject) else NA_real_,
    mean_high_count_influence = mean(x$high_count_influence),
    mean_empty_cell_fraction = mean(x$empty_cell_fraction),
    gate_id = NA_character_, passed = NA, observed = NA_real_,
    threshold = NA_character_, stringsAsFactors = FALSE
  )
})
summary_table <- do.call(rbind, summary_rows)

method_value <- function(data, method, column) {
  values <- data[data$method == method, column]
  stats::median(values, na.rm = TRUE)
}
nonnull <- replicates$signal != "null"
hard_noise <- replicates$noise %in% c("ar", "unequal_precision") & nonnull
high_imbalance <- (replicates$noise == "unequal_precision" |
                     replicates$missingness == "response") & nonnull
iid_nonnull <- replicates$noise == "iid" & nonnull

gate <- function(id, passed, observed, threshold) {
  data.frame(
    row_type = "gate", scenario = NA_integer_, signal = NA_character_,
    noise = NA_character_, missingness = NA_character_, method = NA_character_,
    n_replicates = nrow(replicates), median_subspace_angle_deg = NA_real_,
    median_inverse_count_correlation = NA_real_, contrast_bias = NA_real_,
    contrast_rmse = NA_real_, interval_coverage = NA_real_,
    false_positive_rate = NA_real_, mean_high_count_influence = NA_real_,
    mean_empty_cell_fraction = NA_real_, gate_id = id, passed = passed,
    observed = observed, threshold = threshold, stringsAsFactors = FALSE
  )
}

legacy_angle <- method_value(replicates[hard_noise, ], "legacy", "subspace_angle_deg")
cov_angle <- method_value(replicates[hard_noise, ], "covariance_analytic", "subspace_angle_deg")
legacy_rmse <- sqrt(mean(replicates$contrast_error[hard_noise & replicates$method == "legacy"]^2,
                         na.rm = TRUE))
cov_rmse <- sqrt(mean(replicates$contrast_error[hard_noise &
                                                 replicates$method == "covariance_analytic"]^2,
                      na.rm = TRUE))
gates <- list(gate("G1_angle", cov_angle <= 0.90 * legacy_angle,
                   cov_angle / legacy_angle, "ratio <= 0.90"),
              gate("G1_rmse", cov_rmse <= 0.90 * legacy_rmse,
                   cov_rmse / legacy_rmse, "ratio <= 0.90"))

for (method in c("explicit_precision", "random_effects",
                 "covariance_analytic", "split_half")) {
  legacy <- method_value(replicates[nonnull, ], "legacy",
                         "inverse_count_correlation")
  corrected <- method_value(replicates[nonnull, ], method,
                            "inverse_count_correlation")
  gates[[length(gates) + 1L]] <- gate(
    paste0("G2_", method), corrected <= 0.10 || corrected <= legacy - 0.05,
    corrected - legacy, "difference <= -0.05 or corrected <= 0.10"
  )
}
for (method in c("covariance_analytic", "split_half")) {
  legacy <- abs(method_value(replicates[nonnull, ], "legacy", "contrast_error"))
  corrected <- abs(method_value(replicates[nonnull, ], method, "contrast_error"))
  gates[[length(gates) + 1L]] <- gate(
    paste0("G3_", method), corrected <= legacy,
    corrected - legacy, "absolute median bias difference <= 0"
  )
  null_rows <- replicates$signal == "null" & replicates$method == method
  fpr <- mean(replicates$null_reject[null_rows])
  gates[[length(gates) + 1L]] <- gate(
    paste0("G4_", method), fpr >= 0.01 && fpr <= 0.12,
    fpr, "0.01 <= FPR <= 0.12"
  )
  coverage <- mean(replicates$interval_covered[nonnull & replicates$method == method])
  gates[[length(gates) + 1L]] <- gate(
    paste0("G5_", method), coverage >= 0.80 && coverage <= 1,
    coverage, "0.80 <= coverage <= 1.00"
  )
}
count_influence <- method_value(replicates[high_imbalance, ], "count",
                                "high_count_influence")
random_influence <- method_value(replicates[high_imbalance, ], "random_effects",
                                 "high_count_influence")
gates[[length(gates) + 1L]] <- gate(
  "G6_random_effects_influence",
  random_influence <= 0.90 * count_influence,
  random_influence / count_influence, "ratio <= 0.90"
)
iid_difference <- abs(
  method_value(replicates[iid_nonnull, ], "iid_analytic", "subspace_angle_deg") -
    method_value(replicates[iid_nonnull, ], "covariance_analytic",
                 "subspace_angle_deg")
)
gates[[length(gates) + 1L]] <- gate(
  "G7_iid_equivalence", iid_difference <= 5, iid_difference,
  "absolute difference <= 5 degrees"
)
summary_table <- rbind(summary_table, do.call(rbind, gates))

dir.create("inst/extdata", recursive = TRUE, showWarnings = FALSE)
utils::write.csv(replicates,
                 "inst/extdata/dkge-unbalanced-trialwise-replicates.csv",
                 row.names = FALSE)
utils::write.csv(summary_table,
                 "inst/extdata/dkge-unbalanced-trialwise-summary.csv",
                 row.names = FALSE)
message(sprintf("completed %d replicate-method rows; %d/%d gates passed",
                nrow(replicates), sum(vapply(gates, `[[`, logical(1), "passed")),
                length(gates)))
