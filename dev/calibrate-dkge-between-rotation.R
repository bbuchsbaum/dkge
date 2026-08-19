#!/usr/bin/env Rscript

# Frozen calibration runner for residual-space rotation inference.
# The design, seed families, gates, and promotion rule are declared in
# data-raw/dkge-between-rotation-plan.md. Do not tune them from these results.

if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("This calibration runner requires devtools.", call. = FALSE)
}

repo_root <- normalizePath(".", mustWork = TRUE)
required <- c(
  "DESCRIPTION",
  "data-raw/dkge-between-rotation-plan.md",
  "tests/testthat/helper-between-rotation.R"
)
if (!all(file.exists(file.path(repo_root, required)))) {
  stop("Run this script from the dkge repository root.", call. = FALSE)
}

devtools::load_all(repo_root, quiet = TRUE)
source(file.path(repo_root, "tests/testthat/helper-between-rotation.R"),
       local = TRUE)

terms_tested <- c("group", "trait", "group:trait")
n_subjects <- 20L
n_features <- 24L
rrr_rank <- 2L
B_null <- 399L
B_power <- 199L

one_result <- function(arm, error, fixture, term, method, batch, replicate,
                       data_seed, resample_seed, B) {
  inference <- dkge_between_permute(
    fixture$fit,
    terms = term,
    method = method,
    B = B,
    seed = resample_seed,
    scope = "global",
    parallel = FALSE
  )
  data.frame(
    arm = arm,
    error = error,
    term = term,
    method = method,
    batch = batch,
    replicate = replicate,
    data_seed = data_seed,
    resample_seed = resample_seed,
    n = n_subjects,
    p_features = n_features,
    rank = rrr_rank,
    B = B,
    p = inference$summary$p,
    reject_005 = inference$summary$p <= 0.05,
    stringsAsFactors = FALSE
  )
}

run_global_family <- function(arm, error, first_base, batches = 5L,
                              per_batch = 200L) {
  rows <- vector("list", batches * per_batch * length(terms_tested))
  out_index <- 0L
  for (batch in seq_len(batches)) {
    base <- first_base + (batch - 1L) * 100000L
    for (replicate in seq_len(per_batch)) {
      data_seed <- base + replicate
      fixture <- dkge_test_rotation_simulate(
        seed = data_seed, n = n_subjects, p = n_features,
        rank = rrr_rank, error = error
      )
      for (term_index in seq_along(terms_tested)) {
        term <- terms_tested[[term_index]]
        resample_seed <- base + 50000L + replicate + term_index * 1000L
        out_index <- out_index + 1L
        rows[[out_index]] <- one_result(
          arm, error, fixture, term, "rotation", batch, replicate,
          data_seed, resample_seed, B_null
        )
      }
    }
    message(sprintf("completed %s batch %d/%d", arm, batch, batches))
  }
  do.call(rbind, rows)
}

run_partial_family <- function(term, first_base, batches = 5L,
                               per_batch = 100L) {
  rows <- vector("list", batches * per_batch)
  out_index <- 0L
  for (batch in seq_len(batches)) {
    base <- first_base + (batch - 1L) * 100000L
    for (replicate in seq_len(per_batch)) {
      data_seed <- base + replicate
      resample_seed <- base + 50000L + replicate
      fixture <- dkge_test_rotation_simulate(
        seed = data_seed, n = n_subjects, p = n_features,
        rank = rrr_rank, error = "gaussian", null_term = term
      )
      out_index <- out_index + 1L
      rows[[out_index]] <- one_result(
        "gaussian_partial_null", "gaussian", fixture, term, "rotation",
        batch, replicate, data_seed, resample_seed, B_null
      )
    }
    message(sprintf("completed partial-null %s batch %d/%d",
                    term, batch, batches))
  }
  do.call(rbind, rows)
}

run_power_family <- function(first_base = 5300000L, total = 300L,
                             per_batch = 100L) {
  rows <- vector("list", total * 2L)
  out_index <- 0L
  for (index in seq_len(total)) {
    batch <- (index - 1L) %/% per_batch + 1L
    replicate <- (index - 1L) %% per_batch + 1L
    base <- first_base + (batch - 1L) * 100000L
    data_seed <- base + replicate
    resample_seed <- base + 50000L + replicate
    fixture <- dkge_test_rotation_simulate(
      seed = data_seed, n = n_subjects, p = n_features,
      rank = rrr_rank, error = "gaussian",
      power_term = "group:trait", effect_size = 1
    )
    for (method in c("rotation", "freedman_lane")) {
      out_index <- out_index + 1L
      rows[[out_index]] <- one_result(
        "strong_interaction_power", "gaussian", fixture, "group:trait",
        method, batch, replicate, data_seed, resample_seed, B_power
      )
    }
    if (replicate == per_batch || index == total) {
      message(sprintf("completed power batch %d/%d",
                      batch, ceiling(total / per_batch)))
    }
  }
  do.call(rbind, rows)
}

wilson_interval <- function(successes, total, level = 0.95) {
  z <- stats::qnorm(1 - (1 - level) / 2)
  proportion <- successes / total
  denominator <- 1 + z^2 / total
  center <- (proportion + z^2 / (2 * total)) / denominator
  half_width <- z * sqrt(
    proportion * (1 - proportion) / total + z^2 / (4 * total^2)
  ) / denominator
  c(lower = center - half_width, upper = center + half_width)
}

summarize_uniform <- function(rows) {
  p <- rows$p
  n_rows <- length(p)
  rejected <- sum(p <= 0.05)
  interval <- wilson_interval(rejected, n_rows)
  quantile_probs <- c(0.05, 0.25, 0.50, 0.75, 0.95)
  empirical_quantiles <- as.numeric(stats::quantile(
    p, probs = quantile_probs, names = FALSE, type = 7
  ))
  ks <- suppressWarnings(stats::ks.test(p, "punif", exact = FALSE))
  batch_sizes <- vapply(split(rows$reject_005, rows$batch), mean, numeric(1))
  arm <- rows$arm[[1L]]
  primary <- arm %in% c("gaussian_global_null", "gaussian_partial_null")
  if (primary) {
    size_gate <- rejected / n_rows >= 0.035 && rejected / n_rows <= 0.065
    wilson_gate <- interval[["lower"]] <= 0.05 &&
      interval[["upper"]] >= 0.05
    ks_gate <- unname(ks$p.value) > 0.01
    quantile_gate <- max(abs(empirical_quantiles - quantile_probs)) < 0.06
    gate_class <- "gaussian_primary"
  } else {
    size_gate <- rejected / n_rows >= 0.025 && rejected / n_rows <= 0.075
    wilson_gate <- interval[["lower"]] <= 0.05
    ks_gate <- unname(ks$p.value) > 0.005
    quantile_gate <- max(abs(empirical_quantiles - quantile_probs)) < 0.08
    gate_class <- "non_gaussian_robustness"
  }
  data.frame(
    arm = arm,
    error = rows$error[[1L]],
    term = rows$term[[1L]],
    method = rows$method[[1L]],
    n_replicates = n_rows,
    B = rows$B[[1L]],
    rejection_005 = rejected / n_rows,
    mcse = sqrt(0.05 * 0.95 / n_rows),
    wilson_lower = interval[["lower"]],
    wilson_upper = interval[["upper"]],
    mean_p = mean(p),
    median_p = stats::median(p),
    ks_D = unname(ks$statistic),
    ks_p = unname(ks$p.value),
    batch_rejection_min = min(batch_sizes),
    batch_rejection_max = max(batch_sizes),
    q05 = empirical_quantiles[[1L]],
    q25 = empirical_quantiles[[2L]],
    q50 = empirical_quantiles[[3L]],
    q75 = empirical_quantiles[[4L]],
    q95 = empirical_quantiles[[5L]],
    max_quantile_deviation = max(abs(empirical_quantiles - quantile_probs)),
    gate_class = gate_class,
    size_gate = size_gate,
    wilson_gate = wilson_gate,
    ks_gate = ks_gate,
    quantile_gate = quantile_gate,
    gate_pass = size_gate && wilson_gate && ks_gate && quantile_gate,
    stringsAsFactors = FALSE
  )
}

summarize_power <- function(rows, rotation_power, fl_power) {
  n_rows <- nrow(rows)
  rejected <- sum(rows$reject_005)
  interval <- wilson_interval(rejected, n_rows)
  method <- rows$method[[1L]]
  gate <- if (identical(method, "rotation")) {
    rotation_power >= 0.75 && rotation_power >= fl_power - 0.12
  } else {
    NA
  }
  data.frame(
    arm = rows$arm[[1L]],
    error = rows$error[[1L]],
    term = rows$term[[1L]],
    method = method,
    n_replicates = n_rows,
    B = rows$B[[1L]],
    rejection_005 = rejected / n_rows,
    mcse = sqrt((rejected / n_rows) * (1 - rejected / n_rows) / n_rows),
    wilson_lower = interval[["lower"]],
    wilson_upper = interval[["upper"]],
    mean_p = mean(rows$p),
    median_p = stats::median(rows$p),
    ks_D = NA_real_,
    ks_p = NA_real_,
    batch_rejection_min = min(vapply(
      split(rows$reject_005, rows$batch), mean, numeric(1)
    )),
    batch_rejection_max = max(vapply(
      split(rows$reject_005, rows$batch), mean, numeric(1)
    )),
    q05 = unname(stats::quantile(rows$p, 0.05)),
    q25 = unname(stats::quantile(rows$p, 0.25)),
    q50 = unname(stats::quantile(rows$p, 0.50)),
    q75 = unname(stats::quantile(rows$p, 0.75)),
    q95 = unname(stats::quantile(rows$p, 0.95)),
    max_quantile_deviation = NA_real_,
    gate_class = "power",
    size_gate = NA,
    wilson_gate = NA,
    ks_gate = NA,
    quantile_gate = NA,
    gate_pass = gate,
    stringsAsFactors = FALSE
  )
}

started_at <- format(Sys.time(), tz = "UTC", usetz = TRUE)

replicates <- rbind(
  run_global_family("gaussian_global_null", "gaussian", 1700000L),
  run_partial_family("group", 2300000L),
  run_partial_family("trait", 2900000L),
  run_partial_family("group:trait", 3500000L),
  run_global_family("t5_global_null", "t5", 4100000L,
                    batches = 5L, per_batch = 100L),
  run_global_family("skew_global_null", "skew", 4700000L,
                    batches = 5L, per_batch = 100L),
  run_power_family()
)

expected_rows <- 3000L + 1500L + 1500L + 1500L + 600L
stopifnot(nrow(replicates) == expected_rows)
stopifnot(!anyNA(replicates$p), all(replicates$p > 0),
          all(replicates$p <= 1))

null_rows <- replicates[replicates$arm != "strong_interaction_power", ,
                        drop = FALSE]
null_groups <- interaction(null_rows$arm, null_rows$error, null_rows$term,
                           null_rows$method, drop = TRUE)
summary_rows <- do.call(rbind, lapply(split(null_rows, null_groups),
                                     summarize_uniform))

power_rows <- replicates[replicates$arm == "strong_interaction_power", ,
                         drop = FALSE]
power_by_method <- vapply(split(power_rows$reject_005, power_rows$method),
                          mean, numeric(1))
power_summary <- do.call(rbind, lapply(split(power_rows, power_rows$method),
                                      summarize_power,
                                      rotation_power = power_by_method[["rotation"]],
                                      fl_power = power_by_method[["freedman_lane"]]))
summary_rows <- rbind(summary_rows, power_summary)
rownames(summary_rows) <- NULL
summary_rows <- summary_rows[order(summary_rows$arm, summary_rows$term,
                                   summary_rows$method), , drop = FALSE]

sha256_file <- function(path) {
  output <- system2("shasum", c("-a", "256", path), stdout = TRUE)
  sub("[[:space:]].*$", "", output[[1L]])
}

plan_path <- file.path(repo_root, "data-raw/dkge-between-rotation-plan.md")
runner_path <- file.path(repo_root, "dev/calibrate-dkge-between-rotation.R")
metadata <- data.frame(
  key = c(
    "started_at_utc", "completed_at_utc", "plan_sha256", "runner_sha256",
    "git_head", "R_version", "dkge_version", "platform",
    "gaussian_primary_all_pass", "non_gaussian_all_pass", "power_pass"
  ),
  value = c(
    started_at,
    format(Sys.time(), tz = "UTC", usetz = TRUE),
    sha256_file(plan_path),
    sha256_file(runner_path),
    system2("git", c("rev-parse", "HEAD"), stdout = TRUE)[[1L]],
    R.version.string,
    as.character(utils::packageVersion("dkge")),
    R.version$platform,
    as.character(all(summary_rows$gate_pass[
      summary_rows$gate_class == "gaussian_primary"
    ])),
    as.character(all(summary_rows$gate_pass[
      summary_rows$gate_class == "non_gaussian_robustness"
    ])),
    as.character(summary_rows$gate_pass[
      summary_rows$gate_class == "power" & summary_rows$method == "rotation"
    ])
  ),
  stringsAsFactors = FALSE
)

replicate_path <- file.path(
  repo_root, "inst/extdata/dkge-between-rotation-replicates.csv"
)
summary_path <- file.path(
  repo_root, "inst/extdata/dkge-between-rotation-summary.csv"
)
metadata_path <- file.path(
  repo_root, "inst/extdata/dkge-between-rotation-metadata.csv"
)
utils::write.csv(replicates, replicate_path, row.names = FALSE)
utils::write.csv(summary_rows, summary_path, row.names = FALSE)
utils::write.csv(metadata, metadata_path, row.names = FALSE)

print(summary_rows[, c("arm", "term", "method", "n_replicates",
                       "rejection_005", "wilson_lower", "wilson_upper",
                       "ks_p", "max_quantile_deviation", "gate_pass")],
      row.names = FALSE)
message("wrote ", replicate_path)
message("wrote ", summary_path)
message("wrote ", metadata_path)
