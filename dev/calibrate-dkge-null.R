#!/usr/bin/env Rscript

if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("The calibration runner requires the suggested package 'devtools'.")
}

devtools::load_all(quiet = TRUE)
source("tests/testthat/helper-between-null.R", local = TRUE)
source("tests/testthat/helper-weights-null.R", local = TRUE)

started <- Sys.time()
alpha <- 0.05
quantile_probs <- c(0.10, 0.25, 0.50, 0.75, 0.90)

between_primary_seeds <- c(800L, 100800L, 200800L, 300800L, 400800L)
between_n40_seeds <- c(600800L, 700800L, 800800L, 900800L, 1000800L)
between_n80_seeds <- c(1100800L, 1200800L, 1300800L, 1400800L, 1500800L)
adaptive_seed_starts <- c(none = 101L, kenergy = 20101L,
                          kenergy_prec = 40101L)

run_between_batches <- function(setting, n, reps_per_batch, seeds) {
  rows <- vector("list", length(seeds))
  for (batch in seq_along(seeds)) {
    message("between ", setting, " batch ", batch, "/", length(seeds))
    pvals <- dkge_test_between_null_uniformity(
      nrep = reps_per_batch,
      n = n,
      p = 24,
      rank = 2,
      B = 399,
      seed = seeds[[batch]]
    )
    long <- do.call(rbind, lapply(names(pvals), function(term) {
      data.frame(
        family = "between_freedman_lane",
        setting = setting,
        estimand = term,
        n = n,
        batch = batch,
        base_seed = seeds[[batch]],
        replicate = seq_len(nrow(pvals)),
        p_value = pvals[[term]],
        legacy_ks_threshold = 0.05,
        legacy_quantile_threshold = 0.09,
        stringsAsFactors = FALSE
      )
    }))
    rows[[batch]] <- long
  }
  do.call(rbind, rows)
}

run_adaptive_batches <- function(method, reps_per_batch = 200L) {
  seeds <- adaptive_seed_starts[[method]] + 100000L * 0:4
  rows <- vector("list", length(seeds))
  for (batch in seq_along(seeds)) {
    message("adaptive ", method, " batch ", batch, "/", length(seeds))
    pvals <- dkge_test_null_uniformity(
      nrep = reps_per_batch,
      nsub = 16,
      V = 128,
      adapt = method,
      n_perm = 499,
      seed = seeds[[batch]]
    )
    rows[[batch]] <- data.frame(
      family = "adaptive_loso_signflip",
      setting = "n16_v128",
      estimand = method,
      n = 16L,
      batch = batch,
      base_seed = seeds[[batch]],
      replicate = seq_along(pvals),
      p_value = pvals,
      legacy_ks_threshold = 0.05,
      legacy_quantile_threshold = 0.08,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

wilson_interval <- function(successes, total, level = 0.95) {
  z <- stats::qnorm(1 - (1 - level) / 2)
  phat <- successes / total
  denom <- 1 + z^2 / total
  center <- (phat + z^2 / (2 * total)) / denom
  half <- z * sqrt(phat * (1 - phat) / total + z^2 / (4 * total^2)) / denom
  c(lower = center - half, upper = center + half)
}

summarize_group <- function(x) {
  p <- x$p_value
  total <- length(p)
  successes <- sum(p <= alpha)
  size <- successes / total
  interval <- wilson_interval(successes, total)
  ks <- suppressWarnings(stats::ks.test(p, "punif"))
  empirical_quantiles <- as.numeric(stats::quantile(
    p, probs = quantile_probs, names = FALSE, type = 7
  ))
  max_deviation <- max(abs(empirical_quantiles - quantile_probs))
  batch_sizes <- vapply(split(p, x$batch), function(z) mean(z <= alpha), numeric(1))
  legacy_gate <- unname(ks$p.value) > x$legacy_ks_threshold[[1]] &&
    max_deviation < x$legacy_quantile_threshold[[1]]
  classification <- if (interval[["lower"]] > alpha) {
    "anti_conservative"
  } else if (interval[["upper"]] < alpha) {
    "conservative"
  } else if (legacy_gate) {
    "calibrated"
  } else {
    "inconclusive"
  }
  data.frame(
    family = x$family[[1]],
    setting = x$setting[[1]],
    estimand = x$estimand[[1]],
    n = x$n[[1]],
    replicates = total,
    permutations = if (x$family[[1]] == "between_freedman_lane") 399L else 499L,
    empirical_size_0.05 = size,
    size_mcse = sqrt(size * (1 - size) / total),
    size_wilson_lower = interval[["lower"]],
    size_wilson_upper = interval[["upper"]],
    batch_size_min = min(batch_sizes),
    batch_size_max = max(batch_sizes),
    mean_p = mean(p),
    median_p = stats::median(p),
    ks_d = unname(ks$statistic),
    ks_p = unname(ks$p.value),
    q10 = empirical_quantiles[[1]],
    q25 = empirical_quantiles[[2]],
    q50 = empirical_quantiles[[3]],
    q75 = empirical_quantiles[[4]],
    q90 = empirical_quantiles[[5]],
    max_quantile_deviation = max_deviation,
    legacy_ks_threshold = x$legacy_ks_threshold[[1]],
    legacy_quantile_threshold = x$legacy_quantile_threshold[[1]],
    legacy_gate_passed = legacy_gate,
    classification = classification,
    stringsAsFactors = FALSE
  )
}

replicates <- rbind(
  run_between_batches("primary_n20", 20L, 200L, between_primary_seeds),
  run_between_batches("diagnostic_n40", 40L, 100L, between_n40_seeds),
  run_between_batches("diagnostic_n80", 80L, 100L, between_n80_seeds),
  run_adaptive_batches("none"),
  run_adaptive_batches("kenergy"),
  run_adaptive_batches("kenergy_prec")
)

group_key <- interaction(replicates$family, replicates$setting,
                         replicates$estimand, drop = TRUE, lex.order = TRUE)
summary_rows <- lapply(split(replicates, group_key), summarize_group)
summary <- do.call(rbind, summary_rows)
rownames(summary) <- NULL

dir.create("inst/extdata", recursive = TRUE, showWarnings = FALSE)
utils::write.csv(
  replicates,
  "inst/extdata/dkge-null-calibration-replicates.csv",
  row.names = FALSE
)
utils::write.csv(
  summary,
  "inst/extdata/dkge-null-calibration-summary.csv",
  row.names = FALSE
)

plan_hash <- tryCatch(
  system2(
    "shasum",
    c("-a", "256", "data-raw/dkge-null-calibration-plan.md"),
    stdout = TRUE
  ),
  error = function(e) NA_character_
)
plan_hash <- if (length(plan_hash) == 1L && !is.na(plan_hash)) {
  sub("[[:space:]].*$", "", plan_hash)
} else {
  NA_character_
}
git_head <- tryCatch(
  system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
  error = function(e) NA_character_
)
metadata <- data.frame(
  calibration_version = "null-calibration-v1",
  started_utc = format(started, tz = "UTC", usetz = TRUE),
  completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
  plan_sha256 = plan_hash,
  git_head = git_head[[1]],
  package_version = as.character(utils::packageVersion("dkge")),
  r_version = R.version.string,
  platform = R.version$platform,
  stringsAsFactors = FALSE
)
utils::write.csv(
  metadata,
  "inst/extdata/dkge-null-calibration-metadata.csv",
  row.names = FALSE
)

print(summary, row.names = FALSE)
message("calibration completed in ", round(metadata$elapsed_seconds, 1), " seconds")

