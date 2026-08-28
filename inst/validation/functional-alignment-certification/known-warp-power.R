#!/usr/bin/env Rscript

dktf_power_wilson <- function(x, n, level = 0.95) {
  z <- stats::qnorm(1 - (1 - level) / 2)
  p <- x / n
  denominator <- 1 + z^2 / n
  center <- (p + z^2 / (2 * n)) / denominator
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) /
    denominator
  c(lower = max(0, center - half), upper = min(1, center + half))
}

dktf_power_signs <- function(S) {
  # Representatives modulo the global sign action have the first subject fixed
  # positive. Integer zero is [+,-,...,-]; the final integer is the all-positive
  # identity. Enumerate every representative except that identity because the
  # observed statistic is supplied separately by the add-one test.
  indices <- 0:(2^(S - 1L) - 2L)
  bits <- vapply(indices, function(i) {
    as.integer(intToBits(i))[seq_len(S - 1L)]
  }, integer(S - 1L))
  rbind(1, ifelse(bits == 1L, 1, -1))
}

dktf_power_test <- function(Y, signs) {
  Y <- as.matrix(Y)
  statistic <- function(x) {
    standard_error <- apply(x, 2L, stats::sd) / sqrt(nrow(x)) + 1e-12
    max(abs(colMeans(x) / standard_error))
  }
  observed <- statistic(Y)
  null <- vapply(seq_len(ncol(signs)), function(b) {
    statistic(signs[, b] * Y)
  }, numeric(1))
  p_value <- (1 + sum(null >= observed)) / (ncol(signs) + 1)
  c(statistic = observed, p_value = p_value, reject = p_value <= 0.05)
}

dktf_power_run_one <- function(seed, signal_scale = 0.35,
                               value_noise = 0.8) {
  config <- dktf_config(seed)
  config$heldout_signal_scale <- signal_scale
  config$heldout_value_noise <- value_noise
  result <- dktf_run(config)
  signs <- dktf_power_signs(config$S)
  rows <- lapply(names(result$aligned_rows), function(arm) {
    test <- dktf_power_test(result$aligned_rows[[arm]], signs)
    metric <- result$metrics[result$metrics$arm == arm, , drop = FALSE]
    data.frame(
      seed = seed,
      arm = arm,
      maxT_statistic = unname(test[["statistic"]]),
      p_value = unname(test[["p_value"]]),
      reject = as.logical(test[["reject"]]),
      correlation = metric$correlation,
      rmse = metric$rmse,
      amplitude_ratio = metric$amplitude_ratio,
      positive_peak_error = metric$positive_peak_error,
      negative_peak_error = metric$negative_peak_error,
      latent_error = metric$latent_error,
      point_spread = metric$point_spread,
      solver_converged = isTRUE(result$solver_converged[[arm]]),
      template_converged = isTRUE(result$template$fitting$converged),
      template_iterations = result$template$fitting$iterations,
      runtime_seconds = result$runtime_seconds,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

dktf_power_summarize <- function(raw) {
  split_rows <- split(raw, raw$arm)
  do.call(rbind, lapply(split_rows, function(x) {
    rejected <- sum(x$reject)
    interval <- dktf_power_wilson(rejected, nrow(x))
    data.frame(
      arm = x$arm[[1]],
      n_cohorts = nrow(x),
      rejections = rejected,
      power = rejected / nrow(x),
      power_ci_lower = interval[["lower"]],
      power_ci_upper = interval[["upper"]],
      mean_correlation = mean(x$correlation),
      mean_rmse = mean(x$rmse),
      mean_amplitude_ratio = mean(x$amplitude_ratio),
      mean_positive_peak_error = mean(x$positive_peak_error),
      mean_negative_peak_error = mean(x$negative_peak_error),
      mean_latent_error = mean(x$latent_error),
      mean_point_spread = mean(x$point_spread),
      template_convergence = mean(x$template_converged),
      stringsAsFactors = FALSE
    )
  }))
}
