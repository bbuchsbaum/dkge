#!/usr/bin/env Rscript

# Deterministic post-processing for the frozen formal court. This script does
# not alter thresholds, seeds, arms, or simulation output.

root <- normalizePath(getwd(), mustWork = TRUE)
outdir <- file.path(root, "data-raw", "functional-alignment-court")
raw <- utils::read.csv(file.path(outdir, "formal-raw.csv"),
                       stringsAsFactors = FALSE)
summary_rows <- utils::read.csv(file.path(outdir, "formal-summary.csv"),
                                stringsAsFactors = FALSE)

same_data <- raw[raw$arm %in% c(
  "l2_fold_loading", "kernel_image_residual_prototype"
), , drop = FALSE]
same_data$excess_rejection <- as.numeric(same_data$reject_fwer) - 0.05

excess_model <- stats::lm(
  excess_rejection ~ latent_span_predictor + arm,
  data = same_data
)
logistic_model <- stats::glm(
  reject_fwer ~ latent_span_predictor + arm,
  family = stats::binomial(),
  data = same_data
)

write_coefficients <- function(model, path) {
  coefficients <- as.data.frame(summary(model)$coefficients)
  coefficients$term <- rownames(coefficients)
  rownames(coefficients) <- NULL
  coefficients <- coefficients[, c("term", setdiff(names(coefficients), "term")),
                               drop = FALSE]
  utils::write.csv(coefficients, path, row.names = FALSE)
}

write_coefficients(
  excess_model,
  file.path(outdir, "latent-span-excess-rejection-model.csv")
)
write_coefficients(
  logistic_model,
  file.path(outdir, "latent-span-logistic-model.csv")
)

required_numeric <- setdiff(
  names(raw)[vapply(raw, is.numeric, logical(1))],
  "exact_elapsed_seconds"
)
audit <- data.frame(
  check = c(
    "raw_rows", "expected_raw_rows", "unique_cohorts",
    "expected_unique_cohorts", "unique_seeds",
    "nonfinite_required_raw", "nonfinite_summary",
    "minimum_cell_arm_convergence", "exact_oracle_gate",
    "negative_control_gate", "candidate_valid"
  ),
  value = c(
    nrow(raw), 21L * 40L * 6L + 6L * 40L,
    length(unique(paste(raw$cell, raw$seed))), 21L * 40L,
    length(unique(raw$seed)),
    sum(!is.finite(as.matrix(raw[, required_numeric, drop = FALSE]))),
    sum(!is.finite(as.matrix(summary_rows[
      , vapply(summary_rows, is.numeric, logical(1)), drop = FALSE
    ]))),
    min(summary_rows$convergence),
    TRUE, FALSE, FALSE
  ),
  stringsAsFactors = FALSE
)
utils::write.csv(audit, file.path(outdir, "formal-audit.csv"), row.names = FALSE)

