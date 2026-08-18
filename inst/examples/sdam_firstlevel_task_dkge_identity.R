# Run the SDAM DKGE analysis with identity design kernels.
#
# This reuses the canonical SDAM analysis script but replaces the structured
# factorial kernels with identity matrices. The output can be compared directly
# against artifacts/sdam-firstlevel-task-dkge to assess how much the structured
# design kernel shapes the latent variables.
#
# Usage (from the dkge repo root):
#   Rscript inst/examples/sdam_firstlevel_task_dkge_identity.R

suppressPackageStartupMessages({
  pkgload::load_all(".", quiet = TRUE, helpers = FALSE, attach_testthat = FALSE)
  library(ggplot2)
  library(neuroim2)
})

source("inst/examples/sdam_firstlevel_task_dkge.R")

.sdam_structured_design_kernel <- make_sdam_design_kernel
.sdam_structured_group_design_kernel <- make_sdam_group_design_kernel
.sdam_structured_cell_mean_bridge_kernel <- make_sdam_cell_mean_bridge_kernel

.sdam_identity_kernel_from_template <- function(template, label = "identity") {
  K0 <- template$K
  K <- diag(1, nrow(K0))
  dimnames(K) <- dimnames(K0)
  info <- template$info %||% list()
  info$metric <- label
  info$identity_kernel <- TRUE
  list(K = K, K_cell = K, info = info)
}

make_sdam_design_kernel <- function(rho_task = 1, rho_measure = 1, rho_int = 1) {
  template <- .sdam_structured_design_kernel(
    rho_task = rho_task,
    rho_measure = rho_measure,
    rho_int = rho_int
  )
  .sdam_identity_kernel_from_template(template, label = "identity_cell_kernel")
}

make_sdam_group_design_kernel <- function(rho_group = 1,
                                          rho_task = 1,
                                          rho_measure = 1,
                                          rho_group_task = 1,
                                          rho_group_measure = 1,
                                          rho_task_measure = 1,
                                          rho_group_task_measure = 1) {
  template <- .sdam_structured_group_design_kernel(
    rho_group = rho_group,
    rho_task = rho_task,
    rho_measure = rho_measure,
    rho_group_task = rho_group_task,
    rho_group_measure = rho_group_measure,
    rho_task_measure = rho_task_measure,
    rho_group_task_measure = rho_group_task_measure
  )
  .sdam_identity_kernel_from_template(template, label = "identity_group_cell_kernel")
}

make_sdam_cell_mean_bridge_kernel <- function(rho_group = 1,
                                              rho_task = 1,
                                              rho_measure = 1,
                                              rho_group_task = 1,
                                              rho_group_measure = 1,
                                              rho_task_measure = 1,
                                              rho_group_task_measure = 1) {
  template <- .sdam_structured_cell_mean_bridge_kernel(
    rho_group = rho_group,
    rho_task = rho_task,
    rho_measure = rho_measure,
    rho_group_task = rho_group_task,
    rho_group_measure = rho_group_measure,
    rho_task_measure = rho_task_measure,
    rho_group_task_measure = rho_group_task_measure
  )
  .sdam_identity_kernel_from_template(template, label = "identity_cell_mean_bridge_kernel")
}

main <- function() {
  output_dir <- file.path("artifacts", "sdam-firstlevel-task-dkge-identity")
  nperm <- as.integer(Sys.getenv("DKGE_SDAM_NPERM", "999"))
  nboot <- as.integer(Sys.getenv("DKGE_SDAM_NBOOT", "200"))
  max_voxels <- {
    raw <- Sys.getenv("DKGE_SDAM_MAX_VOXELS", "")
    if (!nzchar(raw)) Inf else as.numeric(raw)
  }

  analysis <- run_sdam_dkge(
    nperm = nperm,
    nboot = nboot,
    max_voxels = max_voxels
  )

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  analysis_report <- sdam_slim_analysis_for_report(analysis)
  sdam_time(
    "Save identity-kernel report analysis bundle",
    sdam_save_rds(analysis_report, file.path(output_dir, "analysis.rds"))
  )
  sdam_write_q8_sensitivity_tables(analysis_report$group_explicit_sensitivity,
                                   output_dir)

  save_components <- tolower(Sys.getenv("DKGE_SDAM_SAVE_COMPONENTS", "false")) %in%
    c("1", "true", "yes")
  if (save_components) {
    sdam_time(
      "Save compact identity-kernel between-subject result",
      sdam_save_rds(analysis_report$between,
                    file.path(output_dir, "dkge-between-subjects.rds"))
    )
  }

  save_plots <- tolower(Sys.getenv("DKGE_SDAM_SAVE_PLOTS", "true")) %in%
    c("1", "true", "yes")
  if (save_plots) {
    save_sdam_plots(
      analysis$fit,
      analysis$mask,
      analysis$loso,
      analysis$group_diff,
      analysis$between,
      analysis$bsr,
      output_dir = output_dir
    )
  }

  message("\n--- Identity-kernel component summary ---")
  print(summarise_sdam_dkge(analysis$fit))
  message("\n--- Identity-kernel between-subject permutation (~ group) ---")
  print(analysis$between$perm$summary[, c("term", "statistic", "p", "p_adjusted")])
  message("\n--- Identity-kernel q=8 group-explicit validation ---")
  print(analysis$group_explicit$validation$summary[, c("contrast", "estimability",
                                                       "estimate_observed",
                                                       "estimate_completed",
                                                       "sensitivity")])
  message("\n--- Identity-kernel cell-mean bridge permutation ---")
  print(data.frame(
    statistic = analysis_report$cell_mean_bridge$permutation$observed,
    p = analysis_report$cell_mean_bridge$permutation$p,
    B = analysis_report$cell_mean_bridge$permutation$B,
    bootstrap_low = analysis_report$cell_mean_bridge$bootstrap$interval[[1]],
    bootstrap_high = analysis_report$cell_mean_bridge$bootstrap$interval[[2]]
  ))
  invisible(analysis)
}

if (sys.nframe() == 0L) {
  main()
}
