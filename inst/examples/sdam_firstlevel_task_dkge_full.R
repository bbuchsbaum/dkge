# Run the full SDAM DKGE analysis (no voxel cap) and cache results to disk
# for the Quarto report (`sdam_firstlevel_task_dkge_report.qmd`).
#
# Usage (from the dkge repo root):
#   Rscript inst/examples/sdam_firstlevel_task_dkge_full.R

suppressPackageStartupMessages({
  pkgload::load_all(".", quiet = TRUE, helpers = FALSE, attach_testthat = FALSE)
  library(ggplot2)
  library(neuroim2)
})

source("inst/examples/sdam_firstlevel_task_dkge.R")

output_dir <- file.path("artifacts", "sdam-firstlevel-task-dkge")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(SDAM_SEED)
t0 <- Sys.time()
analysis <- run_sdam_dkge(
  rank        = 3L,
  rho_int     = 1.0,
  nperm       = 999L,
  nboot       = 1000L,
  max_voxels  = Inf,
  contrast_name = "task"
)
t1 <- Sys.time()
message(sprintf("Primary analysis done in %.1f s", as.numeric(t1 - t0, units = "secs")))

# Compute BSR for all three contrasts so the report can show each one.
C <- sdam_contrast_matrix()
bsr_all <- list()
for (nm in c("task", "measure", "task:measure")) {
  bsr_all[[nm]] <- sdam_bsr_voxels(analysis$fit, C[, nm], B = 1000L,
                                   seed = SDAM_SEED + which(colnames(C) == nm))
  message(sprintf("BSR for %-12s done", nm))
}
analysis$bsr_all <- bsr_all

# Slim the report bundle BEFORE saving. Without this the rds carries the
# stratified per-group fits, full LOSO per-subject voxel vectors, and the
# fit's $scores_matrix / $KU / $Chat_sym, ballooning to >10 GB. The slim
# version is a few MB and is what the Quarto report consumes.
analysis_report <- sdam_slim_analysis_for_report(analysis)
saveRDS(analysis_report, file.path(output_dir, "analysis.rds"))
message("Saved slim analysis bundle to ", file.path(output_dir, "analysis.rds"))

cat("\n--- Component summary ---\n")
print(summarise_sdam_dkge(analysis$fit))
cat("\n--- Between-subject permutation (~ group) ---\n")
print(analysis$between$perm$summary[, c("term", "statistic", "p", "p_adjusted")])
cat("\n--- BSR z-map ranges per contrast ---\n")
for (nm in names(bsr_all)) {
  z <- bsr_all[[nm]]$z
  cat(sprintf("%-14s  min=%+.2f  max=%+.2f  |z|>=3 = %d voxels\n",
              nm, min(z, na.rm = TRUE), max(z, na.rm = TRUE),
              sum(abs(z) >= 3, na.rm = TRUE)))
}
