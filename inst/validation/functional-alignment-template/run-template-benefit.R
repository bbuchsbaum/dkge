#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)
} else {
  normalizePath("inst/validation/functional-alignment-template/run-template-benefit.R",
                mustWork = TRUE)
}
root <- normalizePath(file.path(dirname(script), "..", "..", ".."),
                      mustWork = TRUE)
source(file.path(
  root, "inst", "validation", "functional-alignment-certification",
  "evidence-utils.R"
), local = TRUE)
source(file.path(root, "inst", "validation", "functional-alignment-template",
                 "template-benefit.R"))

source_loader <- dkfa_load_source_package(root, quiet = TRUE)

result <- dktf_run(dktf_config(911L))
output <- file.path(root, "data-raw", "functional-alignment-template")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(result$metrics,
                 file.path(output, "known-warp-metrics.csv"),
                 row.names = FALSE)
jsonlite::write_json(
  list(
    protocol = result$protocol,
    config = result$config,
    gates = result$gates,
    all_gates_pass = result$all_gates_pass,
    runtime_seconds = result$runtime_seconds,
    template_converged = result$template$fitting$converged,
    template_iterations = result$template$fitting$iterations,
    template_eligibility = result$template$eligibility$status,
    template_hash = result$template$structural_hash,
    source_loader = source_loader,
    source_fingerprint = digest::digest(list(
      template_code = readLines(file.path(root, "R", "dkge-template.R")),
      validation_code = readLines(file.path(
        root, "inst", "validation", "functional-alignment-template",
        "template-benefit.R"
      ))
    ), algo = "sha256")
  ),
  file.path(output, "known-warp-protocol.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)
saveRDS(result[c(
  "protocol", "config", "metrics", "gates", "all_gates_pass",
  "runtime_seconds"
)], file.path(output, "known-warp-result.rds"))
if (!isTRUE(result$all_gates_pass)) quit(status = 1L)
