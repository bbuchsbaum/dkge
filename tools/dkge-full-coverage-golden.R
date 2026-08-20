#!/usr/bin/env Rscript

baseline_ref <- "953d2a64ff356d55f4a6fd3278a7f11db018fea6"
fixture_seed <- 20260820L
weight_seed <- 8172026L

arg_value <- function(prefix, args, default = NULL) {
  hit <- args[startsWith(args, prefix)]
  if (!length(hit)) return(default)
  sub(prefix, "", hit[[1L]], fixed = TRUE)
}

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Run this file with Rscript.", call. = FALSE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1L]]))
repo_root <- normalizePath(file.path(dirname(script_path), ".."))
args <- commandArgs(trailingOnly = TRUE)

golden_fixture <- function() {
  set.seed(fixture_seed)
  S <- 6L
  q <- 5L
  P <- 13L
  Tn <- 24L
  effects <- paste0("e", seq_len(q))
  betas <- replicate(S, {
    B <- matrix(stats::rnorm(q * P), q, P)
    rownames(B) <- effects
    B
  }, simplify = FALSE)
  designs <- replicate(S, {
    X <- matrix(stats::rnorm(Tn * q), Tn, q)
    colnames(X) <- effects
    X
  }, simplify = FALSE)
  A <- matrix(stats::rnorm(q * q), q, q)
  K <- crossprod(A) / q + diag(q)
  dimnames(K) <- list(effects, effects)
  list(
    S = S,
    betas = betas,
    designs = designs,
    K = K,
    contrast = stats::setNames(c(-1, 1, 0.5, -0.5, 0), effects)
  )
}

run_worker <- function(library, output, source_sha, source_digest) {
  .libPaths(c(library, .libPaths()))
  suppressPackageStartupMessages(library(dkge))
  fixture <- golden_fixture()

  fit_one <- function(default) {
    set.seed(weight_seed)
    if (default) {
      dkge_fit(fixture$betas, fixture$designs, K = fixture$K, rank = 3L)
    } else {
      dkge_fit(
        fixture$betas, fixture$designs, K = fixture$K, rank = 3L,
        w_method = "none"
      )
    }
  }
  extract <- function(fit) {
    loso <- lapply(seq_len(fixture$S), function(s) {
      dkge_loso_contrast(fit, s = s, contrasts = fixture$contrast)$v
    })
    list(
      K = unname(fit$K),
      Khalf = unname(fit$Khalf),
      Btil = lapply(fit$Btil, unname),
      weights = unname(fit$weights),
      Chat = unname(fit$Chat),
      loso = lapply(loso, unname)
    )
  }

  dependency_names <- c(
    "dkge", "Matrix", "RSpectra", "irlba", "Rcpp", "RcppArmadillo",
    "future.apply", "ggplot2"
  )
  installed <- dependency_names[vapply(
    dependency_names, requireNamespace, logical(1), quietly = TRUE
  )]
  result <- list(
    source_sha = source_sha,
    source_digest = source_digest,
    fixture_seed = fixture_seed,
    weight_seed = weight_seed,
    default = extract(fit_one(TRUE)),
    none = extract(fit_one(FALSE)),
    R = R.version.string,
    platform = R.version$platform,
    dependencies = stats::setNames(vapply(
      installed,
      function(pkg) as.character(utils::packageVersion(pkg)),
      character(1)
    ), installed)
  )
  saveRDS(result, output, version = 3)
  return(invisible(NULL))
}

if ("--worker" %in% args) {
  run_worker(
    library = arg_value("--library=", args),
    output = arg_value("--output=", args),
    source_sha = arg_value("--source-sha=", args),
    source_digest = arg_value("--source-digest=", args)
  )
  quit(save = "no", status = 0L)
}

run_checked <- function(command, command_args, log) {
  status <- system2(command, command_args, stdout = log, stderr = log)
  if (!identical(status, 0L)) {
    stop(
      sprintf(
        "Command failed (%s %s). See %s.",
        command, paste(command_args, collapse = " "), log
      ),
      call. = FALSE
    )
  }
}

copy_tracked_candidate <- function(destination) {
  tracked <- system2(
    "git", c("-C", repo_root, "ls-files"), stdout = TRUE, stderr = TRUE
  )
  if (!length(tracked)) stop("git ls-files returned no candidate files.")
  for (relative in tracked) {
    source <- file.path(repo_root, relative)
    if (!file.exists(source)) next
    target <- file.path(destination, relative)
    dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
    if (!file.copy(source, target, overwrite = TRUE, copy.mode = TRUE)) {
      stop("Could not copy candidate file: ", relative, call. = FALSE)
    }
  }
  tracked[file.exists(file.path(repo_root, tracked))]
}

source_digest <- function(root, tracked) {
  hashes <- tools::md5sum(file.path(root, tracked))
  manifest <- tempfile("dkge-source-manifest-", fileext = ".txt")
  on.exit(unlink(manifest), add = TRUE)
  writeLines(paste(tracked, unname(hashes), sep = "\t"), manifest)
  unname(tools::md5sum(manifest))
}

output_dir <- arg_value(
  "--output-dir=", args,
  file.path(repo_root, "inst", "extdata")
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

temp_root <- tempfile("dkge-full-coverage-golden-")
dir.create(temp_root)
on.exit(unlink(temp_root, recursive = TRUE, force = TRUE), add = TRUE)
baseline_source <- file.path(temp_root, "baseline")
candidate_source <- file.path(temp_root, "candidate")
baseline_library <- file.path(temp_root, "lib-baseline")
candidate_library <- file.path(temp_root, "lib-candidate")
dir.create(baseline_source)
dir.create(candidate_source)
dir.create(baseline_library)
dir.create(candidate_library)

archive <- file.path(temp_root, "baseline.tar")
run_checked(
  "git",
  c("-C", repo_root, "archive", "--format=tar", "-o", archive, baseline_ref),
  file.path(temp_root, "archive.log")
)
utils::untar(archive, exdir = baseline_source)
candidate_files <- copy_tracked_candidate(candidate_source)

candidate_sha <- system2(
  "git", c("-C", repo_root, "rev-parse", "HEAD"),
  stdout = TRUE, stderr = TRUE
)[[1L]]
baseline_tree <- system2(
  "git", c("-C", repo_root, "rev-parse", paste0(baseline_ref, "^{tree}")),
  stdout = TRUE, stderr = TRUE
)[[1L]]
candidate_digest <- source_digest(candidate_source, candidate_files)

r_binary <- file.path(R.home("bin"), "R")
rscript_binary <- file.path(R.home("bin"), "Rscript")
install <- function(source, library, label) {
  run_checked(
    r_binary,
    c(
      "CMD", "INSTALL", "--no-multiarch", "--with-keep.source",
      "-l", library, source
    ),
    file.path(temp_root, paste0("install-", label, ".log"))
  )
}
install(baseline_source, baseline_library, "baseline")
install(candidate_source, candidate_library, "candidate")

baseline_rds <- file.path(temp_root, "baseline.rds")
candidate_rds <- file.path(temp_root, "candidate.rds")
worker <- function(library, output, sha, digest, label) {
  run_checked(
    rscript_binary,
    c(
      script_path, "--worker", paste0("--library=", library),
      paste0("--output=", output), paste0("--source-sha=", sha),
      paste0("--source-digest=", digest)
    ),
    file.path(temp_root, paste0("worker-", label, ".log"))
  )
}
worker(baseline_library, baseline_rds, baseline_ref, baseline_tree, "baseline")
worker(candidate_library, candidate_rds, candidate_sha, candidate_digest, "candidate")

baseline <- readRDS(baseline_rds)
candidate <- readRDS(candidate_rds)
max_abs <- function(reference, observed) {
  max(abs(unlist(reference, use.names = FALSE) - unlist(observed, use.names = FALSE)))
}
relative_error <- function(reference, observed) {
  reference <- unlist(reference, use.names = FALSE)
  observed <- unlist(observed, use.names = FALSE)
  sqrt(sum((reference - observed)^2)) / sqrt(sum(reference^2))
}

metric_specs <- data.frame(
  artifact = c("K", "Khalf", "Btil", "weights", "Chat", "loso"),
  metric = c("max_abs", "max_abs", "max_abs", "max_abs", "relative", "max_abs"),
  threshold = c(0, 0, 0, 0, 1e-15, 1e-14),
  stringsAsFactors = FALSE
)
metric_rows <- lapply(c("default", "none"), function(path) {
  do.call(rbind, lapply(seq_len(nrow(metric_specs)), function(i) {
    spec <- metric_specs[i, ]
    value <- if (identical(spec$metric, "relative")) {
      relative_error(
        baseline[[path]][[spec$artifact]], candidate[[path]][[spec$artifact]]
      )
    } else {
      max_abs(baseline[[path]][[spec$artifact]], candidate[[path]][[spec$artifact]])
    }
    data.frame(
      path = path,
      artifact = spec$artifact,
      metric = spec$metric,
      value = value,
      threshold = spec$threshold,
      pass = is.finite(value) && value <= spec$threshold,
      stringsAsFactors = FALSE
    )
  }))
})
metrics <- do.call(rbind, metric_rows)
row.names(metrics) <- NULL

format_dependencies <- function(x) {
  paste(paste(names(x), unname(x), sep = "="), collapse = ";")
}
command_record <- paste(
  "Rscript tools/dkge-full-coverage-golden.R",
  paste0("--output-dir=", output_dir)
)
metadata <- data.frame(
  role = c("baseline", "candidate"),
  source_sha = c(baseline$source_sha, candidate$source_sha),
  source_digest = c(baseline$source_digest, candidate$source_digest),
  R = c(baseline$R, candidate$R),
  platform = c(baseline$platform, candidate$platform),
  dependencies = c(
    format_dependencies(baseline$dependencies),
    format_dependencies(candidate$dependencies)
  ),
  fixture_seed = fixture_seed,
  weight_seed = weight_seed,
  command = command_record,
  stringsAsFactors = FALSE
)

metrics_path <- file.path(output_dir, "dkge-full-coverage-golden.csv")
metadata_path <- file.path(output_dir, "dkge-full-coverage-golden-metadata.csv")
utils::write.csv(metrics, metrics_path, row.names = FALSE, quote = TRUE)
utils::write.csv(metadata, metadata_path, row.names = FALSE, quote = TRUE)

print(metrics, row.names = FALSE, digits = 17)
cat("Metadata:", metadata_path, "\n")
cat("Metrics:", metrics_path, "\n")
if (!all(metrics$pass)) quit(save = "no", status = 1L)
