#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[[1L]]) else normalizePath(".")
old_wd <- setwd(root)
on.exit(setwd(old_wd), add = TRUE)

read_int_env <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(default)
  parsed <- suppressWarnings(as.integer(value))
  if (is.na(parsed) || parsed < 1L) stop(name, " must be a positive integer")
  parsed
}

formal_run <- identical(Sys.getenv("DKGE_BETA_CAL_MODE", unset = "formal"),
                        "formal")
n_rep <- read_int_env("DKGE_BETA_CAL_NREP", if (formal_run) 200L else 5L)
n_perm <- read_int_env("DKGE_BETA_CAL_NPERM", 199L)
if (formal_run && (n_rep != 200L || n_perm != 199L)) {
  stop("Formal settings are frozen at 200 replicates and 199 randomizations")
}

git_output <- function(...) {
  out <- system2("git", c("-C", root, ...), stdout = TRUE, stderr = TRUE)
  status <- attr(out, "status") %||% 0L
  if (!identical(status, 0L)) stop(paste(out, collapse = "\n"))
  trimws(paste(out, collapse = "\n"))
}

`%||%` <- function(x, y) if (is.null(x)) y else x

source_sha <- git_output("rev-parse", "HEAD")
if (!grepl("^[0-9a-f]{40}$", source_sha)) stop("Could not resolve exact Git SHA")
if (formal_run) {
  dirty <- git_output("status", "--porcelain", "--untracked-files=all")
  if (nzchar(dirty)) {
    stop("Formal calibration requires a clean source commit")
  }
}

if (!requireNamespace("digest", quietly = TRUE) ||
    !requireNamespace("pkgload", quietly = TRUE)) {
  stop("Calibration requires digest and pkgload")
}

plan_path <- file.path(root, "data-raw", "dkge-beta-full-pipeline-null-plan.md")
runner_path <- file.path(root, "tools", "calibration",
                         "run-beta-full-pipeline-null.R")
file_sha256 <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}
plan_sha <- file_sha256(plan_path)
runner_sha <- file_sha256(runner_path)

suppressPackageStartupMessages(pkgload::load_all(root, quiet = TRUE))

wilson_interval <- function(successes, total, z = stats::qnorm(0.975)) {
  p <- successes / total
  denom <- 1 + z^2 / total
  center <- (p + z^2 / (2 * total)) / denom
  half <- z * sqrt(p * (1 - p) / total + z^2 / (4 * total^2)) / denom
  c(lower = center - half, upper = center + half)
}

make_null_fit <- function(seed) {
  set.seed(seed)
  subject_count <- 8L
  q <- 3L
  parcel_count <- rep(c(3L, 4L), length.out = subject_count)
  effects <- paste0("effect", seq_len(q))
  subject_ids <- paste0("subject", seq_len(subject_count))
  betas <- lapply(seq_len(subject_count), function(subject) {
    value <- matrix(rnorm(q * parcel_count[[subject]]), q,
                    parcel_count[[subject]])
    rownames(value) <- effects
    colnames(value) <- paste0("parcel", seq_len(ncol(value)))
    value
  })
  names(betas) <- subject_ids
  designs <- lapply(seq_len(subject_count), function(subject) {
    value <- diag(q)
    colnames(value) <- effects
    value
  })
  names(designs) <- subject_ids
  centroids <- lapply(seq_len(subject_count), function(subject) {
    parcel <- seq_len(parcel_count[[subject]])
    cbind(
      x = parcel + subject / 20,
      y = sin(parcel + subject / 10),
      z = cos(parcel - subject / 10)
    )
  })
  K <- diag(c(1, 1, 0))
  dimnames(K) <- list(effects, effects)
  fit <- dkge_fit(betas, designs, K = K, rank = 3)
  list(fit = fit, betas = betas, designs = designs, K = K,
       centroids = centroids)
}

run_replicate <- function(replicate_id) {
  seed <- 8600000L + replicate_id * 1000L
  fixture <- make_null_fit(seed)
  fit <- fixture$fit

  set.seed(seed + 101L)
  inference <- dkge_infer(
    fit,
    contrasts = c(effect1 = 1, effect2 = -1, effect3 = 0),
    method = "loso",
    inference = "signflip",
    correction = "maxT",
    n_perm = n_perm,
    alpha = 0.05,
    align = TRUE,
    transport = list(
      centroids = fixture$centroids,
      medoid = 1L,
      mapper = dkge_mapper_spec(
        "sinkhorn", epsilon = 0.2,
        lambda_emb = 0, lambda_spa = 1,
        warm_start = FALSE
      ),
      betas = fixture$betas,
      provenance = dkge_transport_provenance("geometry_only")
    )
  )
  p_transport <- min(unlist(inference$p_adjusted, use.names = FALSE))
  fold_ranks <- vapply(
    inference$contrasts$metadata$bases,
    ncol,
    integer(1)
  )
  alignment_present <- !is.null(inference$contrasts$metadata$procrustes) &&
    length(inference$contrasts$metadata$aligned_bases) == length(fit$Btil)

  targets <- rbind(class1 = c(1, 0, 0), class2 = c(0, 1, 0))
  colnames(targets) <- rownames(fit$K)
  randomization_cache <- new.env(parent = emptyenv())
  full_pipeline_recompute <- function(
      labels, row_data, target, fit, fold_assignments, method, mode,
      lambda, metric, class_weights, standardize_within_fold) {
    cache_key <- paste(labels, collapse = "|")
    cached <- randomization_cache[[cache_key]]
    if (!is.null(cached)) return(cached)

    randomized_betas <- lapply(seq_along(fixture$betas), function(subject) {
      subject_rows <- which(row_data$subject_idx == subject)
      assigned <- as.character(labels[subject_rows])
      source_position <- match(target$class_labels, assigned)
      if (anyNA(source_position)) stop("Randomized labels lost a target class")
      beta <- fixture$betas[[subject]]
      beta[seq_along(target$class_labels), ] <-
        beta[source_position, , drop = FALSE]
      beta
    })
    names(randomized_betas) <- names(fixture$betas)
    randomized_fit <- dkge_fit(
      randomized_betas, fixture$designs, K = fixture$K, rank = 3
    )
    randomized_folds <- as_dkge_folds(
      fold_assignments, fit_or_data = randomized_fit
    )
    randomized_result <- dkge_classify(
      randomized_fit,
      targets = target,
      method = method,
      mode = mode,
      folds = randomized_folds,
      lambda = lambda,
      metric = metric,
      n_perm = 0L,
      class_weights = class_weights,
      standardize_within_fold = standardize_within_fold
    )
    result <- list(metrics = randomized_result$results[[1L]]$metrics)
    randomization_cache[[cache_key]] <- result
    result
  }
  set.seed(seed + 202L)
  classification <- dkge_classify(
    fit,
    targets = targets,
    method = "lda",
    mode = "cell_cross",
    lambda = 0.001,
    metric = "logloss",
    n_perm = n_perm,
    seed = seed + 203L,
    control = list(randomization_recompute = full_pipeline_recompute)
  )
  class_result <- classification$results[[1L]]
  p_classification <- unname(class_result$p_values[["logloss"]])
  p_family <- min(1, 2 * min(p_transport, p_classification))

  data.frame(
    replicate = replicate_id,
    seed = seed,
    p_transport = p_transport,
    reject_transport = p_transport <= 0.05,
    p_classification = p_classification,
    reject_classification = p_classification <= 0.05,
    p_family_bonferroni = p_family,
    reject_family = p_family <= 0.05,
    kernel_rank = fit$kernel_rank,
    effective_rank = fit$effective_rank,
    minimum_fold_rank = min(fold_ranks),
    alignment_present = alignment_present,
    transport_provenance = inference$metadata$transport_provenance$class,
    transport_exactness = inference$exactness,
    classification_claim_scope = class_result$claim_scope,
    classification_representation_scope = class_result$representation_scope,
    lambda = class_result$lambda,
    stringsAsFactors = FALSE
  )
}

started <- Sys.time()
records <- vector("list", n_rep)
for (replicate_id in seq_len(n_rep)) {
  records[[replicate_id]] <- run_replicate(replicate_id)
  if (replicate_id %% 10L == 0L || replicate_id == n_rep) {
    message(sprintf("completed %d/%d", replicate_id, n_rep))
  }
}
replicates <- do.call(rbind, records)
completed <- Sys.time()

summarize_arm <- function(name, reject, p_value, lower_gate, upper_gate) {
  successes <- sum(reject)
  interval <- wilson_interval(successes, length(reject))
  rate <- mean(reject)
  data.frame(
    arm = name,
    replicates = length(reject),
    rejections_005 = successes,
    rejection_005 = rate,
    wilson_lower = interval[["lower"]],
    wilson_upper = interval[["upper"]],
    mean_p = mean(p_value),
    median_p = stats::median(p_value),
    gate_lower = lower_gate,
    gate_upper = upper_gate,
    rate_gate_pass = rate >= lower_gate && rate <= upper_gate,
    wilson_contains_005 = interval[["lower"]] <= 0.05 &&
      interval[["upper"]] >= 0.05,
    stringsAsFactors = FALSE
  )
}

summary <- rbind(
  summarize_arm("transport_maxT", replicates$reject_transport,
                replicates$p_transport, 0.015, 0.085),
  summarize_arm("classification_cell_cross",
                replicates$reject_classification,
                replicates$p_classification, 0.015, 0.085),
  summarize_arm("bonferroni_family", replicates$reject_family,
                replicates$p_family_bonferroni, 0, 0.075)
)
summary$gate_pass <- with(summary,
  ifelse(
    arm == "bonferroni_family",
    rate_gate_pass & wilson_lower <= 0.05,
    rate_gate_pass & wilson_contains_005
  )
)

contract_pass <- all(
  replicates$kernel_rank == 2L,
  replicates$effective_rank == 2L,
  replicates$minimum_fold_rank == 2L,
  replicates$alignment_present,
  replicates$transport_provenance == "geometry_only",
  replicates$transport_exactness ==
    "randomization_exact_sign_invariant_operator",
  replicates$classification_claim_scope == "prospective_heldout_subject",
  replicates$classification_representation_scope ==
    "basis_cross_fitted_without_heldout_subject",
  replicates$lambda == 0.001,
  is.finite(replicates$p_transport),
  is.finite(replicates$p_classification)
)

out_dir <- if (formal_run) {
  file.path(root, "inst", "extdata")
} else {
  Sys.getenv("DKGE_BETA_CAL_OUT", unset = file.path(tempdir(), "dkge-beta-cal-pilot"))
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
prefix <- file.path(out_dir, "dkge-beta-full-pipeline-null")
replicate_path <- paste0(prefix, "-replicates.csv.gz")
summary_path <- paste0(prefix, "-summary.csv")
metadata_path <- paste0(prefix, "-metadata.csv")

con <- gzfile(replicate_path, open = "wt")
utils::write.csv(replicates, con, row.names = FALSE)
close(con)
utils::write.csv(summary, summary_path, row.names = FALSE)

metadata <- data.frame(
  key = c(
    "status", "started_at_utc", "completed_at_utc", "runtime_seconds",
    "source_git_sha", "plan_sha256", "runner_sha256", "R_version",
    "dkge_version", "platform", "n_replicates", "n_randomizations",
    "base_seed", "lambda_policy", "lambda", "transport_provenance",
    "contract_pass", "all_statistical_gates_pass",
    "replicate_artifact_sha256", "summary_artifact_sha256"
  ),
  value = c(
    if (formal_run) "formal" else "pilot_noncertifying",
    format(started, tz = "UTC", usetz = TRUE),
    format(completed, tz = "UTC", usetz = TRUE),
    sprintf("%.3f", as.numeric(difftime(completed, started, units = "secs"))),
    source_sha, plan_sha, runner_sha, R.version.string,
    as.character(utils::packageVersion("dkge")), R.version$platform,
    n_rep, n_perm, 8600000L, "externally_preselected_fixed", 0.001,
    "geometry_only", contract_pass, all(summary$gate_pass) && contract_pass,
    file_sha256(replicate_path), file_sha256(summary_path)
  ),
  stringsAsFactors = FALSE
)
utils::write.csv(metadata, metadata_path, row.names = FALSE)

print(summary)
message("contract_pass: ", contract_pass)
message("artifacts: ", out_dir)
if (formal_run && (!contract_pass || !all(summary$gate_pass))) quit(status = 2L)
