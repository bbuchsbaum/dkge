#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
worker_mode <- length(args) >= 1L && identical(args[[1]], "--worker")

if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("The benchmark requires the suggested package 'devtools'.")
}
devtools::load_all(quiet = TRUE)

q <- 60L
make_kernel <- function() {
  K <- outer(seq_len(q), seq_len(q), function(i, j) exp(-abs(i - j) / 8))
  dimnames(K) <- list(paste0("e", seq_len(q)), paste0("e", seq_len(q)))
  K + diag(0.05, q)
}

make_group_subjects <- function(seed, P = 10000L, S = 18L) {
  set.seed(seed)
  effects <- paste0("e", seq_len(q))
  effect_profile <- scale(sin(seq(0, 3 * pi, length.out = q)),
                          center = TRUE, scale = FALSE)[, 1]
  spatial <- cos(seq(0, 5 * pi, length.out = P))
  lapply(seq_len(S), function(s) {
    counts <- 3 + rpois(q, lambda = exp(seq(log(1), log(14), length.out = q)))
    observed <- which(runif(q) > 0.15)
    B <- tcrossprod(effect_profile, spatial) +
      matrix(rnorm(q * P), q, P) / sqrt(counts)
    rownames(B) <- effects
    colnames(B) <- paste0("v", seq_len(P))
    X <- diag(q)
    colnames(X) <- effects
    covariance <- diag(1 / counts, q)
    dimnames(covariance) <- list(effects, effects)
    suppressWarnings(dkge_subject(
      B, X, id = paste0("s", s), observed_rows = observed,
      effect_n = stats::setNames(counts, effects),
      effect_precision = stats::setNames(counts, effects),
      effect_noise_cov = covariance,
      residual_variance = rep(1, P), residual_df = 120
    ))
  })
}

profile_workload <- function(workload, seed) {
  profile <- tempfile("dkge-rprofmem-")
  Rprofmem(profile)
  start <- proc.time()[["elapsed"]]
  design_bytes <- 0
  P_work <- NA_integer_

  object <- switch(
    workload,
    trial_chunks = {
      Tn <- 240L
      P_work <- 100000L
      chunk_width <- 5000L
      X <- matrix(0, Tn, q)
      cell <- rep(seq_len(q), each = Tn / q)
      X[cbind(seq_len(Tn), cell)] <- 1
      colnames(X) <- paste0("e", seq_len(q))
      design_bytes <- as.numeric(object.size(X))
      source <- function(i) {
        first <- (i - 1L) * chunk_width + 1L
        if (first > P_work) return(NULL)
        width <- min(chunk_width, P_work - first + 1L)
        set.seed(seed + i)
        beta <- matrix(rnorm(q * width), q, width)
        Y <- X %*% beta + matrix(rnorm(Tn * width, sd = 0.5), Tn, width)
        colnames(Y) <- paste0("v", first + seq_len(width) - 1L)
        Y
      }
      suppressWarnings(dkge_trial_subject_chunks(source, X, id = "imaging"))
    },
    group_legacy = {
      P_work <- 10000L
      subjects <- make_group_subjects(seed, P = P_work)
      design_bytes <- sum(vapply(subjects, function(x) as.numeric(object.size(x$design)),
                                 numeric(1)))
      dkge_fit(dkge_data(subjects), K = make_kernel(), rank = 3,
               w_method = "none")
    },
    group_advanced = {
      P_work <- 10000L
      subjects <- make_group_subjects(seed, P = P_work)
      design_bytes <- sum(vapply(subjects, function(x) as.numeric(object.size(x$design)),
                                 numeric(1)))
      suppressWarnings(dkge_fit(
        dkge_data(subjects), K = make_kernel(), rank = 3,
        w_method = "none", effect_weights = dkge_effect_weights("precision"),
        debias = "analytic", missingness = "rescale"
      ))
    },
    fold_refit = {
      P_work <- 10000L
      subjects <- make_group_subjects(seed, P = P_work)
      design_bytes <- sum(vapply(subjects, function(x) as.numeric(object.size(x$design)),
                                 numeric(1)))
      fit <- suppressWarnings(dkge_fit(
        dkge_data(subjects), K = make_kernel(), rank = 3,
        w_method = "none", effect_weights = dkge_effect_weights("precision"),
        debias = "analytic", missingness = "rescale"
      ))
      getFromNamespace(".dkge_fold_weight_context", "dkge")(
        fit, train_ids = seq_len(17)
      )
    },
    poisson_bootstrap = {
      P_work <- 2000L
      subjects <- make_group_subjects(seed, P = P_work)
      design_bytes <- sum(vapply(subjects, function(x) as.numeric(object.size(x$design)),
                                 numeric(1)))
      fit <- suppressWarnings(dkge_fit(
        dkge_data(subjects), K = make_kernel(), rank = 3,
        w_method = "none", effect_weights = dkge_effect_weights("precision"),
        debias = "analytic", missingness = "rescale"
      ))
      set.seed(seed + 5000L)
      result <- NULL
      for (b in seq_len(20)) {
        multiplicity <- rpois(length(subjects), 1)
        if (sum(multiplicity) < 2L) multiplicity[1:2] <- 1L
        result <- getFromNamespace(".dkge_repool_fit", "dkge")(
          fit, sample_weights = multiplicity
        )
      }
      result
    },
    stop("Unknown workload: ", workload)
  )
  elapsed <- proc.time()[["elapsed"]] - start
  Rprofmem(NULL)
  allocation_lines <- readLines(profile, warn = FALSE)
  unlink(profile)
  allocation <- suppressWarnings(as.numeric(sub(" .*", "", allocation_lines)))
  allocation <- allocation[is.finite(allocation)]
  total_allocation <- sum(allocation)
  largest_allocation <- if (length(allocation)) max(allocation) else 0
  object_bytes <- as.numeric(object.size(object))
  forbidden_bytes <- as.numeric(P_work)^2 * 8
  p_by_p_detected <- largest_allocation >= forbidden_bytes
  fields <- c(
    workload = workload,
    seed = seed,
    elapsed_seconds = elapsed,
    total_allocation_bytes = total_allocation,
    largest_allocation_bytes = largest_allocation,
    object_bytes = object_bytes,
    design_bytes = design_bytes,
    P = P_work,
    q = q,
    p_by_p_detected = p_by_p_detected
  )
  message("BENCH_RESULT ", paste(names(fields), fields, sep = "=", collapse = ","))
  invisible(object)
}

if (worker_mode) {
  if (length(args) != 3L) stop("Worker usage: --worker WORKLOAD SEED")
  profile_workload(args[[2]], as.integer(args[[3]]))
  quit(save = "no", status = 0)
}

workloads <- c("trial_chunks", "group_legacy", "group_advanced",
               "fold_refit", "poisson_bootstrap")
seeds <- 19421:19423
script_path <- normalizePath("dev/benchmark-dkge-q60.R", mustWork = TRUE)
rscript <- file.path(R.home("bin"), "Rscript")

parse_fields <- function(line) {
  pieces <- strsplit(sub("^BENCH_RESULT ", "", line), ",", fixed = TRUE)[[1]]
  values <- strsplit(pieces, "=", fixed = TRUE)
  stats::setNames(vapply(values, `[[`, character(1), 2L),
                  vapply(values, `[[`, character(1), 1L))
}

rows <- list()
cursor <- 0L
for (workload in workloads) {
  for (seed in seeds) {
    message("benchmark ", workload, " seed ", seed)
    output <- system2(
      "/usr/bin/time",
      c("-l", rscript, script_path, "--worker", workload, as.character(seed)),
      stdout = TRUE, stderr = TRUE
    )
    result_line <- grep("^BENCH_RESULT ", output, value = TRUE)
    if (length(result_line) != 1L) {
      stop("Benchmark worker failed for ", workload, ":\n",
           paste(output, collapse = "\n"))
    }
    maximum_rss <- grep("maximum resident set size", output, value = TRUE)
    rss <- if (length(maximum_rss)) {
      as.numeric(sub("^[[:space:]]*([0-9]+).*$", "\\1", maximum_rss[[1]]))
    } else {
      NA_real_
    }
    fields <- parse_fields(result_line)
    cursor <- cursor + 1L
    rows[[cursor]] <- data.frame(
      row_type = "run", workload = workload, seed = seed,
      elapsed_seconds = as.numeric(fields[["elapsed_seconds"]]),
      peak_rss_bytes = rss,
      total_allocation_bytes = as.numeric(fields[["total_allocation_bytes"]]),
      largest_allocation_bytes = as.numeric(fields[["largest_allocation_bytes"]]),
      object_bytes = as.numeric(fields[["object_bytes"]]),
      design_bytes = as.numeric(fields[["design_bytes"]]),
      P = as.integer(fields[["P"]]), q = as.integer(fields[["q"]]),
      p_by_p_detected = as.logical(fields[["p_by_p_detected"]]),
      elapsed_iqr = NA_real_, passed = NA,
      threshold = NA_character_, stringsAsFactors = FALSE
    )
  }
}
runs <- do.call(rbind, rows)

limits <- list(
  trial_chunks = c(elapsed = 20, rss = 1.25e9),
  group_legacy = c(elapsed = 20, rss = 2e9),
  group_advanced = c(elapsed = 30, rss = 2e9),
  fold_refit = c(elapsed = 30, rss = 2e9),
  poisson_bootstrap = c(elapsed = 60, rss = 1.5e9)
)
summaries <- lapply(workloads, function(workload) {
  x <- runs[runs$workload == workload, ]
  elapsed <- stats::median(x$elapsed_seconds)
  rss <- stats::median(x$peak_rss_bytes, na.rm = TRUE)
  elapsed_iqr <- stats::IQR(x$elapsed_seconds)
  extra <- if (workload == "trial_chunks") {
    stats::median(x$object_bytes) <= 100e6 &&
      stats::median(x$design_bytes) <= 0.25e6 &&
      stats::median(x$largest_allocation_bytes) <= 100e6
  } else {
    TRUE
  }
  passed <- elapsed <= limits[[workload]][["elapsed"]] &&
    rss <= limits[[workload]][["rss"]] && elapsed_iqr <= elapsed &&
    !any(x$p_by_p_detected) && extra
  data.frame(
    row_type = "summary", workload = workload, seed = NA_integer_,
    elapsed_seconds = elapsed, peak_rss_bytes = rss,
    total_allocation_bytes = stats::median(x$total_allocation_bytes),
    largest_allocation_bytes = stats::median(x$largest_allocation_bytes),
    object_bytes = stats::median(x$object_bytes),
    design_bytes = stats::median(x$design_bytes),
    P = x$P[[1]], q = x$q[[1]], p_by_p_detected = any(x$p_by_p_detected),
    elapsed_iqr = elapsed_iqr, passed = passed,
    threshold = sprintf("elapsed<=%gs; peak_rss<=%g bytes; IQR<=median",
                        limits[[workload]][["elapsed"]],
                        limits[[workload]][["rss"]]),
    stringsAsFactors = FALSE
  )
})
result <- rbind(runs, do.call(rbind, summaries))
dir.create("inst/extdata", recursive = TRUE, showWarnings = FALSE)
utils::write.csv(result, "inst/extdata/dkge-q60-performance.csv",
                 row.names = FALSE)
message(sprintf("completed %d runs; %d/%d workload gates passed",
                nrow(runs), sum(vapply(summaries, `[[`, logical(1), "passed")),
                length(summaries)))
