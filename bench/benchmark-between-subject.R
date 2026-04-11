#!/usr/bin/env Rscript

if (!requireNamespace("dkge", quietly = TRUE)) {
  if (!requireNamespace("devtools", quietly = TRUE)) {
    stop("Install either the dkge package or devtools to run this benchmark script.",
         call. = FALSE)
  }
  suppressPackageStartupMessages(devtools::load_all(".", quiet = TRUE))
} else {
  suppressPackageStartupMessages(library(dkge))
}

bench_once <- function(expr, iterations = 3L) {
  times <- numeric(iterations)
  for (i in seq_len(iterations)) {
    gc(FALSE)
    times[i] <- unname(system.time(force(expr))[["elapsed"]])
  }
  c(
    iterations = iterations,
    min_sec = min(times),
    median_sec = stats::median(times),
    max_sec = max(times)
  )
}

make_rrr_fixture <- function(S = 24L, V = 2000L, seed = 1L) {
  set.seed(seed)
  dat <- data.frame(
    subject_id = paste0("s", seq_len(S)),
    group = factor(rep(c("A", "B"), length.out = S)),
    trait = scale(rnorm(S), center = TRUE, scale = FALSE)[, 1],
    age = scale(rnorm(S), center = TRUE, scale = FALSE)[, 1]
  )
  design <- dkge_subject_model(~ group * trait + age, dat)
  X <- design$X
  beta <- matrix(rnorm(ncol(X) * V), ncol(X), V)
  Y <- X %*% beta + matrix(rnorm(S * V, sd = 0.15), S, V)
  rownames(Y) <- dat$subject_id
  target <- dkge_make_target(Y = Y, subject_ids = dat$subject_id)
  list(target = target, design = design)
}

make_transport_input <- function(P_vec = c(10, 12, 9, 11, 8, 10), q = 3L, seed = 2L) {
  set.seed(seed)
  S <- length(P_vec)
  effects <- paste0("e", seq_len(q))
  betas <- lapply(seq_len(S), function(s) {
    B <- matrix(rnorm(q * P_vec[s]), q, P_vec[s])
    rownames(B) <- effects
    B
  })
  designs <- lapply(seq_len(S), function(s) {
    X <- qr.Q(qr(matrix(rnorm(80 * q), 80, q)))
    colnames(X) <- effects
    X
  })
  centroids <- lapply(seq_len(S), function(s) matrix(rnorm(P_vec[s] * 3), P_vec[s], 3))
  list(betas = betas, designs = designs, K = diag(q), S = S, P = P_vec, centroids = centroids)
}

make_transport_fixture <- function(seed = 2L) {
  set.seed(seed)
  data <- make_transport_input(P_vec = c(10, 12, 9, 11, 8, 10), q = 3, seed = seed)
  ids <- paste0("s", seq_len(data$S))
  fit <- dkge_fit(
    dkge_data(data$betas, designs = data$designs, subject_ids = ids),
    K = data$K,
    rank = 2
  )
  transport <- dkge_transport_spec(
    centroids = data$centroids,
    medoid = 1L,
    method = "sinkhorn",
    epsilon = 0.1
  )
  contrast <- matrix(c(1, -1, 0), ncol = 1, dimnames = list(NULL, "c1"))
  list(fit = fit, transport = transport, contrast = contrast, ids = ids)
}

run_between_benchmarks <- function(rrr_iterations = 3L,
                                   perm_iterations = 3L,
                                   transport_iterations = 3L) {
  rrr_fixture <- make_rrr_fixture()
  fit_rrr <- dkge_between_rrr(rrr_fixture$target, rrr_fixture$design, rank = 2)
  transport_fixture <- make_transport_fixture()

  results <- rbind(
    data.frame(
      task = "rrr_s_lt_lt_v",
      t(bench_once(
        dkge_between_rrr(rrr_fixture$target, rrr_fixture$design, rank = 2),
        iterations = rrr_iterations
      )),
      row.names = NULL,
      check.names = FALSE
    ),
    data.frame(
      task = "rrr_permutation_refits",
      t(bench_once(
        dkge_between_permute(fit_rrr, terms = c("group", "trait", "group:trait"), B = 49, seed = 91),
        iterations = perm_iterations
      )),
      row.names = NULL,
      check.names = FALSE
    ),
    data.frame(
      task = "transported_target_build",
      t(bench_once(
        dkge_make_target(
          transport_fixture$fit,
          type = "transported_maps",
          contrast = transport_fixture$contrast,
          transport = transport_fixture$transport,
          crossfit = "analytic"
        ),
        iterations = transport_iterations
      )),
      row.names = NULL,
      check.names = FALSE
    )
  )

  results
}

if (sys.nframe() == 0L) {
  print(run_between_benchmarks())
}
