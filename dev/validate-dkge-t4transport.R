#!/usr/bin/env Rscript

if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("The validation runner requires the suggested package 'devtools'.")
}
if (!requireNamespace("T4transport", quietly = TRUE)) {
  stop("T4transport is unavailable; install it in an isolated library first.")
}

devtools::load_all(quiet = TRUE)

validate_fixture <- function(name, C, mu, nu, epsilon,
                             max_iter = 5000L, tol = 1e-12) {
  ours <- dkge:::.dkge_sinkhorn_plan(
    C, mu, nu, epsilon = epsilon, max_iter = max_iter, tol = tol
  )
  reference <- T4transport::sinkhornD(
    C, p = 1, wx = mu, wy = nu, lambda = epsilon,
    maxiter = max_iter, abstol = tol
  )
  ours_objective <- sum(ours * C)
  ref_objective <- reference$distance
  data.frame(
    fixture = name,
    source_points = nrow(C),
    target_points = ncol(C),
    epsilon = epsilon,
    max_iter = max_iter,
    tolerance = tol,
    plan_max_abs_diff = max(abs(unname(ours) - unname(reference$plan))),
    ours_objective = ours_objective,
    reference_objective = ref_objective,
    objective_abs_diff = abs(ours_objective - ref_objective),
    ours_row_marginal_max_abs_error = max(abs(rowSums(ours) - mu)),
    ours_col_marginal_max_abs_error = max(abs(colSums(ours) - nu)),
    reference_row_marginal_max_abs_error =
      max(abs(rowSums(reference$plan) - mu)),
    reference_col_marginal_max_abs_error =
      max(abs(colSums(reference$plan) - nu)),
    t4transport_version = as.character(utils::packageVersion("T4transport")),
    dkge_version = as.character(utils::packageVersion("dkge")),
    r_version = R.version.string,
    platform = R.version$platform,
    os = paste(Sys.info()[c("sysname", "release", "machine")], collapse = " "),
    stringsAsFactors = FALSE
  )
}

set.seed(123)
x1 <- matrix(stats::rnorm(3 * 2), 3, 2)
y1 <- matrix(stats::rnorm(4 * 2) + 0.1, 4, 2)
c1 <- as.matrix(stats::dist(rbind(x1, y1)))[seq_len(3), 3 + seq_len(4)]

set.seed(456)
x2 <- matrix(stats::runif(5 * 3), 5, 3)
y2 <- matrix(stats::runif(6 * 3) + 0.3, 6, 3)
c2 <- as.matrix(stats::dist(rbind(x2, y2)))[seq_len(5), 5 + seq_len(6)]
mu2 <- stats::runif(5)
mu2 <- mu2 / sum(mu2)
nu2 <- stats::runif(6)
nu2 <- nu2 / sum(nu2)

evidence <- rbind(
  validate_fixture(
    "uniform_gaussian_blobs", c1,
    rep(1 / 3, 3), rep(1 / 4, 4), epsilon = 0.1
  ),
  validate_fixture(
    "nonuniform_three_dimensional", c2,
    mu2, nu2, epsilon = 0.2
  )
)
evidence$passed <- with(
  evidence,
  plan_max_abs_diff <= 1e-6 &
    objective_abs_diff <= 1e-6 &
    ours_row_marginal_max_abs_error <= 1e-6 &
    ours_col_marginal_max_abs_error <= 1e-6
)

dir.create("inst/extdata", recursive = TRUE, showWarnings = FALSE)
utils::write.csv(
  evidence,
  "inst/extdata/dkge-t4transport-validation.csv",
  row.names = FALSE
)
print(evidence, row.names = FALSE)
if (!all(evidence$passed)) {
  stop("At least one T4transport differential fixture failed.")
}

