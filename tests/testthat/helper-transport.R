# Helper utilities for transport-related tests

create_mismatched_data <- function(P_vec = c(3, 5, 4, 6, 4), q = 3, seed = 999) {
  set.seed(seed)
  S <- length(P_vec)
  betas <- lapply(seq_len(S), function(s) {
    matrix(rnorm(q * P_vec[s]), q, P_vec[s])
  })
  designs <- lapply(seq_len(S), function(s) {
    qr.Q(qr(matrix(rnorm(80 * q), 80, q)))
  })
  centroids <- lapply(seq_len(S), function(s) {
    matrix(rnorm(P_vec[s] * 3), P_vec[s], 3)
  })
  list(betas = betas, designs = designs, K = diag(q), S = S, P = P_vec,
       q = q, centroids = centroids)
}

make_simple_transforms <- function(P_vec, target_dim = 2L) {
  lapply(P_vec, function(Ps) {
    mat <- matrix(0, nrow = target_dim, ncol = Ps)
    mat[1, ] <- 1 / Ps
    mat[2, ] <- seq_len(Ps) / sum(seq_len(Ps))
    mat
  })
}

skip_if_no_T4transport <- function() {
  if (!requireNamespace("T4transport", quietly = TRUE)) {
    skip("T4transport package not available")
  }
}
