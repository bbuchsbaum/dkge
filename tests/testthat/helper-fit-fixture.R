# Shared fixtures for dkge tests

make_small_fit <- function(S = 3, q = 3, P = 4, T = 12, rank = 2, seed = 42) {
  set.seed(seed)
  betas <- replicate(S, matrix(rnorm(q * P), q, P), simplify = FALSE)
  designs <- replicate(S, {
    X <- matrix(rnorm(T * q), T, q)
    qr.Q(qr(X))
  }, simplify = FALSE)
  fit <- dkge_fit(betas, designs, K = diag(q), rank = rank)
  list(fit = fit, betas = betas, q = q)
}
