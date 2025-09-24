# Helper utilities for LOSO null permutation tests with adaptive voxel weights.

null_betas <- function(nsub, V, Q = 2, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  lapply(seq_len(nsub), function(s) matrix(rnorm(Q * V), nrow = Q, ncol = V))
}

dkge_test_reliability_prior <- function(Q = 2, V = 120, nsub = 12, seed = 11) {
  set.seed(seed)
  runs1 <- lapply(seq_len(nsub), function(s) matrix(rnorm(Q * V), nrow = Q, ncol = V))
  runs2 <- lapply(seq_len(nsub), function(s) matrix(rnorm(Q * V), nrow = Q, ncol = V))
  r <- sapply(seq_len(V), function(v) {
    x <- vapply(seq_len(nsub), function(s) mean(runs1[[s]][, v]), numeric(1))
    y <- vapply(seq_len(nsub), function(s) mean(runs2[[s]][, v]), numeric(1))
    suppressWarnings(stats::cor(x, y))
  })
  r[!is.finite(r)] <- 0
  (abs(r))^2 + 1e-6
}

# Internal worker that returns the LOSO mean statistic p-value for a single
# synthetic experiment under the null. Optionally accepts a prior vector and a
# mixing weight (log-space convex combination when combine = "product").
dkge_test_loso_delta_pvalue <- function(B_list,
                                        adapt = c("kenergy_prec", "kenergy", "none"),
                                        w_prior = NULL,
                                        mix = 1,
                                        n_perm = 199,
                                        lambda = 1e-3) {
  adapt <- match.arg(adapt)
  stopifnot(length(B_list) >= 6)
  Q <- nrow(B_list[[1L]])
  V <- ncol(B_list[[1L]])
  if (Q < 2) stop("Require Q >= 2 for the two-row delta statistic")

  S <- length(B_list)
  Iq <- diag(Q)

  if (is.null(w_prior)) w_prior <- rep(1, V)
  stopifnot(length(w_prior) == V)

  adapt_fun  <- dkge:::.dkge_adapt_weights
  combine_fn <- dkge:::.dkge_combine_weights

  deltas <- vector("list", S)
  precis <- vector("list", S)

  for (s in seq_len(S)) {
    train_ids <- setdiff(seq_len(S), s)
    B_train <- B_list[train_ids]

    w_adapt <- switch(adapt,
      none = rep(1, V),
      kenergy = adapt_fun(B_list = B_train, adapt = "kenergy", K = Iq,
                          winsor = 0.9999),
      kenergy_prec = {
        ke <- adapt_fun(B_list = B_train, adapt = "kenergy", K = Iq,
                        winsor = 0.9999)
        pr <- Reduce("+", lapply(B_train, function(B) 1 / (colMeans(B^2) + 1e-8))) / length(B_train)
        ke * (pr / mean(pr))
      }
    )

    w_total <- combine_fn(w_prior, w_adapt,
                          combine = "product", mix = mix,
                          shrink  = list(alpha = 1, winsor = 0.9999, normalize = "mean", roi_smooth = FALSE))

    D_train <- vapply(train_ids, function(i) (B_list[[i]][1, ] - B_list[[i]][2, ]) * sqrt(w_total), numeric(V))
    D_train <- t(D_train)
    var_train <- colMeans(D_train^2)
    prec_train <- 1 / (var_train + lambda)

    deltas[[s]] <- (B_list[[s]][1, ] - B_list[[s]][2, ]) * sqrt(w_total)
    precis[[s]] <- prec_train
  }

  g_obs <- numeric(S)
  for (s in seq_len(S)) {
    train_ids <- setdiff(seq_len(S), s)
    D_train <- t(vapply(train_ids, function(i) deltas[[i]], numeric(V)))
    mu_train <- colMeans(D_train)
    d_obs <- precis[[s]] * mu_train
    g_obs[s] <- sum(deltas[[s]] * d_obs)
  }
  T_obs <- mean(g_obs)

  T_perm <- numeric(n_perm)
  for (p in seq_len(n_perm)) {
    sgn <- sample(c(-1, 1), size = S, replace = TRUE)
    g_perm <- numeric(S)
    for (s in seq_len(S)) {
      train_ids <- setdiff(seq_len(S), s)
      D_train <- t(vapply(train_ids, function(i) deltas[[i]] * sgn[i], numeric(V)))
      mu_p <- colMeans(D_train)
      d_p  <- precis[[s]] * mu_p
      g_perm[s] <- sum((deltas[[s]] * sgn[s]) * d_p)
    }
    T_perm[p] <- mean(g_perm)
  }

  (1 + sum(abs(T_perm) >= abs(T_obs))) / (n_perm + 1)
}

# Convenience: run many replicates and return the collection of permutation
# p-values for uniformity checks.
dkge_test_null_uniformity <- function(nrep = 30,
                                      nsub = 12,
                                      V = 120,
                                      adapt = c("kenergy_prec", "kenergy", "none"),
                                      w_prior = NULL,
                                      mix = 1,
                                      n_perm = 199,
                                      seed = 123) {
  adapt <- match.arg(adapt)
  rng <- seed
  pvals <- numeric(nrep)
  for (r in seq_len(nrep)) {
    rng <- rng + 1L
    B_list <- null_betas(nsub = nsub, V = V, seed = rng)
    pvals[r] <- dkge_test_loso_delta_pvalue(B_list,
                                            adapt = adapt,
                                            w_prior = w_prior,
                                            mix = mix,
                                            n_perm = n_perm)
  }
  pvals
}
