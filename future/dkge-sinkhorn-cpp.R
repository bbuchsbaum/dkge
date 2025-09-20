
# dkge-sinkhorn-cpp.R
# Thin R wrappers over the C++ log-domain Sinkhorn for fast transport.

#' Compute Sinkhorn plan in C++ (log-domain)
#' @export
dkge_sinkhorn_plan <- function(C, mu, nu, epsilon = 0.05, max_iter = 300, tol = 1e-6) {
  sinkhorn_plan_log_cpp(C, mu, nu, epsilon, max_iter, tol)
}

#' Pushforward values through Sinkhorn plan without storing T (fast path)
#' @param C cost matrix (n×m)
#' @param mu, nu probability vectors (sum=1)
#' @param values length-n source values
#' @export
dkge_sinkhorn_push <- function(C, mu, nu, values, epsilon = 0.05, max_iter = 300, tol = 1e-6) {
  sinkhorn_push_cpp(C, mu, nu, values, epsilon, max_iter, tol)
}

#' Transport cluster values to a medoid (fast C++ Sinkhorn)
#' @export
dkge_transport_to_medoid_sinkhorn_cpp <- function(v_list, A_list, centroids, sizes = NULL,
                                                  medoid,
                                                  lambda_emb = 1, lambda_spa = 0.5, sigma_mm = 15,
                                                  epsilon = 0.05, max_iter = 300, tol = 1e-6,
                                                  return_plans = FALSE) {
  S <- length(v_list)
  if (is.null(sizes)) sizes <- lapply(v_list, function(v) rep(1, length(v)))
  Ahat <- lapply(A_list, function(A) A / pmax(sqrt(rowSums(A^2)), 1e-8))
  Q <- nrow(A_list[[medoid]])
  Y <- matrix(NA_real_, S, Q)
  plans <- if (return_plans) vector("list", S) else NULL

  sizes_ref <- sizes[[medoid]]; sizes_ref <- sizes_ref / sum(sizes_ref)

  for (s in seq_len(S)) {
    if (s == medoid) {
      Y[s,] <- v_list[[s]]
      if (return_plans) plans[[s]] <- diag(sizes_ref, nrow = Q)
      next
    }
    As <- Ahat[[s]]; Xs <- centroids[[s]]; szs <- sizes[[s]]; szs <- szs / sum(szs)
    Aref <- Ahat[[medoid]]; Xref <- centroids[[medoid]]

    # cost in R (optionally swap to C++ pairwise)
    embC <- pairwise_sqdist_cpp(As, Aref)
    spaC <- pairwise_sqdist_cpp(Xs / sigma_mm, Xref / sigma_mm)
    C <- lambda_emb * embC + lambda_spa * spaC

    if (return_plans) {
      sp <- dkge_sinkhorn_plan(C, szs, sizes_ref, epsilon, max_iter, tol)
      Tplan <- sp$T
      plans[[s]] <- Tplan
      num <- as.numeric(t(v_list[[s]]) %*% Tplan)
      den <- pmax(colSums(Tplan), 1e-12)
      Y[s,] <- num / den
    } else {
      pf <- dkge_sinkhorn_push(C, szs, sizes_ref, v_list[[s]], epsilon, max_iter, tol)
      Y[s,] <- pf$y
    }
  }
  list(value = apply(Y, 2, stats::median, na.rm = TRUE),
       subj_values = Y,
       plans = plans)
}
