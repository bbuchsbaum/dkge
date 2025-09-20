
# dkge-pipeline.R (v0.4)
# High-level orchestration: fit -> LOSO -> transport -> inference -> write

#' End-to-end DKGE pipeline (fit → LOSO contrasts → transport → signflip inference)
#'
#' @param B_list list of q×P_s betas
#' @param X_list list of T_s×q designs
#' @param K q×q kernel
#' @param contrasts list of q-vectors (named)
#' @param rank,ridge,w_method,w_tau DKGE fit options
#' @param transport list with fields: centroids, sizes (optional), medoid (index), method ("sinkhorn_cpp"|"sinkhorn_R"), and cost params
#' @param inference list with fields: B (permutations), tail, center
#' @return list with per-contrast group maps, p-values, and diagnostics
#' @export
dkge_run <- function(B_list, X_list, K, contrasts,
                     rank = 6, ridge = 1e-6, w_method = "mfa_sigma1", w_tau = 0.3,
                     transport = list(method = "sinkhorn_cpp", lambda_emb = 1, lambda_spa = 0.5, sigma_mm = 15,
                                      epsilon = 0.05, max_iter = 300, tol = 1e-6, medoid = 1L,
                                      centroids = NULL, sizes = NULL, return_plans = FALSE),
                     inference = list(B = 2000, tail = "two.sided", center = "mean")) {
  # 1) Fit DKGE
  fit <- dkge_fit(B_list, X_list, K, w_method = w_method, w_tau = w_tau, ridge = ridge, rank = rank)

  # 2) For each contrast: LOSO values then transport
  A_list <- lapply(seq_along(B_list), function(s) t(fit$Btil[[s]]) %*% fit$K %*% fit$U)
  # centroids & sizes required
  stopifnot(!is.null(transport$centroids))
  med <- transport$medoid %||% 1L

  res <- list()
  for (nm in names(contrasts)) {
    c <- contrasts[[nm]]
    loso <- lapply(seq_along(B_list), function(s) dkge_loso_contrast(fit, s, c, ridge = ridge))
    v_list <- lapply(loso, `[[`, "v")

    if (transport$method == "sinkhorn_cpp") {
      ot <- dkge_transport_to_medoid_sinkhorn_cpp(
        v_list, A_list, transport$centroids, sizes = transport$sizes,
        medoid = med,
        lambda_emb = transport$lambda_emb, lambda_spa = transport$lambda_spa, sigma_mm = transport$sigma_mm,
        epsilon = transport$epsilon, max_iter = transport$max_iter, tol = transport$tol,
        return_plans = transport$return_plans %||% FALSE
      )
    } else {
      ot <- dkge_transport_to_medoid_sinkhorn(
        v_list, A_list, transport$centroids, sizes = transport$sizes,
        medoid = med,
        lambda_emb = transport$lambda_emb, lambda_spa = transport$lambda_spa, sigma_mm = transport$sigma_mm,
        epsilon = transport$epsilon, max_iter = transport$max_iter, tol = transport$tol
      )
    }

    # 3) Inference: sign-flip maxT on transported subject maps
    inf <- dkge_signflip_maxT(ot$subj_values, B = inference$B %||% 2000,
                              tail = inference$tail %||% "two.sided",
                              center = inference$center %||% "mean")

    res[[nm]] <- list(group = ot$value, subj = ot$subj_values, p = inf$p, stat = inf$stat,
                      plans = ot$plans %||% NULL, flips = inf$flips, maxnull = inf$maxnull)
  }
  res
}
