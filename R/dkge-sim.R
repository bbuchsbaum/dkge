#' Simulate toy DKGE datasets with known factorial structure
#'
#' Generates subject-level beta matrices \eqn{B_s} (`q x P_s`) in effect space with
#' user-selected active design terms, controlled signal-to-noise ratio (SNR), and
#' per-subject noise levels. The routine also returns the per-subject design
#' 'rulers' (typically identity matrices) and the design kernel \eqn{K} so the
#' synthetic dataset can be fed directly to [dkge_fit()]. Components are planted
#' inside chosen effect blocks, then \eqn{K}-orthonormalised to provide a ground
#' truth subspace.
#'
#' Construction follows the model
#' \deqn{B_s = U_{true} M_s^\top + E_s,}
#' where `U_true` lies inside the requested design-term blocks, `M_s` encodes the
#' subject-specific spatial loadings, and `E_s` is Gaussian noise adjusted to hit
#' the requested Frobenius SNR.
#'
#' @param factors Named list describing the experimental factors (as for
#'   [design_kernel()]), e.g. `list(A = list(L = 2), B = list(L = 3))`.
#' @param terms Optional list of character vectors specifying which terms to
#'   include in the kernel. Defaults to the full factorial set used by
#'   [design_kernel()].
#' @param active_terms Character vector of term names to activate (e.g.,
#'   `c("A", "B", "A:B")`). These must be present in the kernel.
#' @param r_per_term Named integer vector giving how many latent columns to draw
#'   within each active term. Defaults to one column per term.
#' @param S Number of subjects.
#' @param P Either a single integer (clusters/voxels per subject) or a length-`S`
#'   integer vector.
#' @param snr Target Frobenius SNR (signal / noise) per subject. Scalar or
#'   length-`S`.
#' @param noise_scales Optional length-`S` multipliers applied to the noise
#'   standard deviation for each subject.
#' @param seed RNG seed for reproducibility.
#' @param contrasts_type Factor-contrast system, one of `"helmert"` or `"sum"`.
#' @return A list with the following entries:
#'   \describe{
#'     \item{B_list}{List of `q x P_s` beta matrices.}
#'     \item{X_list}{List of subject-level design rulers (identity matrices).}
#'     \item{K}{`q x q` design kernel.}
#'     \item{info}{Design metadata returned by [design_kernel()].}
#'     \item{U_true}{Ground-truth `q x r_true` component basis (\eqn{K}-orthonormal).}
#'     \item{M_list}{Per-subject spatial patterns (`P_s x r_true`).}
#'     \item{active_cols}{Named list mapping each active term to the selected
#'       column indices in the effect basis.}
#'     \item{subject_ids}{Character vector of subject identifiers.}
#'   }
#' @export
dkge_sim_toy <- function(factors,
                         terms = NULL,
                         active_terms,
                         r_per_term = NULL,
                         S = 3,
                         P = 10,
                         snr = 8,
                         noise_scales = NULL,
                         seed = 1L,
                         contrasts_type = c("helmert", "sum")) {
  set.seed(seed)
  contrasts_type <- match.arg(contrasts_type)

  Ls <- vapply(factors, function(f) as.integer(f$L), integer(1))
  if (is.null(names(Ls))) names(Ls) <- names(factors)
  contrast_list <- if (identical(contrasts_type, "helmert")) {
    helmert_contrasts(Ls)
  } else {
    sum_contrasts(Ls)
  }

  Kobj <- design_kernel(factors,
                        terms = terms,
                        basis = "effect",
                        contrasts = contrast_list,
                        include_intercept = TRUE,
                        rho0 = 1e-8)
  K <- Kobj$K
  info <- Kobj$info
  q <- nrow(K)

  term_names <- info$term_names
  block_index <- setNames(seq_along(term_names), term_names)
  if (is.null(r_per_term)) {
    r_per_term <- setNames(rep(1L, length(active_terms)), active_terms)
  }

  active_cols <- list()
  cols <- integer(0)
  for (tname in active_terms) {
    if (!tname %in% term_names) {
      stop("Active term ", tname, " not found in kernel terms.")
    }
    bidx <- block_index[[tname]]
    available <- info$blocks[[bidx]]
    pick <- sample(available, size = r_per_term[[tname]], replace = FALSE)
    active_cols[[tname]] <- pick
    cols <- c(cols, pick)
  }
  r_true <- length(cols)

  U0 <- matrix(0, q, r_true)
  for (j in seq_along(cols)) {
    U0[cols[j], j] <- 1
  }
  U_true <- dkge_k_orthonormalize(U0, K)

  if (length(P) == 1L) P <- rep(P, S)
  if (is.null(noise_scales)) noise_scales <- rep(1, S)
  if (length(snr) == 1L) snr <- rep(snr, S)

  B_list <- vector("list", S)
  M_list <- vector("list", S)
  X_list <- vector("list", S)

  for (s in seq_len(S)) {
    P_s <- P[s]
    M_s <- matrix(stats::rnorm(P_s * r_true), P_s, r_true)
    Signal <- U_true %*% t(M_s)
    sigF <- sqrt(sum(Signal^2))
    target_noiseF <- if (snr[s] > 0) sigF / snr[s] else 0
    sigma <- if (q > 0 && P_s > 0) (target_noiseF / sqrt(q * P_s)) * noise_scales[s] else 0
    Noise <- matrix(stats::rnorm(q * P_s, sd = sigma), q, P_s)

    B_list[[s]] <- Signal + Noise
    M_list[[s]] <- M_s
    X_list[[s]] <- diag(q)
  }

  list(B_list = B_list,
       X_list = X_list,
       K = K,
       info = info,
       U_true = U_true,
       M_list = M_list,
       active_cols = active_cols,
       subject_ids = paste0("sub", seq_len(S)))
}

#' Principal-angle cosines in the K-metric
#'
#' @param U Matrix with `q` rows.
#' @param V Matrix with `q` rows.
#' @param K Symmetric positive (semi-)definite matrix defining the metric.
#' @return Numeric vector of singular values in [0, 1].
#' @export
dkge_cosines_K <- function(U, V, K) {
  s <- svd(t(U) %*% K %*% V, nu = 0, nv = 0)$d
  s[s > 1] <- 1
  s[s < 0] <- 0
  s
}
