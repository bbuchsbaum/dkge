
# dkge-mixed.R (v0.5)
# Optional mixed-effects analysis over transported subject maps.

#' Cluster-wise linear mixed-effects using lme4 (optional)
#'
#' Fits per-cluster LME: y ~ Z_fixed + (1|subject), using 'lme4' (and 'lmerTest' if installed)
#' for Satterthwaite df. Returns fixed-effect estimates, SEs, t/z, p-values (per cluster),
#' and an optional multiplicity adjustment via BH FDR.
#'
#' @param Y S×Q matrix of subject maps (transported to medoid)
#' @param Z_fixed S×p fixed-effects design (columns include intercept if desired)
#' @param coef index or name of the fixed effect to test
#' @param ids optional subject IDs (factor length S); if NULL, uses 1..S
#' @param p_adjust NULL | "BH"
#' @param parallel logical; if TRUE and 'future.apply' is available, runs in parallel
#' @return data.frame with per-cluster stats
#' @export
dkge_lmer_map <- function(Y, Z_fixed, coef = 2, ids = NULL, p_adjust = "BH",
                          parallel = FALSE) {
  if (!requireNamespace("lme4", quietly = TRUE))
    stop("Install 'lme4' to use dkge_lmer_map().")
  use_lmerTest <- requireNamespace("lmerTest", quietly = TRUE)

  Y <- as.matrix(Y); S <- nrow(Y); Q <- ncol(Y)
  if (is.null(ids)) ids <- factor(seq_len(S))
  dat <- data.frame(ids = ids, Z_fixed)
  colnames(dat) <- c("ids", paste0("x", seq_len(ncol(Z_fixed))))

  # worker for one column (cluster)
  worker <- function(j) {
    dat$y <- Y[, j]
    fml <- as.formula(paste("y ~", paste(colnames(dat)[-c(1, ncol(dat))], collapse = " + "), "+ (1|ids)"))
    fit <- lme4::lmer(fml, data = dat, REML = TRUE)
    if (use_lmerTest) {
      s <- suppressWarnings(lmerTest::summary(fit))
      cn <- colnames(s$coefficients)
      row <- s$coefficients[coef, , drop = FALSE]
      est <- row[1, "Estimate"]; se <- row[1, "Std. Error"]
      stat <- row[1, "t value"]; p <- row[1, "Pr(>|t|)"]; df <- row[1, "df"]
    } else {
      sm <- summary(fit)
      est <- sm$coefficients[coef, "Estimate"]; se <- sm$coefficients[coef, "Std. Error"]
      stat <- est / (se + 1e-12); p <- 2 * stats::pnorm(-abs(stat)); df <- NA_real_
    }
    c(estimate = est, se = se, stat = stat, p = p, df = df)
  }

  if (parallel && requireNamespace("future.apply", quietly = TRUE)) {
    res <- future.apply::future_sapply(seq_len(Q), worker)
  } else {
    res <- sapply(seq_len(Q), worker)
  }
  res <- t(res); colnames(res) <- c("estimate","se","stat","p","df")
  res <- as.data.frame(res)
  if (!is.null(p_adjust)) res$p_adj <- p.adjust(res$p, method = p_adjust)
  res$cluster <- seq_len(Q)
  res
}

#' Random-effects meta-analysis across subjects (per cluster)
#'
#' Uses DerSimonian-Laird estimator for tau^2, returns pooled mean and z/p-value.
#' This is a lightweight alternative when mixed-effects fitting is unavailable.
#'
#' @param Y S×Q matrix of subject maps
#' @param se optional S×Q matrix of within-subject SEs; if NULL, uses sample SD across subjects
#' @param p_adjust NULL|"BH"
#' @export
dkge_meta_random <- function(Y, se = NULL, p_adjust = "BH") {
  Y <- as.matrix(Y); S <- nrow(Y); Q <- ncol(Y)
  if (is.null(se)) {
    se <- matrix(stats::sd(Y), nrow = S, ncol = Q, byrow = TRUE) / sqrt(S)
  }
  W <- 1 / (se^2 + 1e-12)  # S×Q
  wsum <- colSums(W)
  mu <- colSums(W * Y) / (wsum + 1e-12)
  # Q statistic and tau^2 (DL)
  Qstat <- colSums(W * (Y - rep(mu, each = S))^2)
  c_k <- wsum - colSums(W^2) / (wsum + 1e-12)
  tau2 <- pmax((Qstat - (S - 1)) / (c_k + 1e-12), 0)
  Wstar <- 1 / (se^2 + rep(tau2, each = S) + 1e-12)
  mu_hat <- colSums(Wstar * Y) / (colSums(Wstar) + 1e-12)
  se_hat <- sqrt(1 / (colSums(Wstar) + 1e-12))
  z <- mu_hat / (se_hat + 1e-12)
  p <- 2 * stats::pnorm(-abs(z))
  out <- data.frame(cluster = seq_len(Q), estimate = mu_hat, se = se_hat, stat = z, p = p)
  if (!is.null(p_adjust)) out$p_adj <- p.adjust(out$p, method = p_adjust)
  out
}
