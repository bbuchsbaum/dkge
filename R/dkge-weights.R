# dkge-weights.R
# Voxel-level weighting specification for DKGE fits and cross-fitting.

#' Create a DKGE voxel-weight specification
#'
#' Constructs a lightweight object describing how voxel-level weights are
#' resolved and applied when building weighted second moments inside DKGE. The
#' specification is attached to fits and reused by fold builders to ensure
#' cross-fitting remains leak-free.
#'
#' @param prior Optional prior weights: numeric vector of length V, logical mask,
#'   integer indices, or ROI labels (factor/integer) per voxel. Helpers
#'   [dkge_weights_prior_mask()] and [dkge_weights_prior_roi()] ease construction.
#' @param adapt Adaptive weighting rule, one of `"none"`, `"kenergy"`,
#'   `"precision"`, `"kenergy_prec"`, `"reliability"` (reserved for future use).
#' @param combine How prior and adaptive sources combine: `"product"` (default),
#'   `"sum"`, `"override_adapt"`, or `"override_prior"`.
#' @param mix Numeric in [0,1] controlling the relative influence of the adaptive
#'   component. Interpreted in log-space for `combine = "product"`.
#' @param shrink List with fields `alpha` (shrink towards uniform), `winsor`
#'   (upper quantile cap), `normalize` (`"mean"` or `"sum"`), and optional
#'   `roi_smooth = TRUE` to median-smooth within ROIs.
#' @param scope Either `"fold"` (default: compute adapt weights from training
#'   subjects within each fold) or `"subject"` (per-subject adaptive weights
#'   averaged for fold pooling).
#' @param k_weight Optional effect-space kernel for k-energy rules. When `NULL`
#'   we reuse the kernel stored in the fit, with optional factor `collapse`.
#' @param collapse Optional list describing factor collapses (e.g.,
#'   `list(time = "mean")` or `list(time = list(method = "mean", window = 3:8))`).
#' @param roi Optional ROI labels used when `shrink$roi_smooth = TRUE`.
#' @return Object of class `"dkge_weights"`.
#' @export
dkge_weights <- function(prior = NULL,
                         adapt = c("none", "kenergy", "precision", "kenergy_prec", "reliability"),
                         combine = c("product", "sum", "override_adapt", "override_prior"),
                         mix = 0.6,
                         shrink = list(alpha = 0.5, winsor = 0.99, normalize = "mean", roi_smooth = FALSE),
                         scope = c("fold", "subject"),
                         k_weight = NULL,
                         collapse = NULL,
                         roi = NULL) {
  adapt <- match.arg(adapt)
  combine <- match.arg(combine)
  scope <- match.arg(scope)

  shrink <- .dkge_weights_shrink_defaults(shrink)

  structure(list(
    prior = prior,
    adapt = adapt,
    combine = combine,
    mix = .dkge_clamp01(mix),
    shrink = shrink,
    scope = scope,
    k_weight = k_weight,
    collapse = collapse,
    roi = roi
  ), class = "dkge_weights")
}

#' @export
print.dkge_weights <- function(x, ...) {
  cat("dkge weight specification\n")
  cat("  prior   :", .dkge_weights_prior_summary(x$prior), "\n")
  cat("  adapt   :", x$adapt, "(scope =", x$scope, ")\n")
  cat("  combine :", x$combine, "(mix =", sprintf("%.2f", x$mix), ")\n")
  cat("  shrink  : alpha =", sprintf("%.2f", x$shrink$alpha),
      ", winsor =", sprintf("%.3f", x$shrink$winsor),
      ", normalize =", x$shrink$normalize,
      if (isTRUE(x$shrink$roi_smooth)) ", roi_smooth = TRUE" else "",
      "\n", sep = "")
  if (!is.null(x$collapse)) {
    cat("  collapse:", .dkge_collapse_summary(x$collapse), "\n")
  }
  invisible(x)
}

# -------------------------------------------------------------------------
# Prior helpers ------------------------------------------------------------
# -------------------------------------------------------------------------

#' Build prior weights from a mask
#' @param mask Logical or numeric vector identifying voxels in the ROI.
#' @param value_in,value_out Weights inside/outside the ROI (defaults 1/0).
#' @return Numeric vector of prior weights.
#' @export
dkge_weights_prior_mask <- function(mask, value_in = 1, value_out = 0) {
  if (is.null(mask)) stop("`mask` cannot be NULL")
  m <- as.numeric(mask)
  m[is.na(m)] <- 0
  inside <- m != 0
  m[inside] <- value_in
  m[!inside] <- value_out
  m
}

#' Build prior weights from ROI labels
#' @param labels Integer/factor ROI labels per voxel.
#' @param roi_values Optional named numeric multipliers per ROI.
#' @return Numeric vector of prior weights.
#' @export
dkge_weights_prior_roi <- function(labels, roi_values = NULL) {
  if (is.null(labels)) stop("`labels` cannot be NULL")
  lab <- as.integer(as.factor(labels))
  w <- rep(1, length(lab))
  if (!is.null(roi_values)) {
    if (is.null(names(roi_values))) {
      stop("`roi_values` must be a named vector matching ROI labels")
    }
    roi_map <- match(names(roi_values), levels(as.factor(labels)))
    for (i in seq_along(roi_values)) {
      idx <- which(lab == roi_map[i])
      if (length(idx)) w[idx] <- as.numeric(roi_values[[i]])
    }
  }
  w
}

# -------------------------------------------------------------------------
# Internals: shrink/normalize helpers -------------------------------------
# -------------------------------------------------------------------------

.dkge_clamp01 <- function(x) {
  x[is.na(x)] <- 0
  pmax(0, pmin(1, x))
}

.dkge_weights_shrink_defaults <- function(sh) {
  sh$alpha <- if (is.null(sh$alpha)) 0.5 else .dkge_clamp01(sh$alpha)
  sh$winsor <- if (is.null(sh$winsor)) 0.99 else max(0.9, min(0.9999, sh$winsor))
  sh$normalize <- if (is.null(sh$normalize)) "mean" else match.arg(sh$normalize, c("mean", "sum"))
  sh$roi_smooth <- isTRUE(sh$roi_smooth)
  sh
}

.dkge_norm_vec <- function(x, method = c("mean", "sum"), eps = 1e-12) {
  method <- match.arg(method)
  finite <- is.finite(x)
  if (!any(finite)) return(rep(1, length(x)))
  if (method == "mean") {
    m <- mean(x[finite], na.rm = TRUE)
    if (!is.finite(m) || m < eps) m <- 1
    x / m
  } else {
    s <- sum(x[finite], na.rm = TRUE)
    if (!is.finite(s) || s < eps) s <- 1
    x / s
  }
}

.dkge_winsor <- function(x, upper = 0.99) {
  finite <- x[is.finite(x)]
  if (!length(finite)) return(x[])  # return empty vector if needed
  thresh <- stats::quantile(finite, probs = upper, names = FALSE, na.rm = TRUE)
  x[x > thresh] <- thresh
  x
}

.dkge_weights_prior_summary <- function(prior) {
  if (is.null(prior)) return("uniform")
  cls <- class(prior)[1]
  paste0(cls, " length=", length(prior))
}

.dkge_collapse_summary <- function(collapse) {
  paste(vapply(names(collapse), function(nm) {
    rule <- collapse[[nm]]
    if (is.character(rule)) {
      sprintf("%s=%s", nm, paste(rule, collapse = ","))
    } else if (is.list(rule) && !is.null(rule$window)) {
      rng <- range(rule$window)
      sprintf("%s=%s[%d:%d]", nm, rule$method %||% "mean", rng[1], rng[2])
    } else if (is.numeric(rule)) {
      sprintf("%s=weights[%d]", nm, length(rule))
    } else {
      nm
    }
  }, character(1)), collapse = "; ")
}

.dkge_symmetrize <- function(M) {
  0.5 * (M + t(M))
}

# -------------------------------------------------------------------------
# Collapse helpers --------------------------------------------------------
# -------------------------------------------------------------------------

.dkge_collapse_weights <- function(L, rule) {
  if (is.character(rule)) {
    if (length(rule) != 1L || rule != "mean") stop("Unsupported collapse rule" )
    rep(1 / L, L)
  } else if (is.list(rule)) {
    method <- tolower(rule$method %||% "mean")
    if (method != "mean") stop("Only method = 'mean' implemented for collapse")
    window <- rule$window %||% seq_len(L)
    window <- intersect(window, seq_len(L))
    if (!length(window)) stop("Collapse window empty")
    w <- rep(0, L)
    w[window] <- 1 / length(window)
    w
  } else if (is.numeric(rule)) {
    if (length(rule) != L) stop("Numeric collapse weights must match length L")
    if (sum(rule) <= 0) stop("Collapse weights must sum to positive")
    as.numeric(rule) / sum(rule)
  } else {
    stop("Unsupported collapse specification")
  }
}

.dkge_weights_resolve_k <- function(weights, kernel_info) {
  if (!inherits(weights, "dkge_weights")) stop("`weights` must be dkge_weights()")
  K <- weights$k_weight
  if (!is.null(K)) return(.dkge_symmetrize(as.matrix(K)))

  K <- kernel_info$K %||% NULL
  if (is.null(K)) stop("Kernel information required to resolve k-energy weights")
  K <- as.matrix(K)

  collapse <- weights$collapse
  if (is.null(collapse)) return(.dkge_symmetrize(K))

  map <- kernel_info$map %||% kernel_info$info$map %||% NULL
  levels <- NULL
  if (!is.null(kernel_info$info)) levels <- kernel_info$info$levels %||% NULL
  if (is.null(map) || is.null(levels)) {
    warning("Collapse requested but kernel info lacks mapping; using base kernel")
    return(.dkge_symmetrize(K))
  }

  kron_list <- vector("list", length(levels))
  names(kron_list) <- names(levels)
  for (nm in names(levels)) {
    L <- levels[[nm]]
    rule <- collapse[[nm]]
    if (is.null(rule)) {
      kron_list[[nm]] <- diag(L)
    } else {
      w <- .dkge_collapse_weights(L, rule)
      kron_list[[nm]] <- tcrossprod(w)
    }
  }
  K_cell <- Reduce(function(a, b) kronecker(a, b), kron_list)
  Keff <- t(map) %*% K_cell %*% map
  .dkge_symmetrize(Keff)
}

# -------------------------------------------------------------------------
# Adaptive weights --------------------------------------------------------
# -------------------------------------------------------------------------

.dkge_eval_prior <- function(prior, V) {
  if (is.null(prior)) return(rep(1, V))
  if (is.logical(prior)) {
    w <- as.numeric(prior)
  } else if (is.numeric(prior)) {
    w <- as.numeric(prior)
    if (length(w) != V) stop("Numeric prior must have length matching voxels")
  } else if (is.factor(prior)) {
    ids <- as.integer(prior)
    w <- rep(1, length(ids))
    if (length(w) != V) stop("Factor prior length mismatch")
  } else if (is.integer(prior)) {
    w <- rep(0, V)
    w[prior] <- 1
  } else {
    stop("Unsupported prior type")
  }
  if (length(w) != V) stop("Prior length does not match voxel count")
  w[!is.finite(w) | w < 0] <- 0
  if (!any(w > 0)) w[] <- 1
  .dkge_norm_vec(w, "mean")
}

.dkge_adapt_weights <- function(B_list,
                                adapt = "kenergy",
                                K = NULL,
                                sigma2_list = NULL,
                                winsor = 0.99) {
  V <- ncol(B_list[[1]])
  if (adapt == "none") return(rep(1, V))

  if (adapt %in% c("kenergy", "kenergy_prec") && is.null(K)) {
    stop("k-energy weighting requires an effect-space kernel")
  }

  energy_acc <- NULL
  prec_acc <- NULL
  n <- 0L
  for (idx in seq_along(B_list)) {
    B <- B_list[[idx]]
    if (adapt %in% c("kenergy", "kenergy_prec")) {
      KB <- K %*% B
      e <- colSums(B * KB)
      energy_acc <- if (is.null(energy_acc)) e else (energy_acc + e)
    }
    if (adapt %in% c("precision", "kenergy_prec")) {
      if (!is.null(sigma2_list)) {
        sig2 <- sigma2_list[[idx]]
      } else {
        sig2 <- colMeans(B * B)
      }
      p <- 1 / (sig2 + 1e-8)
      prec_acc <- if (is.null(prec_acc)) p else (prec_acc + p)
    }
    n <- n + 1L
  }

  w <- switch(adapt,
              kenergy = energy_acc / n,
              precision = prec_acc / n,
              kenergy_prec = (energy_acc / n) * (prec_acc / n),
              reliability = stop("adapt = 'reliability' not implemented"))
  w[!is.finite(w) | w < 0] <- 0
  if (!any(w > 0)) w[] <- 1
  .dkge_winsor(w, upper = winsor)
}

.dkge_combine_weights <- function(w_prior,
                                  w_adapt,
                                  combine = "product",
                                  mix = 0.6,
                                  shrink = list(alpha = 0.5, winsor = 0.99, normalize = "mean", roi_smooth = FALSE),
                                  roi_labels = NULL) {
  w_prior <- .dkge_norm_vec(w_prior, "mean")
  w_adapt <- .dkge_norm_vec(w_adapt, "mean")

  w <- switch(combine,
              product = {
                log_w <- (1 - mix) * log(pmax(w_prior, 1e-12)) + mix * log(pmax(w_adapt, 1e-12))
                exp(log_w)
              },
              sum = (1 - mix) * w_prior + mix * w_adapt,
              override_adapt = ifelse(is.finite(w_adapt) & w_adapt > 0, w_adapt, w_prior),
              override_prior = ifelse(is.finite(w_prior) & w_prior > 0, w_prior, w_adapt),
              stop("Unknown combine rule"))

  w <- .dkge_winsor(w, upper = shrink$winsor)

  if (shrink$roi_smooth && !is.null(roi_labels) && length(roi_labels) == length(w)) {
    groups <- as.integer(as.factor(roi_labels))
    for (g in unique(groups)) {
      idx <- which(groups == g)
      if (length(idx)) w[idx] <- stats::median(w[idx], na.rm = TRUE)
    }
  }

  w <- (1 - shrink$alpha) + shrink$alpha * w
  .dkge_norm_vec(w, shrink$normalize)
}

# -------------------------------------------------------------------------
# Public default ----------------------------------------------------------
# -------------------------------------------------------------------------

#' Default DKGE weight specification
#' @return `dkge_weights()` object with k-energy precision combination.
#' @export
dkge_weights_auto <- function() {
  dkge_weights(
    prior = NULL,
    adapt = "kenergy_prec",
    combine = "product",
    mix = 0.6,
    shrink = list(alpha = 0.5, winsor = 0.99, normalize = "mean", roi_smooth = FALSE),
    scope = "fold",
    k_weight = NULL,
    collapse = NULL,
    roi = NULL
  )
}

# -------------------------------------------------------------------------
# High-level evaluators ---------------------------------------------------
# -------------------------------------------------------------------------

.dkge_resolve_voxel_weights <- function(weights,
                                        B_list,
                                        kernel_info,
                                        sigma2_list = NULL) {
  stopifnot(inherits(weights, "dkge_weights"))
  V <- ncol(B_list[[1]])
  w_prior <- .dkge_eval_prior(weights$prior, V)

  K_effects <- NULL
  kernel_failed <- FALSE
  if (weights$adapt %in% c("kenergy", "kenergy_prec")) {
    K_effects <- tryCatch(
      .dkge_weights_resolve_k(weights, kernel_info),
      error = function(e) {
        warning(sprintf("k-energy weights unavailable (%s); using uniform weights", e$message))
        NULL
      }
    )
    if (is.null(K_effects)) kernel_failed <- TRUE
  }

  if (kernel_failed) {
    w_adapt <- rep(1, V)
    w_subject <- NULL
    w_total_subject <- NULL
  } else if (weights$scope == "fold" || length(B_list) <= 1L) {
    w_adapt <- .dkge_adapt_weights(B_list,
                                   adapt = weights$adapt,
                                   K = K_effects,
                                   sigma2_list = sigma2_list,
                                   winsor = weights$shrink$winsor)
    w_subject <- NULL
    w_total_subject <- NULL
  } else {
    w_subject <- lapply(B_list, function(B) {
      .dkge_adapt_weights(list(B),
                          adapt = weights$adapt,
                          K = K_effects,
                          sigma2_list = NULL,
                          winsor = weights$shrink$winsor)
    })
    w_adapt <- Reduce(`+`, w_subject) / length(w_subject)
    w_total_subject <- lapply(w_subject, function(ws) {
      .dkge_combine_weights(w_prior,
                            ws,
                            combine = weights$combine,
                            mix = weights$mix,
                            shrink = weights$shrink,
                            roi_labels = weights$roi)
    })
  }

  if (is.null(w_total_subject)) {
    w_total <- .dkge_combine_weights(w_prior,
                                     w_adapt,
                                     combine = weights$combine,
                                     mix = weights$mix,
                                     shrink = weights$shrink,
                                     roi_labels = weights$roi)
  } else {
    w_total <- Reduce(`+`, w_total_subject) / length(w_total_subject)
  }

  list(prior = w_prior,
       adapt = w_adapt,
       total = w_total,
       total_subject = w_total_subject,
       subject = w_subject,
       kernel = K_effects)
}

.dkge_weight_kernel_payload <- function(K, info) {
  levels <- NULL
  map <- NULL
  if (!is.null(info)) {
    map <- info$map %||% info$info$map %||% NULL
    if (!is.null(info$levels)) {
      levels <- info$levels
    } else if (!is.null(info$info$levels)) {
      levels <- info$info$levels
    }
  }
  list(
    K = as.matrix(K),
    map = map,
    info = list(levels = levels)
  )
}

.dkge_fold_weight_context <- function(fit, train_ids, weight_spec = NULL, ridge = 0) {
  stopifnot(inherits(fit, "dkge"))
  weight_spec <- weight_spec %||% fit$weight_spec %||% dkge_weights(adapt = "none")
  stopifnot(inherits(weight_spec, "dkge_weights"))

  kernel_payload <- .dkge_weight_kernel_payload(fit$K, fit$kernel_info)
  B_train <- fit$Btil[train_ids]
  Omega_train <- fit$Omega[train_ids]
  subject_weights <- fit$weights[train_ids]

  weight_eval <- .dkge_resolve_voxel_weights(weight_spec, B_train, kernel_payload)
  voxel_weights_train <- weight_eval$total_subject %||% weight_eval$total

  accum <- .dkge_accumulate_chat(B_train, Omega_train, fit$Khalf, subject_weights,
                                 voxel_weights = voxel_weights_train)
  Chat <- accum$Chat
  if (ridge > 0) Chat <- Chat + ridge * diag(nrow(Chat))
  Chat <- (Chat + t(Chat)) / 2

  list(
    Chat = Chat,
    weights = weight_eval,
    train_ids = train_ids,
    weight_spec = weight_spec
  )
}

#' Refit a DKGE object with a new voxel-weight specification
#'
#' Convenience helper that rebuilds a DKGE fit using the original inputs but a
#' different [dkge_weights()] specification. The original fit must have been
#' constructed with `keep_inputs = TRUE` (the default for [dkge()]), so the
#' underlying `dkge_data` bundle is available.
#'
#' @param fit A fitted `dkge` object.
#' @param weights A `dkge_weights` specification to apply. When `NULL`, the
#'   existing weight spec stored in `fit$weight_spec` is reused.
#' @return A new `dkge` object with updated voxel weights.
#' @export
dkge_update_weights <- function(fit, weights = NULL) {
  stopifnot(inherits(fit, "dkge"))
  weight_spec <- weights %||% fit$weight_spec %||% dkge_weights(adapt = "none")
  stopifnot(inherits(weight_spec, "dkge_weights"))

  data_bundle <- fit$input
  if (is.null(data_bundle)) {
    stop("fit does not retain its input data (set keep_inputs = TRUE); refit manually with dkge().")
  }

  cpca <- fit$cpca %||% list(part = "none")
  new_fit <- dkge_fit(
    data = data_bundle,
    K = fit$K,
    Omega_list = fit$Omega,
    w_method = fit$w_method %||% "mfa_sigma1",
    w_tau = fit$w_tau %||% 0.3,
    ridge = fit$ridge_input %||% 0,
    rank = fit$rank_requested %||% fit$rank,
    keep_X = !is.null(fit$X_concat),
    cpca_blocks = cpca$blocks %||% NULL,
    cpca_T = cpca$T %||% NULL,
    cpca_part = cpca$part %||% "none",
    cpca_ridge = cpca$ridge %||% 0,
    weights = weight_spec
  )
  new_fit$input <- data_bundle
  new_fit
}

