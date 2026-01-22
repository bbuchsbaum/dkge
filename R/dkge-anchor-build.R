# dkge-anchor-build.R
# Anchor selection strategies and kernel alignment utilities for feature-anchored DKGE.

#' Pairwise RBF kernel between feature matrices
#'
#' Computes \eqn{k(x, y) = \exp(-\|x - y\|^2 / (2\sigma^2))} using an optional
#' bandwidth. When `sigma` is `NULL`, the median distance heuristic is applied on
#' the rows of `X`.
#'
#' @keywords internal
#' @noRd
.dKGE_anchor_rbf_kernel <- function(X, Y = NULL, sigma = NULL) {
  stopifnot(is.matrix(X))
  if (is.null(Y)) {
    Y <- X
  } else {
    stopifnot(is.matrix(Y), ncol(Y) == ncol(X))
  }

  if (is.null(sigma)) {
    sigma <- .dkge_anchor_median_distance(X)
  }
  sigma <- as.numeric(sigma)
  if (!is.finite(sigma) || sigma <= 0) sigma <- 1

  XX <- rowSums(X^2)
  YY <- rowSums(Y^2)
  D2 <- outer(XX, YY, "+") - 2 * (X %*% t(Y))
  exp(-D2 / (2 * sigma^2))
}

#' Fast median pairwise distance heuristic
#'
#' Uses an optional random subsample to keep cost manageable for large pools.
#'
#' @keywords internal
#' @noRd
.dkge_anchor_median_distance <- function(X, max_samples = 2000L) {
  n <- nrow(X)
  if (n <= 1L) {
    return(1)
  }
  if (n > max_samples) {
    idx <- sample.int(n, max_samples)
    X <- X[idx, , drop = FALSE]
  }
  dists <- stats::dist(X)
  med <- stats::median(dists)
  if (!is.finite(med) || med <= 0) 1 else med
}

#' Greedy MAP k-DPP selection on an RBF kernel
#'
#' @keywords internal
#' @noRd
.dkge_anchor_dpp_greedy <- function(K, L0, min_gain = 1e-12) {
  n <- nrow(K)
  stopifnot(ncol(K) == n, L0 >= 1L, L0 <= n)
  diag_resid <- pmax(diag(K), 0)
  selected <- integer(0)
  U <- matrix(0, n, L0)
  gains <- numeric(L0)

  for (t in seq_len(L0)) {
    j <- which.max(diag_resid)
    gain <- diag_resid[j]
    gains[t] <- gain
    if (gain < min_gain) {
      gains <- gains[seq_len(t - 1L)]
      break
    }
    selected <- c(selected, j)
    v <- K[, j]
    if (t > 1L) {
      coeff <- matrix(U[j, seq_len(t - 1L), drop = FALSE], ncol = 1)
      proj <- U[, seq_len(t - 1L), drop = FALSE] %*% coeff
      v <- v - proj
    }
    norm_v <- sqrt(sum(v * v))
    if (norm_v <= 0) {
      diag_resid[j] <- 0
      next
    }
    U[, t] <- v / norm_v
    diag_resid <- pmax(diag_resid - U[, t]^2, 0)
  }

  list(indices = unique(selected), gains = gains[gains > 0])
}

#' Farthest-first coverage fill from a seed set
#'
#' @keywords internal
#' @noRd
.dkge_anchor_fill_farthest <- function(X, seed_idx, L) {
  stopifnot(is.matrix(X), length(seed_idx) >= 1L)
  n <- nrow(X)
  d <- ncol(X)
  sel <- unique(pmax(1L, pmin(seed_idx, n)))
  if (!length(sel)) sel <- sample.int(n, 1L)
  sqmin <- rep(Inf, n)
  for (j in sel) {
    diff <- X - matrix(X[j, ], n, d, byrow = TRUE)
    sqmin <- pmin(sqmin, rowSums(diff * diff))
  }
  while (length(sel) < L) {
    j <- which.max(sqmin)
    sel <- c(sel, j)
    diff <- X - matrix(X[j, ], n, d, byrow = TRUE)
    sqmin <- pmin(sqmin, rowSums(diff * diff))
    sqmin[sel] <- 0
  }
  unique(sel)
}

#' k-means++ style D^2 sampling fill
#'
#' @keywords internal
#' @noRd
.dkge_anchor_fill_kpp <- function(X, seed_idx, L, seed = NULL) {
  stopifnot(is.matrix(X))
  n <- nrow(X)
  d <- ncol(X)
  sel <- unique(seed_idx)
  if (!length(sel)) sel <- sample.int(n, 1L)
  if (!is.null(seed)) set.seed(seed)
  sqmin <- rep(Inf, n)
  for (j in sel) {
    diff <- X - matrix(X[j, ], n, d, byrow = TRUE)
    sqmin <- pmin(sqmin, rowSums(diff * diff))
  }
  while (length(sel) < L) {
    prob <- sqmin
    prob[sel] <- 0
    tot <- sum(prob)
    if (!is.finite(tot) || tot <= 0) {
      cand <- setdiff(seq_len(n), sel)
      if (!length(cand)) break
      j <- sample(cand, 1L)
    } else {
      prob <- prob / tot
      j <- sample.int(n, 1L, prob = prob)
    }
    sel <- c(sel, j)
    diff <- X - matrix(X[j, ], n, d, byrow = TRUE)
    sqmin <- pmin(sqmin, rowSums(diff * diff))
    sqmin[sel] <- 0
  }
  unique(sel)
}

#' Determinantal k-means++ (d-kpp) anchor selection
#'
#' Combines a greedy k-DPP seeding phase with either farthest-first or D^2
#' sampling to reach `L` anchors. The kernel bandwidth defaults to the median
#' distance heuristic computed on the training pool.
#'
#' @param X_train Matrix of training feature vectors pooled across subjects.
#' @param L Target number of anchors.
#' @param rho Fraction of anchors chosen by the DPP stage (0 < rho <= 1).
#' @param sigma Optional RBF bandwidth. When `NULL`, the median heuristic is
#'   applied on `X_train`.
#' @param fill Strategy for completing the anchor set after the DPP stage:
#'   `"kcenter"` (farthest-first) or `"kmeanspp"`.
#' @param seed Integer seed for reproducibility.
#' @param min_gain Minimum marginal gain tolerated during the DPP phase.
#'
#' @return List with `indices` (row indices in `X_train`), `sigma`, and
#'   selection diagnostics.
#' @examples
#' X_train <- matrix(rnorm(20 * 3), 20, 3)
#' sel <- dkpp_select_anchors(X_train, L = 6, seed = 1)
#' length(sel$indices)
#' @export
dkpp_select_anchors <- function(X_train,
                                L = 128,
                                rho = 0.5,
                                sigma = NULL,
                                fill = c("kcenter", "kmeanspp"),
                                seed = 1,
                                min_gain = 1e-12) {
  stopifnot(is.matrix(X_train), nrow(X_train) >= 1L)
  fill <- match.arg(fill)
  L <- as.integer(L)
  L <- max(1L, min(L, nrow(X_train)))
  rho <- min(1, max(rho, 0))
  if (is.null(sigma)) {
    sigma <- .dkge_anchor_median_distance(X_train)
  }
  if (!is.null(seed)) set.seed(seed)

  K <- .dKGE_anchor_rbf_kernel(X_train, sigma = sigma)
  L0 <- max(1L, min(L, round(rho * L)))
  dpp <- .dkge_anchor_dpp_greedy(K, L0, min_gain = min_gain)
  seeds <- dpp$indices
  if (!length(seeds)) {
    seeds <- sample.int(nrow(X_train), 1L)
  }
  fill_idx <- if (fill == "kcenter") {
    .dkge_anchor_fill_farthest(X_train, seeds, L)
  } else {
    # Offset seed to decorrelate the kmeans++ stage from the greedy DPP seed
    .dkge_anchor_fill_kpp(X_train, seeds, L, seed = seed + 13L)
  }
  list(indices = unique(fill_idx), sigma = sigma, seeds = seeds, gains = dpp$gains)
}

#' Fold-aware anchor kernel construction
#'
#' Projects per-subject item kernels onto a shared anchor basis derived from a
#' pooled training feature set inside each cross-validation fold. The resulting
#' aligned kernels all have dimension `L x L` with identical ordering.
#'
#' @param features_list List of length `S`; element `s` is an `n_s x d` feature
#'   matrix for subject `s`.
#' @param K_item_list List of length `S`; element `s` is an `n_s x n_s` PSD item
#'   kernel aligned with `features_list[[s]]`.
#' @param folds Optional fold specification. Accepts `NULL` (single context), a
#'   list of held-out subject indices, or a `dkge_folds` object.
#' @param L Number of anchors.
#' @param method Anchor selection strategy (`"kcenter"`, `"kmeanspp"`,
#'   `"random"`, or `"dkpp"`).
#' @param seed Integer seed applied per fold (incremented by fold index).
#' @param sigma Optional bandwidth; when `NULL`, fold-specific heuristics are
#'   used.
#' @param rho Fraction used by the DPP stage when `method = "dkpp"`.
#' @param fill Completion strategy for `method = "dkpp"`.
#' @param center Logical; when `TRUE`, subtract column means of the feature
#'   response matrix before forming the projection.
#' @param whiten Logical; apply whitening with the anchor Gram inverse square
#'   root.
#' @param eps Diagonal jitter used during whitening.
#' @param unit_trace Logical; trace-normalise each subject kernel to maintain
#'   comparable scale.
#' @param item_weights Optional list of numeric vectors providing per-item
#'   weights (`length == n_s`).
#'
#' @return A list indexed by folds; each element contains the anchors, bandwidth,
#'   aligned kernels, anchor identifiers, and the fold's train/test indices.
#' @examples
#' set.seed(1)
#' features_list <- list(
#'   s1 = matrix(rnorm(30 * 5), 30, 5),
#'   s2 = matrix(rnorm(40 * 5), 40, 5),
#'   s3 = matrix(rnorm(35 * 5), 35, 5)
#' )
#' K_item_list <- lapply(features_list, function(X) {
#'   Z <- matrix(rnorm(nrow(X) * 4), nrow(X), 4)
#'   tcrossprod(Z)
#' })
#' built <- dkge_build_anchor_kernels(features_list, K_item_list, L = 8, method = "dkpp")
#' dim(built[[1]]$K_aligned[[1]])
#' @export
dkge_build_anchor_kernels <- function(features_list,
                                       K_item_list,
                                       folds = NULL,
                                       L = 128,
                                       method = c("kcenter", "kmeanspp", "random", "dkpp"),
                                       seed = 1,
                                       sigma = NULL,
                                       rho = 0.5,
                                       fill = c("kcenter", "kmeanspp"),
                                       center = TRUE,
                                       whiten = TRUE,
                                       eps = 1e-6,
                                       unit_trace = TRUE,
                                       item_weights = NULL) {
  method <- match.arg(method)
  fill <- match.arg(fill)
  stopifnot(is.list(features_list), is.list(K_item_list),
            length(features_list) == length(K_item_list))
  S <- length(features_list)
  if (!is.null(item_weights)) {
    stopifnot(is.list(item_weights), length(item_weights) == S)
  }

  subject_ids <- names(features_list)
  if (is.null(subject_ids)) {
    subject_ids <- paste0("subject", seq_len(S))
  }

  for (s in seq_len(S)) {
    Fs <- features_list[[s]]
    Ks <- K_item_list[[s]]
    stopifnot(is.matrix(Fs), nrow(Fs) == nrow(Ks), nrow(Ks) == ncol(Ks))
    storage.mode(Fs) <- "double"
    storage.mode(Ks) <- "double"
  }

  fold_ctx <- .dkge_anchor_expand_folds(folds, S)
  out <- vector("list", length(fold_ctx))
  names(out) <- names(fold_ctx)

  anchor_ids <- function(L) paste0("anchor_", seq_len(L))

  for (f in seq_along(fold_ctx)) {
    test_idx <- fold_ctx[[f]]$test_idx
    train_idx <- fold_ctx[[f]]$train_idx
    if (!length(train_idx)) {
      stop("Each fold must include at least one training subject.", call. = FALSE)
    }

    X_train <- do.call(rbind, features_list[train_idx])
    if (!is.matrix(X_train)) {
      X_train <- matrix(X_train, ncol = ncol(features_list[[train_idx[[1]]]]))
    }

    select_seed <- seed + f - 1L
    sel <- switch(method,
      dkpp = dkpp_select_anchors(X_train, L = L, rho = rho, sigma = sigma,
                                 fill = fill, seed = select_seed),
      kcenter = {
        idx0 <- sample.int(nrow(X_train), 1L)
        idx <- .dkge_anchor_fill_farthest(X_train, idx0, L)
        list(indices = idx, sigma = sigma %||% .dkge_anchor_bandwidth(X_train, X_train[idx, , drop = FALSE]))
      },
      kmeanspp = {
        km <- stats::kmeans(X_train, centers = L, nstart = 10, iter.max = 100)
        list(indices = integer(0), sigma = sigma %||% .dkge_anchor_bandwidth(X_train, km$centers), centers = km$centers)
      },
      random = {
        idx <- sample.int(nrow(X_train), min(L, nrow(X_train)))
        list(indices = idx, sigma = sigma %||% .dkge_anchor_bandwidth(X_train, X_train[idx, , drop = FALSE]))
      }
    )

    if (!is.null(sel$centers)) {
      Z <- sel$centers
    } else {
      idx <- unique(sel$indices)
      if (length(idx) < L) {
        filler <- setdiff(seq_len(nrow(X_train)), idx)
        if (length(filler)) {
          idx <- c(idx, sample(filler, min(L - length(idx), length(filler))))
        }
      }
      idx <- unique(idx)[seq_len(min(length(idx), L))]
      Z <- X_train[idx, , drop = FALSE]
    }
    if (nrow(Z) < L) {
      Z <- rbind(Z, matrix(Z[rep(nrow(Z), L - nrow(Z)), ], nrow = L - nrow(Z), ncol = ncol(Z), byrow = FALSE))
    }

    sigma_f <- sel$sigma
    if (is.null(sigma_f)) {
      sigma_f <- .dkge_anchor_bandwidth(X_train, Z)
    }
    sigma_f <- as.numeric(sigma_f)
    if (!is.finite(sigma_f) || sigma_f <= 0) sigma_f <- 1

    A <- .dKGE_anchor_rbf_kernel(Z, sigma = sigma_f)
    if (whiten) {
      egA <- eigen(A + eps * diag(nrow(A)), symmetric = TRUE)
      vals <- pmax(egA$values, eps)
      A_mhalf <- egA$vectors %*% (diag(1 / sqrt(vals), nrow = length(vals))) %*% t(egA$vectors)
    } else {
      A_mhalf <- diag(1, nrow(A))
    }

    K_aligned <- vector("list", S)
    names(K_aligned) <- subject_ids

    for (s in seq_len(S)) {
      Fs <- features_list[[s]]
      Ks <- K_item_list[[s]]
      Cs <- .dKGE_anchor_rbf_kernel(Fs, Z, sigma = sigma_f)
      if (center) {
        Cs <- sweep(Cs, 2L, colMeans(Cs), FUN = "-")
      }
      if (!is.null(item_weights)) {
        ws <- item_weights[[s]]
        if (!is.null(ws)) {
          stopifnot(length(ws) == nrow(Cs))
          Cs <- Cs * sqrt(pmax(as.numeric(ws), 0))
        }
      }
      Phi <- Cs %*% A_mhalf
      if (!is.null(Ks)) {
        K_anchor <- t(Phi) %*% Ks %*% Phi
      } else {
        K_anchor <- t(Phi) %*% Phi
      }
      K_anchor <- 0.5 * (K_anchor + t(K_anchor))
      if (unit_trace) {
        tr <- sum(diag(K_anchor))
        if (tr > 1e-12) {
          K_anchor <- K_anchor / (tr / ncol(K_anchor))
        }
      }
      K_aligned[[s]] <- K_anchor
    }

    out[[f]] <- list(
      anchors = Z,
      sigma = sigma_f,
      A = A,
      K_aligned = K_aligned,
      anchor_ids = anchor_ids(nrow(Z)),
      train_idx = train_idx,
      test_idx = test_idx,
      selection = list(method = method, seed = select_seed, rho = if (method == "dkpp") rho else NA_real_, fill = fill)
    )
  }
  out
}

#' Compute per-row bandwidth relative to anchors
#'
#' @keywords internal
#' @noRd
.dkge_anchor_bandwidth <- function(X, Z) {
  if (!nrow(Z) || !nrow(X)) {
    return(1)
  }
  XX <- rowSums(X^2)
  ZZ <- rowSums(Z^2)
  D2 <- outer(XX, ZZ, "+") - 2 * (X %*% t(Z))
  nearest <- apply(D2, 1L, min)
  med <- stats::median(sqrt(pmax(nearest, 0)))
  if (!is.finite(med) || med <= 0) 1 else med
}

#' Expand user folds into train/test index lists
#'
#' @keywords internal
#' @noRd
.dkge_anchor_expand_folds <- function(folds, S) {
  if (is.null(folds)) {
    return(list(context1 = list(train_idx = seq_len(S), test_idx = integer(0))))
  }

  make_entry <- function(test_idx, nm = NULL) {
    test_idx <- sort(unique(as.integer(test_idx)))
    test_idx <- test_idx[!is.na(test_idx)]
    test_idx <- test_idx[test_idx >= 1L & test_idx <= S]
    list(train_idx = setdiff(seq_len(S), test_idx), test_idx = test_idx, name = nm)
  }

  if (inherits(folds, "dkge_folds")) {
    if (!is.null(folds$assignments)) {
      assign_names <- names(folds$assignments)
      entries <- lapply(seq_along(folds$assignments), function(i) {
        nm <- if (is.null(assign_names)) NULL else assign_names[[i]]
        make_entry(folds$assignments[[i]], nm)
      })
    } else if (!is.null(folds$folds)) {
      fold_names <- names(folds$folds)
      entries <- lapply(seq_along(folds$folds), function(i) {
        nm <- if (is.null(fold_names)) NULL else fold_names[[i]]
        holdout <- folds$folds[[i]]$subjects %||% integer(0)
        make_entry(holdout, nm)
      })
    } else {
      stop("dkge_folds object missing assignments for anchor expansion.", call. = FALSE)
    }
  } else if (is.list(folds)) {
    fold_names <- names(folds)
    entries <- lapply(seq_along(folds), function(i) {
      nm <- if (is.null(fold_names)) NULL else fold_names[[i]]
      make_entry(folds[[i]], nm)
    })
  } else {
    stop("Unsupported fold specification for dkge_build_anchor_kernels().", call. = FALSE)
  }

  res <- lapply(entries, function(entry) {
    list(train_idx = entry$train_idx, test_idx = entry$test_idx)
  })
  name_vec <- vapply(entries, function(entry) {
    nm <- entry$name
    if (length(nm) && !is.na(nm) && nzchar(nm)) nm else NA_character_
  }, character(1))
  names(res) <- name_vec

  missing_names <- is.na(name_vec) | name_vec == ""
  if (any(missing_names)) {
    default_names <- paste0("fold", seq_along(res))
    names(res)[missing_names] <- default_names[missing_names]
  }
  res
}
