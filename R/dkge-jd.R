# dkge-jd.R
# Optional joint diagonalisation solver operating in the K^{1/2} metric.

#' Control parameters for DKGE joint diagonalisation
#'
#' @param maxit Maximum number of optimization iterations.
#' @param step Initial step size for the projected gradient descent.
#' @param step_decay Multiplicative decay applied during backtracking.
#' @param step_min Minimum step size allowed before terminating.
#' @param tol Convergence tolerance applied to both gradient norm and off-diagonal energy.
#' @param linesearch Logical; when `TRUE` perform Armijo backtracking.
#' @param armijo Armijo condition constant used during backtracking.
#' @param verbose Logical; emit per-iteration progress when `TRUE`.
#' @param record Logical; store per-iteration diagnostics.
#' @return A list of control parameters.
#' @examples
#' ctrl <- dkge_jd_control(maxit = 10, verbose = FALSE)
#' names(ctrl)
#' @export
dkge_jd_control <- function(maxit = 200L,
                            step = 0.5,
                            step_decay = 0.5,
                            step_min = 1e-6,
                            tol = 1e-7,
                            linesearch = TRUE,
                            armijo = 1e-4,
                            verbose = FALSE,
                            record = FALSE) {
  list(
    maxit = maxit,
    step = step,
    step_decay = step_decay,
    step_min = step_min,
    tol = tol,
    linesearch = linesearch,
    armijo = armijo,
    verbose = verbose,
    record = record
  )
}

.dkge_jd_default_mask <- function(q) {
  mask <- matrix(1, q, q)
  diag(mask) <- 0
  mask
}

.dkge_jd_check_inputs <- function(A_list, weights, mask_list) {
  stopifnot(length(A_list) >= 1)
  q <- nrow(A_list[[1]])
  for (i in seq_along(A_list)) {
    A <- A_list[[i]]
    stopifnot(is.matrix(A), nrow(A) == q, ncol(A) == q)
    if (!isTRUE(all.equal(A, t(A)))) {
      warning("A_list contains a non-symmetric matrix; symmetrising it.")
      A <- 0.5 * (A + t(A))
      A_list[[i]] <- A
    }
  }
  if (is.null(weights)) {
    weights <- rep(1, length(A_list))
  } else {
    stopifnot(length(weights) == length(A_list))
  }
  if (is.null(mask_list)) {
    mask_list <- replicate(length(A_list), .dkge_jd_default_mask(q), simplify = FALSE)
  } else {
    stopifnot(length(mask_list) == length(A_list))
    mask_list <- lapply(mask_list, function(M) {
      if (is.null(M)) {
        M <- .dkge_jd_default_mask(q)
      } else {
        stopifnot(is.matrix(M), nrow(M) == q, ncol(M) == q)
        M <- 0.5 * (M + t(M))
        diag(M) <- 0
      }
      M
    })
  }
  list(q = q, weights = weights, masks = mask_list, A_list = A_list)
}

.dkge_jd_eval <- function(Q, A_list, weights, mask_list, compute_grad = TRUE) {
  q <- nrow(Q)
  sym <- function(M) 0.5 * (M + t(M))

  value_acc <- 0
  diag_energy <- numeric(q)
  Ge <- matrix(0, q, q)

  for (s in seq_along(A_list)) {
    A <- A_list[[s]]
    w <- weights[s]
    mask <- mask_list[[s]]

    AQ <- A %*% Q
    D <- crossprod(Q, AQ) # Q^T A Q
    diag_vec <- diag(D)
    diag_energy <- diag_energy + w * diag_vec

    G <- D * mask
    G <- sym(G)
    diag(G) <- 0

    energy_s <- sum(G * G)
    value_acc <- value_acc + w * energy_s

    if (compute_grad && energy_s > 0) {
      Ge <- Ge + 2 * w * (AQ %*% G)
    }
  }

  value <- 0.5 * value_acc
  offdiag <- sqrt(value_acc)
  grad <- NULL
  grad_norm <- NA_real_

  if (compute_grad) {
    proj <- Ge - Q %*% sym(crossprod(Q, Ge))
    grad <- proj
    grad_norm <- sqrt(sum(proj * proj))
  }

  list(
    value = value,
    offdiag = offdiag,
    grad = grad,
    grad_norm = grad_norm,
    diag_energy = diag_energy
  )
}

.dkge_jd_retract <- function(Y) {
  sv <- svd(Y)
  sv$u %*% t(sv$v)
}

dkge_jd_solve <- function(A_list,
                          weights = NULL,
                          rank = NULL,
                          Q_init = NULL,
                          mask_list = NULL,
                          Chat = NULL,
                          control = dkge_jd_control()) {
  check <- .dkge_jd_check_inputs(A_list, weights, mask_list)
  q <- check$q
  weights <- check$weights
  mask_list <- check$masks
  A_list <- check$A_list
  rank <- if (is.null(rank)) q else max(1L, min(rank, q))

  Q <- if (is.null(Q_init)) diag(q) else {
    stopifnot(is.matrix(Q_init), nrow(Q_init) == q, ncol(Q_init) == q)
    Q_init
  }

  stats <- .dkge_jd_eval(Q, A_list, weights, mask_list)

  history <- if (control$record) {
    data.frame(iter = integer(), value = numeric(), grad = numeric(), offdiag = numeric())
  } else {
    NULL
  }

  step <- control$step
  iter <- 0L

  while (iter < control$maxit) {
    grad <- stats$grad
    grad_norm <- stats$grad_norm
    if ((is.finite(grad_norm) && grad_norm <= control$tol) ||
        stats$offdiag <= control$tol) {
      break
    }

    if (control$record) {
      history <- rbind(history, data.frame(
        iter = iter,
        value = stats$value,
        grad = grad_norm,
        offdiag = stats$offdiag
      ))
    }

    step_try <- step
    accepted <- FALSE

    while (!accepted && step_try >= control$step_min) {
      Q_trial <- .dkge_jd_retract(Q - step_try * grad)
      stats_trial <- .dkge_jd_eval(Q_trial, A_list, weights, mask_list)
      if (!control$linesearch ||
          stats_trial$value <= stats$value - control$armijo * step_try * grad_norm^2) {
        Q <- Q_trial
        stats <- stats_trial
        step <- step_try
        accepted <- TRUE
      } else {
        step_try <- step_try * control$step_decay
      }
    }

    if (!accepted) {
      if (control$verbose) {
        message("JD: failed to find acceptable step; stopping.")
      }
      break
    }

    iter <- iter + 1L
    if (control$verbose) {
      message(sprintf("JD iter %d: value %.3e, offdiag %.3e, grad %.3e",
                      iter, stats$value, stats$offdiag, stats$grad_norm))
    }
  }

  if (is.null(Chat)) {
    Chat <- matrix(0, q, q)
    for (s in seq_along(A_list)) {
      Chat <- Chat + weights[s] * A_list[[s]]
    }
    Chat <- 0.5 * (Chat + t(Chat))
  }

  Chat_Q <- crossprod(Q, Chat %*% Q)
  diag_vals <- diag(Chat_Q)
  order_idx <- order(diag_vals, decreasing = TRUE)
  Q_sorted <- Q[, order_idx, drop = FALSE]
  diag_sorted <- diag_vals[order_idx]
  Chat_Q_sorted <- Chat_Q[order_idx, order_idx, drop = FALSE]

  list(
    Q = Q_sorted,
    Chat_Q = Chat_Q_sorted,
    diag_vals = diag_sorted,
    value = stats$value,
    offdiag = stats$offdiag,
    grad_norm = stats$grad_norm,
    iterations = iter,
    history = history,
    weights = weights
  )
}
