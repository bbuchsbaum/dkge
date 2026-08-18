# dkge-align-data.R
# Utilities for aligning subject betas/designs to a shared effect union

#' Align dkge_subject records to a shared union of effects
#'
#' Detects partial effect overlap across subjects, embeds each subject's beta
#' and design matrices into the shared union (filling missing effects with
#' zeros), and returns provenance metadata describing observed coverage.
#'
#' @param subjects List of `dkge_subject` objects.
#' @return List with aligned `subjects` and `provenance` (effect IDs, per-subject
#'   observation masks, pairwise coverage counts, and coverage summary).
#' @keywords internal
#' @noRd
.dkge_align_subjects_to_union <- function(subjects) {
  stopifnot(length(subjects) > 0L)

  effect_list <- lapply(subjects, `[[`, "effects")
  subject_ids <- vapply(subjects, function(s) s$id %||% "", character(1))

  # A partial-coverage union must not depend on subject or local row order.
  # Use a deterministic label order; named kernels are subsequently aligned to
  # these labels in .dkge_fit_prepare().
  union_ids <- sort(unique(unlist(effect_list, use.names = FALSE)),
                    method = "radix")
  union_ids <- as.character(union_ids)
  n_union <- length(union_ids)

  embed_beta <- function(beta, current_ids) {
    stopifnot(nrow(beta) == length(current_ids))
    out <- matrix(0, n_union, ncol(beta),
                  dimnames = list(union_ids, colnames(beta)))
    idx <- match(current_ids, union_ids)
    out[idx, ] <- beta
    out
  }

  embed_design <- function(design, current_ids) {
    stopifnot(ncol(design) == length(current_ids))
    out <- matrix(0, nrow(design), n_union,
                  dimnames = list(rownames(design), union_ids))
    idx <- match(current_ids, union_ids)
    out[, idx] <- design
    out
  }

  embed_effect_vector <- function(x, current_ids, fill = 0) {
    if (is.null(x)) return(NULL)
    out <- rep(fill, n_union)
    names(out) <- union_ids
    out[match(current_ids, union_ids)] <- as.numeric(x)
    out
  }

  embed_effect_covariance <- function(x, current_ids) {
    if (is.null(x)) return(NULL)
    out <- matrix(0, n_union, n_union,
                  dimnames = list(union_ids, union_ids))
    idx <- match(current_ids, union_ids)
    out[idx, idx] <- x
    out
  }

  obs_mask <- vector("list", length(subjects))
  observed_rows <- vector("list", length(subjects))
  names(obs_mask) <- subject_ids
  names(observed_rows) <- subject_ids

  pair_counts <- matrix(0L, n_union, n_union,
                        dimnames = list(union_ids, union_ids))

  for (i in seq_along(subjects)) {
    subj <- subjects[[i]]
    current_ids <- as.character(subj$effects)
    local_observed <- rep(FALSE, length(current_ids))
    local_observed[subj$observed_rows %||% seq_along(current_ids)] <- TRUE
    observed_ids <- current_ids[local_observed]
    mask <- union_ids %in% observed_ids
    names(mask) <- union_ids
    obs_mask[[i]] <- mask

    idx <- which(mask)
    observed_rows[[i]] <- unname(idx)
    if (length(idx)) {
      pair_counts[idx, idx] <- pair_counts[idx, idx] + 1L
    }

    beta_local <- subj$beta
    design_local <- subj$design
    if (any(!local_observed)) {
      beta_local[!local_observed, ] <- 0
      design_local[, !local_observed] <- 0
    }
    subjects[[i]]$beta <- embed_beta(beta_local, current_ids)
    subjects[[i]]$design <- embed_design(design_local, current_ids)
    subjects[[i]]$effects <- union_ids
    subjects[[i]]$observed_rows <- unname(idx)
    subjects[[i]]$effect_n <- embed_effect_vector(subj$effect_n, current_ids)
    subjects[[i]]$effect_precision <- embed_effect_vector(subj$effect_precision,
                                                          current_ids)
    subjects[[i]]$effect_noise_cov <- embed_effect_covariance(
      subj$effect_noise_cov, current_ids
    )
    subjects[[i]]$effect_information <- embed_effect_covariance(
      subj$effect_information, current_ids
    )
    subjects[[i]]$effect_score <- if (is.null(subj$effect_score)) NULL else {
      embed_beta(subj$effect_score, current_ids)
    }
    subjects[[i]]$split_betas <- if (is.null(subj$split_betas)) NULL else {
      lapply(subj$split_betas, function(B) {
        if (any(!local_observed)) B[!local_observed, ] <- 0
        embed_beta(B, current_ids)
      })
    }
    subjects[[i]]$split_provenance <- subj$split_provenance
    if (!is.null(subjects[[i]]$split_provenance$effect_counts)) {
      counts_local <- subjects[[i]]$split_provenance$effect_counts
      counts_union <- matrix(
        0, n_union, ncol(counts_local),
        dimnames = list(union_ids, colnames(counts_local))
      )
      counts_union[match(current_ids, union_ids), ] <- counts_local
      subjects[[i]]$split_provenance$effect_counts <- counts_union
    }
    subjects[[i]]$split_reliability <- subj$split_reliability
    if (!is.null(subjects[[i]]$split_reliability)) {
      for (field in c("correlation", "reliability", "min_half_count",
                      "count_factor")) {
        value <- subjects[[i]]$split_reliability[[field]]
        if (!is.null(value)) {
          fill <- if (identical(field, "min_half_count")) 0 else NA_real_
          subjects[[i]]$split_reliability[[field]] <-
            embed_effect_vector(value, current_ids, fill = fill)
        }
      }
    }
    subjects[[i]]$beta <- `rownames<-`(subjects[[i]]$beta, union_ids)
    subjects[[i]]$design <- `colnames<-`(subjects[[i]]$design, union_ids)
  }

  subj_counts <- as.integer(diag(pair_counts))
  coverage <- data.frame(
    effect = union_ids,
    train_subjects = subj_counts,
    subjects = subj_counts,
    stringsAsFactors = FALSE
  )

  # Warn about sparse subjects (>50% missing effects)
  for (i in seq_along(obs_mask)) {
    mask <- obs_mask[[i]]
    pct_present <- 100 * sum(mask) / length(mask)
    if (pct_present < 50) {
      warning(sprintf(
        "Subject '%s': sparse effect coverage (%.1f%% of effects present). Consider excluding.",
        subject_ids[i], pct_present
      ), call. = FALSE)
    }
  }

  list(
    subjects = subjects,
    provenance = list(
      effect_ids = union_ids,
      obs_mask = obs_mask,
      observed_rows = observed_rows,
      pair_counts = pair_counts,
      coverage = coverage
    )
  )
}

#' Build provenance summary for subjects with identical effect labels
#'
#' @param subjects List of `dkge_subject` objects with identical effect sets.
#' @return Provenance list matching the structure returned by
#'   `.dkge_align_subjects_to_union()`.
#' @keywords internal
#' @noRd
.dkge_full_coverage_provenance <- function(subjects) {
  stopifnot(length(subjects) > 0L)
  effect_ids <- as.character(subjects[[1]]$effects)
  subject_ids <- vapply(subjects, function(s) s$id %||% "", character(1))
  n <- length(effect_ids)

  obs_mask <- vector("list", length(subjects))
  observed_rows <- vector("list", length(subjects))
  names(obs_mask) <- subject_ids
  names(observed_rows) <- subject_ids
  pair_counts <- matrix(0L, n, n, dimnames = list(effect_ids, effect_ids))

  for (i in seq_along(subjects)) {
    idx <- subjects[[i]]$observed_rows %||% seq_len(n)
    mask <- rep(FALSE, n)
    mask[idx] <- TRUE
    names(mask) <- effect_ids
    obs_mask[[i]] <- mask
    observed_rows[[i]] <- unname(which(mask))
    if (length(idx)) {
      pair_counts[idx, idx] <- pair_counts[idx, idx] + 1L
    }
  }

  subj_counts <- as.integer(diag(pair_counts))
  coverage <- data.frame(
    effect = effect_ids,
    train_subjects = subj_counts,
    subjects = subj_counts,
    stringsAsFactors = FALSE
  )

  list(
    effect_ids = effect_ids,
    obs_mask = obs_mask,
    observed_rows = observed_rows,
    pair_counts = pair_counts,
    coverage = coverage
  )
}
