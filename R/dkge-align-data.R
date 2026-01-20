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

  # Determine union while preserving the order in which effects first appear
  union_ids <- character(0)
  for (eff in effect_list) {
    eff <- as.character(eff)
    new_ids <- eff[!eff %in% union_ids]
    union_ids <- c(union_ids, new_ids)
  }
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

  obs_mask <- vector("list", length(subjects))
  names(obs_mask) <- subject_ids

  pair_counts <- matrix(0L, n_union, n_union,
                        dimnames = list(union_ids, union_ids))

  for (i in seq_along(subjects)) {
    subj <- subjects[[i]]
    current_ids <- as.character(subj$effects)
    mask <- union_ids %in% current_ids
    names(mask) <- union_ids
    obs_mask[[i]] <- mask

    idx <- which(mask)
    if (length(idx)) {
      pair_counts[idx, idx] <- pair_counts[idx, idx] + 1L
    }

    subjects[[i]]$beta <- embed_beta(subj$beta, current_ids)
    subjects[[i]]$design <- embed_design(subj$design, current_ids)
    subjects[[i]]$effects <- union_ids
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
      pair_counts = pair_counts,
      coverage = coverage
    )
  )
}

#' Build provenance summary assuming full effect coverage
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

  mask_template <- rep(TRUE, n)
  names(mask_template) <- effect_ids
  obs_mask <- setNames(rep(list(mask_template), length(subjects)), subject_ids)

  pair_counts <- matrix(as.integer(length(subjects)), n, n,
                        dimnames = list(effect_ids, effect_ids))

  subj_counts <- rep.int(length(subjects), n)
  coverage <- data.frame(
    effect = effect_ids,
    train_subjects = subj_counts,
    subjects = subj_counts,
    stringsAsFactors = FALSE
  )

  list(
    effect_ids = effect_ids,
    obs_mask = obs_mask,
    pair_counts = pair_counts,
    coverage = coverage
  )
}
