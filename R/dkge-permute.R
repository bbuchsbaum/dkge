# dkge-permute.R
# Permutation helpers for DKGE classification.

.dkge_shuffle_within_groups <- function(values, groups) {
  stopifnot(length(values) == length(groups))
  groups <- as.factor(groups)
  shuffled <- values
  for (g in levels(groups)) {
    idx <- which(groups == g)
    if (length(idx) > 1) {
      shuffled[idx] <- sample(values[idx], length(idx))
    }
  }
  shuffled
}

.dkge_permute_labels <- function(labels,
                                 scope = c("within_subject", "global"),
                                 subjects = NULL,
                                 blocks = NULL) {
  scope <- match.arg(scope)
  if (scope == "global") {
    return(sample(labels, length(labels)))
  }
  if (is.null(subjects)) {
    stop("subjects must be supplied for within-subject permutations")
  }
  if (!is.null(blocks)) {
    if (length(blocks) != length(labels)) {
      stop("blocks vector must match length of labels")
    }
    grp <- interaction(subjects, blocks, drop = TRUE)
    .dkge_shuffle_within_groups(labels, grp)
  } else {
    .dkge_shuffle_within_groups(labels, subjects)
  }
}

.dkge_sign_flips <- function(groups) {
  if (is.null(groups)) {
    return(sample(c(-1, 1), 1, replace = TRUE))
  }
  uniq <- unique(groups)
  flips <- sample(c(-1, 1), length(uniq), replace = TRUE)
  setNames(flips, uniq)
}
