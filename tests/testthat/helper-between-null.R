# Helper utilities for between-subject permutation null-calibration tests.

dkge_test_between_null_experiment <- function(n = 18,
                                              p = 24,
                                              rank = 2,
                                              B = 149,
                                              seed = 123) {
  stopifnot(n >= 8L, p >= 3L, B >= 19L)
  set.seed(seed)

  dat <- data.frame(
    subject_id = paste0("s", seq_len(n)),
    group = factor(rep(c("A", "B"), length.out = n)),
    trait = scale(rnorm(n), center = TRUE, scale = FALSE)[, 1],
    age = scale(rnorm(n), center = TRUE, scale = FALSE)[, 1]
  )
  design <- dkge_subject_model(~ group * trait + age, dat)

  latent <- matrix(rnorm(n * 2L), n, 2L)
  pattern <- matrix(rnorm(2L * p), 2L, p)
  Y <- latent %*% pattern + matrix(rnorm(n * p, sd = 0.35), n, p)
  rownames(Y) <- dat$subject_id
  colnames(Y) <- paste0("f", seq_len(p))

  target <- dkge_make_target(Y = Y, subject_ids = dat$subject_id)
  fit <- dkge_between_rrr(target, design, rank = min(rank, ncol(design$X), ncol(Y)))
  perm <- dkge_between_permute(
    fit,
    terms = c("group", "trait", "group:trait"),
    B = B,
    seed = seed + 10000L
  )

  stats::setNames(perm$summary$p, perm$summary$term)
}

dkge_test_between_null_uniformity <- function(nrep = 30,
                                              n = 18,
                                              p = 24,
                                              rank = 2,
                                              B = 149,
                                              seed = 123) {
  pvals <- vapply(
    seq_len(nrep),
    function(i) {
      dkge_test_between_null_experiment(
        n = n,
        p = p,
        rank = rank,
        B = B,
        seed = seed + i
      )
    },
    FUN.VALUE = c(group = 0, trait = 0, "group:trait" = 0)
  )

  as.data.frame(t(pvals))
}
