# test-data.R
# Diagnostic tests for dkge data constructors and high-level entry point

library(testthat)

make_subject_fixture <- function(q = 3, P = 4, T = 20, seed = 101) {
  set.seed(seed)
  beta <- matrix(rnorm(q * P), q, P,
                 dimnames = list(paste0("eff", seq_len(q)), NULL))
  design <- qr.Q(qr(matrix(rnorm(T * q), T, q)))
  colnames(design) <- c("eff2", "eff1", "eff3")[seq_len(q)]
  list(beta = beta, design = design)
}

test_that("dkge_subject aligns effect names and sets defaults", {
  fx <- make_subject_fixture()
  subj <- dkge_subject(fx$beta, design = fx$design)
  expect_s3_class(subj, "dkge_subject")
  expect_equal(rownames(subj$beta), colnames(subj$design))
  expect_equal(subj$n_clusters, ncol(fx$beta))
  expect_true(all(grepl("cluster_", subj$cluster_ids)))
})

test_that("dkge_subject validates omega", {
  fx <- make_subject_fixture()
  expect_error(dkge_subject(fx$beta, fx$design, omega = 1:2), "length", fixed = FALSE)
  expect_error(dkge_subject(fx$beta, fx$design, omega = matrix(1, 2, 2)), "clusters", fixed = FALSE)
  ok <- dkge_subject(fx$beta, fx$design, omega = rep(1, ncol(fx$beta)))
  expect_equal(ok$omega, rep(1, ncol(fx$beta)))
})

test_that("dkge_subject list method and defaults work", {
  fx <- make_subject_fixture()
  beta <- fx$beta
  rownames(beta) <- NULL
  subj <- dkge_subject(list(beta = beta, design = fx$design, id = "id1"))
  expect_equal(subj$id, "id1")
  expect_equal(rownames(subj$beta), colnames(subj$design))
})

test_that("dkge_subject default method errors", {
  expect_error(dkge_subject(1:5), "Unsupported")
})

test_that("dkge_data bundles raw matrices and normalises ids", {
  fx <- make_subject_fixture()
  betas <- replicate(3, fx$beta, simplify = FALSE)
  designs <- replicate(3, fx$design, simplify = FALSE)
  dat <- dkge_data(betas, designs)
  expect_s3_class(dat, "dkge_data")
  expect_equal(length(dat$betas), 3)
  expect_equal(dat$effects, colnames(fx$design))
  expect_true(all(grepl("sub", dat$subject_ids)))
})

test_that("dkge_data aligns partial effect overlaps and records provenance", {
  beta1 <- matrix(1:4, 2, 2, dimnames = list(c("eff1", "eff2"), NULL))
  beta2 <- matrix(5:10, 2, 3, dimnames = list(c("eff2", "eff3"), NULL))
  design1 <- matrix(seq_len(10), 5, 2)
  design2 <- matrix(seq_len(10), 5, 2)
  colnames(design1) <- c("eff1", "eff2")
  colnames(design2) <- c("eff2", "eff3")
  dat <- dkge_data(list(beta1, beta2), list(design1, design2))
  expect_equal(dat$effects, c("eff1", "eff2", "eff3"))
  expect_equal(nrow(dat$betas[[1]]), 3)
  expect_equal(nrow(dat$betas[[2]]), 3)
  prov <- dat$provenance
  expect_true(is.list(prov))
  expect_equal(length(prov$effect_ids), 3)
  expect_false(prov$obs_mask[[1]][3])
  expect_false(prov$obs_mask[[2]][1])
  expect_equal(unname(diag(prov$pair_counts)), c(1L, 2L, 1L))
})

test_that("dkge_data respects provided subject ids and omega", {
  fx <- make_subject_fixture()
  betas <- replicate(2, fx$beta, simplify = FALSE)
  designs <- replicate(2, fx$design, simplify = FALSE)
  omega <- list(rep(1, ncol(fx$beta)), diag(ncol(fx$beta)))
  dat <- dkge_data(betas, designs, omega = omega, subject_ids = c("A", "B"))
  expect_equal(dat$subject_ids, c("A", "B"))
  expect_equal(dat$omega[[1]], omega[[1]])
  expect_equal(dat$omega[[2]], omega[[2]])
})

test_that("dkge_data errors on mismatched effects", {
  fx <- make_subject_fixture()
  betas <- list(fx$beta, fx$beta)
  designs <- list(fx$design, fx$design)
  colnames(designs[[2]]) <- letters[1:ncol(designs[[2]])]
  expect_error(dkge_data(betas, designs), "Row names of beta matrix must match design column names")
})

make_dkge_fixture <- function(S = 3, q = 3, P = 4, T = 20, seed = 202) {
  set.seed(seed)
  betas <- replicate(S, matrix(rnorm(q * P), q, P,
                               dimnames = list(paste0("eff", seq_len(q)), NULL)), simplify = FALSE)
  designs <- replicate(S, {
    X <- matrix(rnorm(T * q), T, q)
    qr.Q(qr(X))
  }, simplify = FALSE)
  for (i in seq_along(designs)) colnames(designs[[i]]) <- paste0("eff", seq_len(q))
  list(betas = betas, designs = designs)
}

test_that("dkge high-level stores inputs when requested", {
  fx <- make_dkge_fixture()
  kernel <- diag(nrow(fx$betas[[1]]))
  fit1 <- dkge(fx$betas, designs = fx$designs, kernel = kernel, keep_inputs = TRUE, rank = 2)
  fit2 <- dkge(fx$betas, designs = fx$designs, kernel = list(K = kernel, info = "identity"),
               keep_inputs = FALSE, rank = 2)
  expect_true(!is.null(fit1$input))
  expect_null(fit2$input)
  expect_equal(fit2$kernel_info, "identity")
})

test_that("dkge accepts dkge_data and omega overrides", {
  fx <- make_dkge_fixture()
  kernel <- diag(nrow(fx$betas[[1]]))
  data_bundle <- dkge_data(fx$betas, fx$designs)
  override <- lapply(fx$betas, function(b) rep(2, ncol(b)))
  fit <- dkge(data_bundle, kernel = kernel, omega = override, rank = 2)
  expect_equal(length(fit$Omega), length(override))
  expect_equal(fit$Omega[[1]], override[[1]])
})

# -------------------------------------------------------------------------
# Ordering invariance tests -----------------------------------------------
# -------------------------------------------------------------------------

test_that("dkge_data produces identical aligned effects regardless of input effect order", {
  withr::local_seed(999)
  # Create two subjects with permuted effect ordering
  effects_order1 <- c("eff1", "eff2", "eff3")
  effects_order2 <- c("eff3", "eff1", "eff2")  # permuted

  beta1 <- matrix(rnorm(3 * 10), 3, 10, dimnames = list(effects_order1, NULL))
  beta2 <- matrix(rnorm(3 * 10), 3, 10, dimnames = list(effects_order2, NULL))

  design1 <- matrix(rnorm(20 * 3), 20, 3, dimnames = list(NULL, effects_order1))
  design2 <- matrix(rnorm(20 * 3), 20, 3, dimnames = list(NULL, effects_order2))

  data <- dkge_data(list(beta1, beta2), list(design1, design2))

  # All outputs should have consistent effect ordering
  expect_identical(rownames(data$betas[[1]]), rownames(data$betas[[2]]))
  expect_identical(colnames(data$designs[[1]]), colnames(data$designs[[2]]))
  expect_identical(data$effects, rownames(data$betas[[1]]))
})

test_that("dkge_data is invariant to subject list ordering", {
  withr::local_seed(888)
  beta_a <- matrix(rnorm(3 * 10), 3, 10, dimnames = list(c("e1", "e2", "e3"), NULL))
  beta_b <- matrix(rnorm(3 * 10), 3, 10, dimnames = list(c("e1", "e2", "e3"), NULL))
  design_a <- matrix(rnorm(20 * 3), 20, 3, dimnames = list(NULL, c("e1", "e2", "e3")))
  design_b <- matrix(rnorm(20 * 3), 20, 3, dimnames = list(NULL, c("e1", "e2", "e3")))

  data1 <- dkge_data(list(beta_a, beta_b), list(design_a, design_b),
                     subject_ids = c("A", "B"))
  data2 <- dkge_data(list(beta_b, beta_a), list(design_b, design_a),
                     subject_ids = c("B", "A"))

  # Effects should be identical
  expect_identical(data1$effects, data2$effects)
  # Subject A's data should be identical in both
  idx1 <- which(data1$subject_ids == "A")
  idx2 <- which(data2$subject_ids == "A")
  expect_equal(data1$betas[[idx1]], data2$betas[[idx2]])
})

test_that("provenance correctly tracks effect coverage with partial overlap", {
  beta1 <- matrix(1:6, 2, 3, dimnames = list(c("e1", "e2"), NULL))
  beta2 <- matrix(1:6, 2, 3, dimnames = list(c("e2", "e3"), NULL))
  design1 <- matrix(1:10, 5, 2, dimnames = list(NULL, c("e1", "e2")))
  design2 <- matrix(1:10, 5, 2, dimnames = list(NULL, c("e2", "e3")))

  data <- dkge_data(list(beta1, beta2), list(design1, design2))
  prov <- data$provenance

  # e1: only subject 1, e2: both, e3: only subject 2
  expect_equal(prov$pair_counts["e1", "e1"], 1L)
  expect_equal(prov$pair_counts["e2", "e2"], 2L)
  expect_equal(prov$pair_counts["e3", "e3"], 1L)
  expect_equal(prov$pair_counts["e1", "e3"], 0L)  # No subject has both
})

test_that("effect alignment produces correct values after reordering", {
  withr::local_seed(777)
  # Create two subjects with known values in specific effect positions
  beta_orig <- matrix(1:9, 3, 3, dimnames = list(c("A", "B", "C"), NULL))
  # Create the same beta with permuted row order
  beta_perm <- beta_orig[c("C", "A", "B"), ]

  # Second subject (normal beta) to meet minimum 2-subject requirement
  beta_other <- matrix(10:18, 3, 3, dimnames = list(c("A", "B", "C"), NULL))

  design_orig <- matrix(rnorm(20 * 3), 20, 3, dimnames = list(NULL, c("A", "B", "C")))
  design_perm <- design_orig[, c("C", "A", "B")]
  colnames(design_perm) <- c("C", "A", "B")

  data1 <- dkge_data(list(beta_orig, beta_other), list(design_orig, design_orig))
  data2 <- dkge_data(list(beta_perm, beta_other), list(design_perm, design_orig))

  # Effect order follows first subject's design column order
  # Both should preserve all effects, though order may differ
  expect_setequal(data1$effects, data2$effects)
  # Values for each effect should be the same regardless of input order
  expect_equal(data1$betas[[1]]["A", ], data2$betas[[1]]["A", ])
  expect_equal(data1$betas[[1]]["B", ], data2$betas[[1]]["B", ])
  expect_equal(data1$betas[[1]]["C", ], data2$betas[[1]]["C", ])
})
