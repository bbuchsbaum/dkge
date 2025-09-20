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
