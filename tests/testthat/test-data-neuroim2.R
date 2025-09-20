# test-data-neuroim2.R
# Simplified tests for neuroim2 integration

library(testthat)
library(dkge)

# Helper functions for creating neuroim2 test objects
create_test_neurospace <- function(dims = c(4, 4, 4), spacing = c(1, 1, 1)) {
  neuroim2::NeuroSpace(dims, spacing)
}

create_test_mask <- function(dims = c(4, 4, 4)) {
  sp3 <- create_test_neurospace(dims)
  mask_array <- array(FALSE, dims)
  # Create a spherical mask in the center
  center <- dims / 2
  radius_sq <- (min(dims) / 3)^2
  for (i in 1:dims[1]) {
    for (j in 1:dims[2]) {
      for (k in 1:dims[3]) {
        dist_sq <- sum((c(i, j, k) - center)^2)
        if (dist_sq <= radius_sq) {
          mask_array[i, j, k] <- TRUE
        }
      }
    }
  }
  neuroim2::LogicalNeuroVol(mask_array, sp3)
}

test_that("neuroim2 integration: NeuroVec with pre-computed betas works", {
  skip_if_not_installed("neuroim2")

  # Create test data
  dims <- c(4, 4, 4)
  q <- 3  # Number of design effects
  n_time <- 20

  # Create a NeuroVec with pre-computed betas (3D + effects dimension)
  sp_beta <- neuroim2::NeuroSpace(c(dims, q), c(1, 1, 1))
  beta_array <- array(rnorm(prod(dims) * q), c(dims, q))
  beta_vec <- neuroim2::NeuroVec(beta_array, sp_beta)

  # Create design matrix (just for metadata)
  design <- matrix(rnorm(n_time * q), n_time, q)
  colnames(design) <- paste0("eff", seq_len(q))

  # Build subject (won't compute betas since length(dims) != 4)
  subj <- dkge_subject(beta_vec, design = design, id = "neuro_sub01", compute_betas = FALSE)

  expect_s3_class(subj, "dkge_subject")
  expect_equal(nrow(subj$beta), q)
  expect_equal(ncol(subj$beta), prod(dims))
  expect_equal(subj$id, "neuro_sub01")
  expect_equal(subj$effects, paste0("eff", seq_len(q)))
})

test_that("neuroim2 integration: NeuroVec with mask works", {
  skip_if_not_installed("neuroim2")

  # Create test data
  dims <- c(4, 4, 4)
  q <- 3
  n_time <- 20

  # Create mask
  mask <- create_test_mask(dims)
  n_voxels <- sum(as.logical(neuroim2::values(mask)))

  # Create a NeuroVec with pre-computed betas
  sp_beta <- neuroim2::NeuroSpace(c(dims, q), c(1, 1, 1))
  beta_array <- array(rnorm(prod(dims) * q), c(dims, q))
  beta_vec <- neuroim2::NeuroVec(beta_array, sp_beta)

  # Create design matrix (just for metadata)
  design <- matrix(rnorm(n_time * q), n_time, q)
  colnames(design) <- paste0("eff", seq_len(q))

  # Build subject with mask
  subj <- dkge_subject(beta_vec, design = design, id = "masked_sub", mask = mask, compute_betas = FALSE)

  expect_s3_class(subj, "dkge_subject")
  expect_equal(nrow(subj$beta), q)
  expect_equal(ncol(subj$beta), n_voxels)  # Should only include masked voxels
  expect_equal(subj$id, "masked_sub")
})

test_that("neuroim2 integration: dkge_data bundles neuroim2 subjects", {
  skip_if_not_installed("neuroim2")

  n_subjects <- 3
  dims <- c(4, 4, 4)
  q <- 3
  n_time <- 20

  subjects <- list()
  for (s in seq_len(n_subjects)) {
    # Create a NeuroVec with pre-computed betas
    sp_beta <- neuroim2::NeuroSpace(c(dims, q), c(1, 1, 1))
    beta_array <- array(rnorm(prod(dims) * q, sd = s), c(dims, q))
    beta_vec <- neuroim2::NeuroVec(beta_array, sp_beta)

    # Create design matrix
    design <- matrix(rnorm(n_time * q), n_time, q)
    colnames(design) <- paste0("eff", seq_len(q))

    subjects[[s]] <- dkge_subject(beta_vec, design = design, id = paste0("sub", s), compute_betas = FALSE)
  }

  # Bundle subjects
  bundle <- dkge_data(subjects)

  expect_s3_class(bundle, "dkge_data")
  expect_equal(bundle$n_subjects, n_subjects)
  expect_equal(bundle$effects, paste0("eff", seq_len(q)))
  expect_equal(bundle$subject_ids, paste0("sub", seq_len(n_subjects)))
})