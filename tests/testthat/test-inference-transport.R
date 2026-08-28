# test-inference-transport.R
# Ensure inference stack respects transport metadata

library(testthat)

set.seed(2024)

test_that("dkge_infer errors without transport when cluster sizes differ", {
  data <- create_mismatched_data()
  fit <- dkge_fit(data$betas, data$designs, K = data$K, rank = 2)
  expect_error(
    suppressWarnings(dkge_infer(
      fit, c(1, -1, 0), method = "loso", n_perm = 100,
      allow_approximate_alignment = TRUE
    )),
    "Subject cluster counts differ",
    fixed = FALSE
  )
})

test_that("current static fold-loading transport is typed and refused for inference", {
  data <- create_mismatched_data()
  fit <- dkge_fit(data$betas, data$designs, K = data$K, rank = 2)

  transport_cfg <- list(
    centroids = data$centroids,
    medoid = 1L,
    mapper = dkge_mapper_spec("ridge", lambda = 1e-2)
  )

  expect_error(
    suppressWarnings(dkge_infer(
      fit, c(1, -1, 0), method = "loso", n_perm = 100,
      transport = transport_cfg
    )),
    "Transport inside .* is retired|dkge_infer_aligned",
    class = "dkge_alignment_ineligible_error"
  )

  contrast <- suppressWarnings(
    dkge_contrast(fit, c(1, -1, 0), method = "loso")
  )
  transported <- suppressWarnings(dkge_transport_contrasts_to_medoid(
    fit, contrast, medoid = 1L, centroids = data$centroids,
    mapper = transport_cfg$mapper
  ))
  expect_equal(ncol(transported[[1]]$subj_values),
               nrow(data$centroids[[1]]))
  expect_identical(attr(transported, "alignment_mode"), "fold_safe")
  expect_false(attr(transported, "alignment_provenance")$inferentially_eligible)
  expect_s3_class(attr(transported, "aligned_maps"), "dkge_aligned_maps")
  expect_identical(attr(transported, "aligned_maps")$eligibility$status,
                   "ineligible")
})

test_that("inferential transport never falls back to full-fit loadings", {
  data <- create_mismatched_data()
  fit <- dkge_fit(data$betas, data$designs, K = data$K, rank = 2)
  contrast <- suppressWarnings(
    dkge_contrast(fit, c(1, -1, 0), method = "loso")
  )
  contrast$metadata$alignment_receipts <- NULL

  expect_error(
    suppressWarnings(dkge_transport_contrasts_to_medoid(
      fit, contrast, medoid = 1L, centroids = data$centroids,
      mapper = dkge_mapper_spec("ridge", lambda = 1e-2)
    )),
    class = "dkge_alignment_receipt_error"
  )

  descriptive <- suppressWarnings(dkge_transport_contrasts_to_medoid(
    fit, contrast, medoid = 1L, centroids = data$centroids,
    mapper = dkge_mapper_spec("ridge", lambda = 1e-2),
    alignment_mode = "descriptive"
  ))
  expect_false(attr(descriptive, "alignment_provenance")$inferentially_eligible)
})
