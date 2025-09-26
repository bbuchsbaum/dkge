# test-neuroim2.R
# Light-touch tests for neuroim2 integration helpers

library(testthat)

test_that("cluster aggregation works with neuroim2 objects", {
  skip_if_not_installed("neuroim2")
  skip_if_not_installed("fmrireg")

  time_points <- 4L
  spatial_dims <- c(3L, 1L, 1L)

  # Construct deterministic voxel time series (time x voxel)
  voxel_ts <- matrix(
    c(1, 3, 5,
      2, 4, 6,
      3, 5, 7,
      4, 6, 8),
    nrow = time_points,
    byrow = TRUE
  )

  # Create NeuroVec time-series (TxXxYxZ) with time in the final dimension
  ts_array <- array(0, dim = c(spatial_dims, time_points))
  for (t in seq_len(time_points)) {
    ts_array[, , , t] <- array(voxel_ts[t, ], dim = spatial_dims)
  }
  space4d <- neuroim2::NeuroSpace(
    dim = c(spatial_dims, time_points),
    spacing = rep(1, 3),
    origin = rep(0, 3)
  )
  bv <- neuroim2::NeuroVec(ts_array, space4d)

  space3d <- neuroim2::NeuroSpace(
    dim = spatial_dims,
    spacing = rep(1, 3),
    origin = rep(0, 3)
  )

  # Define cluster labels (first two voxels -> cluster 1, third voxel -> cluster 2)
  label_array <- array(c(1, 1, 2), dim = spatial_dims)
  labels <- neuroim2::NeuroVol(label_array, space3d)

  # Expected cluster sums matching dkge_cluster_ts aggregation logic
  ids <- sort(unique(label_array[label_array > 0]))
  keep <- as.vector(label_array) > 0
  y_mat <- t(neuroim2::as.matrix(bv))  # time x voxels
  label_sel <- as.vector(label_array)[keep]
  y_k <- y_mat[, keep, drop = FALSE]
  expected_ts <- sapply(ids, function(id) {
    rowSums(y_k[, label_sel == id, drop = FALSE])
  })
  if (!is.matrix(expected_ts)) {
    expected_ts <- matrix(expected_ts, ncol = length(ids))
  }
  colnames(expected_ts) <- as.character(ids)

  agg_ts <- dkge_cluster_ts(bv, labels)
  expect_equal(agg_ts, expected_ts)

  # Compare dkge_cluster_betas against closed-form OLS solution
  design <- cbind(1, seq_len(time_points))
  betas_expected <- tryCatch(
    fmrireg::fmri_ols_fit(design, expected_ts)$beta,
    error = function(e) {
      if (grepl("lazy-load database", e$message, fixed = TRUE)) {
        skip("fmrireg lazy-load database unavailable in this environment.")
      }
      stop(e)
    }
  )
  betas_dkge <- tryCatch(
    dkge_cluster_betas(bv, design, labels),
    error = function(e) {
      if (grepl("lazy-load database", e$message, fixed = TRUE)) {
        skip("fmrireg lazy-load database unavailable in this environment.")
      }
      stop(e)
    }
  )
  expect_equal(betas_dkge, betas_expected, tolerance = 1e-8)
})

test_that("dkge_write_group_map uses named group values", {
  skip_if_not_installed("neuroim2")
  ns <- asNamespace("neuroim2")
  has_brain_volume <- exists("BrainVolume", envir = ns, inherits = FALSE)
  has_neurovol <- exists("NeuroVol", envir = ns, inherits = FALSE)
  if (!has_brain_volume && !has_neurovol) {
    skip("neuroim2 BrainVolume/NeuroVol constructors not available in this environment")
  }

  spatial_dims <- c(2L, 1L, 1L)
  space3d <- neuroim2::NeuroSpace(dim = spatial_dims,
                                  spacing = rep(1, 3),
                                  origin = rep(0, 3))
  label_array <- array(c(1, 2), dim = spatial_dims)
  medoid_labels <- neuroim2::NeuroVol(label_array, space3d)

  group_values <- c(`2` = 42, `1` = 7)
  vol <- dkge_write_group_map(group_values, medoid_labels, out_file = NULL)
  expect_true(inherits(vol, c("BrainVolume", "NeuroVol")))
  expect_equal(neuroim2::values(vol), array(c(7, 42), dim = spatial_dims))

  dup_values <- c(`1` = 1, `1` = 2)
  expect_error(dkge_write_group_map(dup_values, medoid_labels, out_file = NULL), "names must be unique")
})
