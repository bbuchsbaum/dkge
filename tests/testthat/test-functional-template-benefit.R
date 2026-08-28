library(testthat)

template_benefit_path <- dkge_validation_path(
  "functional-alignment-template", "template-benefit.R"
)
if (!is.na(template_benefit_path)) source(template_benefit_path, local = TRUE)

test_that("iterative template is validated on held-out known-warp maps", {
  skip_on_cran()
  skip_if(is.na(template_benefit_path),
          "functional-template validation sources are unavailable")
  result <- dktf_run(dktf_config(911L))
  expect_true(result$all_gates_pass, info = paste(
    names(result$gates)[!unlist(result$gates)], collapse = ", "
  ))
  expect_setequal(result$metrics$arm, c(
    "iterative_template", "raw_functional_medoid", "geometry_only",
    "ordinary_mni_average"
  ))
  expect_true(all(is.finite(result$metrics$correlation)))
  expect_true(all(is.finite(result$metrics$rmse)))
  expect_true(all(is.finite(result$metrics$latent_error)))
  expect_true(all(is.finite(result$metrics$point_spread)))
  expect_setequal(names(result$aligned_rows), result$metrics$arm)
  expect_true(all(result$solver_converged))
  expect_true(all(vapply(result$aligned_rows, is.matrix, logical(1))))
  expect_true(all(vapply(result$aligned_rows, nrow, integer(1)) ==
                    result$config$S))
  expect_identical(result$template$eligibility$status, "approximate")
  expect_true(result$template$fitting$converged)
  expect_identical(
    result$initializer$provenance$initialization_channel,
    "same_as_alignment_features"
  )
  expect_identical(
    result$initializer$provenance$alignment_features_hash,
    result$template$alignment_features_hash
  )
})
