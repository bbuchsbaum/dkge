test_that("formal public-beta full-pipeline calibration evidence is complete", {
  artifact <- function(name) {
    path <- system.file("extdata", name, package = "dkge")
    expect_true(nzchar(path), info = name)
    path
  }

  replicate_path <- artifact(
    "dkge-beta-full-pipeline-null-replicates.csv.gz"
  )
  summary_path <- artifact("dkge-beta-full-pipeline-null-summary.csv")
  metadata_path <- artifact("dkge-beta-full-pipeline-null-metadata.csv")

  replicates <- utils::read.csv(gzfile(replicate_path))
  summary <- utils::read.csv(summary_path)
  metadata_table <- utils::read.csv(metadata_path, stringsAsFactors = FALSE)
  metadata <- stats::setNames(metadata_table$value, metadata_table$key)

  expect_identical(replicates$replicate, seq_len(200L))
  expect_equal(nrow(summary), 3L)
  expect_setequal(
    summary$arm,
    c("transport_maxT", "classification_cell_cross", "bonferroni_family")
  )
  expect_true(all(summary$rate_gate_pass))
  expect_true(all(summary$wilson_contains_005))
  expect_true(all(summary$gate_pass))
  expect_equal(summary$rejections_005, c(5L, 9L, 5L))
  expect_equal(summary$rejection_005, c(0.025, 0.045, 0.025))

  expect_identical(unname(metadata[["status"]]), "formal")
  expect_match(unname(metadata[["source_git_sha"]]), "^[0-9a-f]{40}$")
  expect_identical(unname(metadata[["n_replicates"]]), "200")
  expect_identical(unname(metadata[["n_randomizations"]]), "199")
  expect_identical(unname(metadata[["contract_pass"]]), "TRUE")
  expect_identical(
    unname(metadata[["all_statistical_gates_pass"]]), "TRUE"
  )

  expect_true(all(replicates$kernel_rank == 2L))
  expect_true(all(replicates$effective_rank == 2L))
  expect_true(all(replicates$minimum_fold_rank == 2L))
  expect_true(all(replicates$alignment_present))
  expect_true(all(
    replicates$transport_provenance == "geometry_only"
  ))
  expect_true(all(
    replicates$transport_exactness ==
      "randomization_exact_sign_invariant_operator"
  ))
  expect_true(all(
    replicates$classification_claim_scope ==
      "prospective_heldout_subject"
  ))
  expect_true(all(
    replicates$classification_representation_scope ==
      "basis_cross_fitted_without_heldout_subject"
  ))
  expect_true(all(replicates$lambda == 0.001))
  expect_true(all(is.finite(replicates$p_transport)))
  expect_true(all(is.finite(replicates$p_classification)))

  expect_identical(
    digest::digest(file = replicate_path, algo = "sha256", serialize = FALSE),
    unname(metadata[["replicate_artifact_sha256"]])
  )
  expect_identical(
    digest::digest(file = summary_path, algo = "sha256", serialize = FALSE),
    unname(metadata[["summary_artifact_sha256"]])
  )
})
