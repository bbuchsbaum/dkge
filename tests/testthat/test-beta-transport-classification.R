as.matrix.beta_transport_contrasts <- function(x, contrast, ...) {
  x$values[[contrast]]
}

beta_transport_contrasts <- function(q = 2L) {
  values <- list(target1 = matrix(
    seq_len(6L * q), nrow = 6L, ncol = q,
    dimnames = list(paste0("s", 1:6), paste0("feature", seq_len(q)))
  ))
  structure(
    list(contrasts = list(target1 = 1L), method = "loso", values = values),
    class = c("beta_transport_contrasts", "list")
  )
}

test_that("loading-derived Sinkhorn transport is sign-sensitive", {
  spec <- dkge_mapper_spec(
    "sinkhorn", epsilon = 0.2,
    lambda_emb = 1, lambda_spa = 0,
    warm_start = FALSE
  )
  source_feat <- matrix(c(-1, 2, 4), ncol = 1)
  target_feat <- matrix(c(-2, 0.5, 3), ncol = 1)
  observed <- fit_mapper(spec, source_feat = source_feat,
                         target_feat = target_feat)
  flipped <- fit_mapper(spec, source_feat = -source_feat,
                        target_feat = target_feat)

  expect_gt(max(abs(observed$fit_info$cost - flipped$fit_info$cost)), 0)
  expect_gt(max(abs(observed$operator - flipped$operator)), 0)
})

test_that("inferential transport rejects sign-sensitive in-sample provenance", {
  provenance <- list(
    class = "data_derived_loading",
    training = "inferential_sample",
    sign_invariant = FALSE,
    recomputed_within_randomization = FALSE
  )
  expect_error(
    .dkge_validate_transport_inference_provenance(provenance),
    "sign-sensitive.*inferential sample.*recomputed",
    class = "dkge_transport_inference_error"
  )
})

test_that("public transport provenance classes have explicit exactness", {
  safe_fixed <- c(
    "independent_training", "prespecified_frozen",
    "geometry_only", "sign_invariant"
  )
  for (provenance_class in safe_fixed) {
    provenance <- dkge_transport_provenance(provenance_class)
    expect_s3_class(provenance, "dkge_transport_provenance")
    expect_match(provenance$exactness, "exact")
    expect_silent(
      .dkge_validate_transport_inference_provenance(provenance)
    )
  }

  recomputed <- dkge_transport_provenance("fully_recomputed")
  expect_silent(
    .dkge_validate_transport_inference_provenance(
      recomputed, randomization_recompute = function(...) NULL
    )
  )
  expect_error(
    .dkge_validate_transport_inference_provenance(recomputed),
    "randomization_recompute.*every randomization",
    class = "dkge_transport_inference_error"
  )

  approximate <- dkge_transport_provenance("conditional_approximate")
  expect_error(
    .dkge_validate_transport_inference_provenance(approximate),
    "conditional/approximate.*descriptively",
    class = "dkge_transport_inference_error"
  )
})

test_that("transported inference requires provenance before fitting", {
  expect_error(
    dkge_infer(
      list(), 1, correction = "none", n_perm = 100,
      transport = list(centroids = list(matrix(0, 1, 1)))
    ),
    "sign-sensitive.*inferential sample.*recomputed",
    class = "dkge_transport_inference_error"
  )
})

test_that("fully recomputed transport runs inside every sign randomization", {
  contrasts <- beta_transport_contrasts(q = 2L)
  calls <- 0L
  recompute <- function(signs, target_index, target_name,
                        contrast_results, observed_values) {
    calls <<- calls + 1L
    signs * observed_values
  }
  set.seed(812)
  result <- .infer_signflip(
    contrasts,
    n_perm = 100,
    correction = "maxT",
    mapped_values = contrasts$values,
    transport_recompute = recompute,
    transport_provenance = dkge_transport_provenance("fully_recomputed")
  )
  expect_equal(calls, 100L)
  expect_identical(result$exactness,
                   "randomization_exact_recomputed_operator")
})

test_that("classification permutations require a preselected scalar lambda", {
  expect_error(
    .dkge_validate_classification_selection(
      n_perm = 99L,
      lambda = NULL,
      lambda_grid = NULL,
      lambda_fun = NULL
    ),
    "preselected scalar `lambda`",
    class = "dkge_classification_inference_error"
  )
  expect_error(
    .dkge_validate_classification_selection(
      n_perm = 99L,
      lambda = NULL,
      lambda_grid = c(0.1, 1),
      lambda_fun = NULL
    ),
    "preselected scalar `lambda`",
    class = "dkge_classification_inference_error"
  )
  expect_error(
    .dkge_validate_classification_selection(
      n_perm = 99L,
      lambda = 0.1,
      lambda_grid = NULL,
      lambda_fun = function(...) 0.1
    ),
    "preselected scalar `lambda`",
    class = "dkge_classification_inference_error"
  )
  expect_equal(
    .dkge_validate_classification_selection(
      n_perm = 99L,
      lambda = 0.1,
      lambda_grid = NULL,
      lambda_fun = NULL
    ),
    0.1
  )

  fake_fit <- structure(list(), class = c("dkge", "list"))
  expect_error(
    dkge_classify(
      fake_fit,
      targets = matrix(c(1, -1), nrow = 2),
      n_perm = 99L,
      control = list(lambda_grid = c(0.1, 1))
    ),
    "preselected scalar `lambda`",
    class = "dkge_classification_inference_error"
  )

  fit_shell <- structure(
    list(Btil = list(matrix(0, 2, 1), matrix(0, 2, 1)),
         subject_ids = c("s1", "s2")),
    class = c("dkge", "list")
  )
  expect_error(
    dkge_classify(
      fit_shell,
      targets = rbind(class1 = c(1, 0), class2 = c(0, 1)),
      mode = "cell_cross", lambda = 0.1, n_perm = 9L
    ),
    "randomization_recompute.*complete.*representation",
    class = "dkge_classification_inference_error"
  )
})
