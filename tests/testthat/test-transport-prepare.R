# test-transport-prepare.R
# Coverage for dkge_prepare_transport helper

library(testthat)

set.seed(99)

test_that("dkge_prepare_transport builds mapper cache without prior operators", {
  fx <- make_small_fit(S = 3, q = 3, P = 5, T = 18, rank = 2, seed = 23)
  fit <- fx$fit
  betas <- fx$betas
  centroids <- replicate(length(betas), matrix(runif(ncol(betas[[1]]) * 3), ncol = 3), simplify = FALSE)
  loadings <- dkge_predict_loadings(fit, betas)

  mapper_spec <- dkge_mapper_spec("ridge", lambda = 1e-3)

  prep <- dkge_prepare_transport(fit,
                                 centroids = centroids,
                                 loadings = loadings,
                                 mapper = mapper_spec,
                                 medoid = 2L)

  expect_s3_class(prep, "dkge_fitted_alignment")
  expect_named(prep$structural_receipt$component_hashes,
               c("features", "centroids", "masses", "reference_subject",
                 "mapper", "preprocessing", "subject_order",
                 "value_semantics", "reference_support",
                 "functional_template", "eligibility",
                 "reference_selection"))
  expect_equal(length(prep$operators), length(loadings))
  expect_equal(prep$medoid, 2L)
  expect_true(is.matrix(prep$feature_ref))
  expect_identical(prep$feature_ref, prep$feature_list[[prep$medoid]])
  expect_true(all(vapply(prep$operators, is.matrix, logical(1))))

  ref_dim <- nrow(prep$feature_ref)
  expect_true(all(vapply(prep$operators, function(op) ncol(op) == ref_dim, logical(1))))
  expect_true(all(vapply(prep$operators, function(op) nrow(op) > 0, logical(1))))

  expect_equal(prep$mapper_spec$strategy, "ridge")
  expect_identical(prep$eligibility$status, "ineligible")
  expect_match(prep$eligibility$reason, "descript", ignore.case = TRUE)
  expect_equal(unname(prep$centroids), centroids)
  expect_identical(names(prep$centroids), fit$subject_ids)

  self_operator <- prep$operators[[prep$reference_subject]]
  expect_false(isTRUE(all.equal(self_operator, diag(1, ref_dim), tolerance = 0)))
  expect_true(prep$diagnostics[[prep$reference_subject]]$self_map)
  expect_s3_class(prep$reference_selection, "dkge_reference_selection")
  expect_identical(prep$reference_method, "explicit")
})

test_that("loose preparation cannot mint inferential provenance", {
  fx <- make_small_fit(S = 3, q = 3, P = 5, T = 18, rank = 2, seed = 25)
  centroids <- replicate(3, matrix(runif(15), ncol = 3), simplify = FALSE)
  loadings <- dkge_predict_loadings(fx$fit, fx$betas)
  spoof <- list(
    source = "independent_spoof",
    feature_source = "independent",
    estimator_source = "independent",
    feature_provenance = list(independent_data_hash = "not-evidence")
  )
  expect_error(
    dkge_prepare_transport(
      fx$fit, centroids = centroids, loadings = loadings,
      mapper = dkge_mapper_spec("ridge"), preprocessing = spoof
    ),
    "typed.*alignment_features|caller-authored|preprocessing",
    class = "dkge_alignment_preprocessing_error"
  )
})

test_that("fitted alignment rejects every mutated structural input", {
  fx <- make_small_fit(S = 3, q = 3, P = 5, T = 18, rank = 2, seed = 24)
  fit <- fx$fit
  centroids <- replicate(3, matrix(runif(15), ncol = 3), simplify = FALSE)
  loadings <- dkge_predict_loadings(fit, fx$betas)
  sizes <- lapply(loadings, function(A) seq_len(nrow(A)))
  mapper <- dkge_mapper_spec("ridge", lambda = 1e-3,
                             value_type = "intensive")
  aligned <- dkge_prepare_transport(
    fit, centroids = centroids, loadings = loadings, sizes = sizes,
    mapper = mapper, medoid = 2L
  )
  values <- lapply(loadings, function(A) A[, 1])

  miss <- dkge:::.dkge_transport_to_medoid(
    mapper, values, loadings, centroids, sizes, 2L,
    subject_ids = fit$subject_ids,
    preprocessing = aligned$preprocessing
  )
  expect_false(miss$cache_provenance$hit)
  expect_s3_class(miss$fitted_alignment, "dkge_fitted_alignment")

  hit <- dkge:::.dkge_transport_to_medoid(
    mapper, values, loadings, centroids, sizes, 2L,
    transport_cache = aligned,
    subject_ids = fit$subject_ids,
    preprocessing = aligned$preprocessing
  )
  expect_true(hit$cache_provenance$hit)
  expect_identical(hit$cache_provenance$validation, "exact_structural_match")

  mutations <- list(
    features = function() {
      x <- loadings; x[[1]][1, 1] <- x[[1]][1, 1] + 0.1; x
    },
    centroids = function() {
      x <- centroids; x[[1]][1, 1] <- x[[1]][1, 1] + 0.1; x
    },
    masses = function() {
      x <- sizes; x[[1]][1] <- x[[1]][1] + 1; x
    }
  )
  for (component in names(mutations)) {
    changed_loadings <- if (component == "features") mutations[[component]]() else loadings
    changed_centroids <- if (component == "centroids") mutations[[component]]() else centroids
    changed_sizes <- if (component == "masses") mutations[[component]]() else sizes
    expected_class <- if (component %in% c("centroids", "masses")) {
      "dkge_reference_selection_error"
    } else {
      "dkge_alignment_cache_mismatch"
    }
    expect_error(
      dkge:::.dkge_transport_to_medoid(
        mapper, values, changed_loadings, changed_centroids, changed_sizes, 2L,
        transport_cache = aligned,
        subject_ids = fit$subject_ids,
        preprocessing = aligned$preprocessing
      ),
      component,
      class = expected_class
    )
  }

  mapper_changed <- dkge_mapper_spec("ridge", lambda = 2e-3,
                                     value_type = "intensive")
  expect_error(
    dkge:::.dkge_transport_to_medoid(
      mapper_changed, values, loadings, centroids, sizes, 2L,
      transport_cache = aligned, subject_ids = fit$subject_ids,
      preprocessing = aligned$preprocessing
    ),
    "mapper", class = "dkge_alignment_cache_mismatch"
  )
  semantics_changed <- dkge_mapper_spec("ridge", lambda = 1e-3,
                                        value_type = "extensive")
  expect_error(
    dkge:::.dkge_transport_to_medoid(
      semantics_changed, values, loadings, centroids, sizes, 2L,
      transport_cache = aligned, subject_ids = fit$subject_ids,
      preprocessing = aligned$preprocessing
    ),
    "value_semantics", class = "dkge_alignment_cache_mismatch"
  )
  expect_error(
    dkge:::.dkge_transport_to_medoid(
      mapper, values, loadings, centroids, sizes, 1L,
      transport_cache = aligned, subject_ids = fit$subject_ids,
      preprocessing = aligned$preprocessing
    ),
    "conflicts", class = "dkge_reference_selection_error"
  )
  expect_error(
    dkge:::.dkge_transport_to_medoid(
      mapper, values, loadings, centroids, sizes, 2L,
      transport_cache = aligned, subject_ids = rev(fit$subject_ids),
      preprocessing = aligned$preprocessing
    ),
    "subject order", class = "dkge_reference_selection_error"
  )
  expect_error(
    dkge:::.dkge_transport_to_medoid(
      mapper, values, loadings, centroids, sizes, 2L,
      transport_cache = aligned, subject_ids = fit$subject_ids,
      preprocessing = c(aligned$preprocessing, list(mutated = TRUE))
    ),
    "preprocessing", class = "dkge_alignment_cache_mismatch"
  )
})
