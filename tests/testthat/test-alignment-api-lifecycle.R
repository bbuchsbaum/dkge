library(testthat)

make_alignment_api_fixture <- local({
  cached <- NULL
  function(seed = 260860L) {
    if (!is.null(cached) && identical(cached$seed, seed)) return(cached)
    set.seed(seed)
    S <- 4L
    q <- 5L
    P <- 7L
    effects <- paste0("e", seq_len(q))
    parcels <- paste0("p", seq_len(P))
    K <- diag(seq(q, 1), q)
    dimnames(K) <- list(effects, effects)
    latent <- lapply(seq_len(S), function(s) {
      matrix(rnorm(q * P), q, P, dimnames = list(effects, parcels))
    })
    subjects <- lapply(seq_len(S), function(s) {
      B <- latent[[s]] + matrix(
        rnorm(q * P, sd = 0.5), q, P,
        dimnames = list(effects, parcels)
      )
      X <- matrix(rnorm((q + 6L) * q), q + 6L, q,
                  dimnames = list(NULL, effects))
      dkge_subject(B, X, id = paste0("s", s))
    })
    fit <- suppressWarnings(dkge_fit(
      dkge_data(subjects), K = K, rank = 2,
      w_method = "none", effect_scaling = "pooled_design"
    ))
    contrast <- suppressWarnings(dkge_contrast(
      fit, c(1, -1, rep(0, q - 2L)), method = "loso", align = FALSE
    ))
    make_channel <- function(offset, data_hash) {
      set.seed(seed + offset)
      betas <- lapply(latent, function(signal) {
        signal + matrix(
          rnorm(length(signal), sd = 0.3), nrow(signal), ncol(signal),
          dimnames = dimnames(signal)
        )
      })
      names(betas) <- fit$subject_ids
      dkge_alignment_features(
        fit, contrast,
        independent_betas = betas,
        independent_data_hash = data_hash
      )
    }
    training <- make_channel(1L, paste0("alignment-api-training-", seed))
    validation <- make_channel(2L, paste0("alignment-api-validation-", seed))
    centroids <- lapply(seq_len(S), function(s) {
      cbind(x = seq_len(P), y = s / 20, z = 0)
    })
    names(centroids) <- fit$subject_ids
    sizes <- stats::setNames(lapply(seq_len(S), function(s) {
      seq(1, 2, length.out = P)
    }), fit$subject_ids)
    mapper <- dkge_mapper_spec(
      "sinkhorn", epsilon = 0.1, lambda_emb = 0.8, lambda_spa = 0.2,
      warm_start = FALSE
    )
    cached <<- list(
      seed = seed, fit = fit, contrast = contrast,
      training = training, validation = validation,
      centroids = centroids, sizes = sizes, mapper = mapper
    )
    cached
  }
})

test_that("strict preparation chooses an auditable reference channel", {
  fx <- make_alignment_api_fixture()
  geometry <- dkge_prepare_alignment(
    fx$fit, fx$training, fx$centroids,
    sizes = fx$sizes, mapper = fx$mapper
  )
  expect_s3_class(geometry, "dkge_fitted_alignment")
  expect_identical(geometry$reference_selection$method, "geometry_only")
  expect_true(geometry$reference_selection$eligibility$eligible)

  functional <- dkge_prepare_alignment(
    fx$fit, fx$training, fx$centroids,
    sizes = fx$sizes, mapper = fx$mapper,
    validation_features = fx$validation
  )
  expect_identical(functional$reference_selection$method,
                   "functional_heldout")
  expect_true(functional$reference_is_medoid)

  explicit <- dkge_prepare_alignment(
    fx$fit, fx$training, fx$centroids,
    sizes = fx$sizes, mapper = fx$mapper,
    reference_subject = "s3"
  )
  expect_identical(explicit$reference_selection$method, "explicit")
  expect_identical(explicit$reference_selection$reference_subject_id, "s3")
  expect_false(explicit$reference_is_medoid)
})

test_that("reference contrast workflow requires typed features and reuses exact cache", {
  fx <- make_alignment_api_fixture()
  alignment <- dkge_prepare_alignment(
    fx$fit, fx$training, fx$centroids,
    sizes = fx$sizes, mapper = fx$mapper,
    validation_features = fx$validation
  )
  transported <- dkge_transport_contrasts_to_reference(
    fx$fit, fx$contrast, fx$training,
    centroids = fx$centroids,
    sizes = fx$sizes,
    transport_cache = alignment
  )
  expect_named(transported, names(fx$contrast$values))
  expect_s3_class(attr(transported, "fitted_alignment"),
                  "dkge_fitted_alignment")
  expect_s3_class(attr(transported, "aligned_maps"), "dkge_aligned_maps")
  expect_identical(
    attr(transported, "reference_selection")$structural_hash,
    alignment$reference_selection$structural_hash
  )
  expect_true(transported[[1]]$cache_provenance$hit)

  expect_error(
    dkge_transport_contrasts_to_reference(
      fx$fit, fx$contrast,
      centroids = fx$centroids,
      sizes = fx$sizes,
      mapper = fx$mapper
    ),
    "alignment_features"
  )
})

test_that("descriptive reference scoring is opt-in and ineligible", {
  fx <- make_alignment_api_fixture()
  alignment <- dkge_prepare_alignment(
    fx$fit, fx$training, fx$centroids,
    sizes = fx$sizes, mapper = fx$mapper,
    selection_method = "descriptive_training"
  )
  expect_identical(alignment$reference_selection$method,
                   "descriptive_training")
  expect_identical(alignment$reference_selection$eligibility$status,
                   "ineligible")
  expect_identical(alignment$eligibility$status, "ineligible")
})

test_that("failed Sinkhorn preparation and calibration fail closed", {
  fx <- make_alignment_api_fixture(seed = 260861L)
  failed_mapper <- dkge_mapper_spec(
    "sinkhorn", epsilon = 1e-4, lambda_emb = 1, lambda_spa = 0.2,
    max_iter = 1L, tol = 1e-14, warm_start = FALSE
  )
  expect_error(
    suppressWarnings(dkge_prepare_alignment(
      fx$fit, fx$training, fx$centroids,
      sizes = fx$sizes, mapper = failed_mapper,
      reference_subject = "s2"
    )),
    "numerically invalid|converge",
    class = "dkge_alignment_numerical_error"
  )

  calibrated <- failed_mapper
  calibrated$params$epsilon_calibration <- list(
    target_effective_points = 2,
    epsilon_grid = c(1e-5, 1e-4),
    calibration_data_hash = "separate-calibration"
  )
  expect_error(
    suppressWarnings(dkge_prepare_alignment(
      fx$fit, fx$training, fx$centroids,
      sizes = fx$sizes, mapper = calibrated,
      reference_subject = "s2"
    )),
    "calibration.*converge|candidate.*converge|numerical",
    class = "dkge_alignment_calibration_error"
  )
})

test_that("legacy medoid contrast API warns and preserves its safe result", {
  fx <- make_alignment_api_fixture()
  modern <- dkge_transport_contrasts_to_reference(
    fx$fit, fx$contrast, fx$training,
    centroids = fx$centroids,
    sizes = fx$sizes,
    mapper = fx$mapper,
    reference_subject = "s2"
  )
  expect_warning(
    legacy <- dkge_transport_contrasts_to_medoid(
      fx$fit, fx$contrast,
      medoid = 2L,
      centroids = fx$centroids,
      sizes = fx$sizes,
      mapper = fx$mapper,
      alignment_features = fx$training
    ),
    "deprecated",
    class = "deprecatedWarning"
  )
  expect_equal(legacy[[1]]$subj_values, modern[[1]]$subj_values,
               tolerance = 0)
  expect_identical(attr(legacy, "alignment_mode"), "independent")
})

test_that("legacy paint API is a deprecation-compatible shim", {
  skip_if_not_installed("neuroim2")
  spatial_dims <- c(2L, 1L, 1L)
  space3d <- neuroim2::NeuroSpace(
    dim = spatial_dims, spacing = rep(1, 3), origin = rep(0, 3)
  )
  labels <- neuroim2::NeuroVol(
    array(c(1, 2), dim = spatial_dims), space3d
  )
  modern <- dkge_paint_reference_map(c(`1` = 7, `2` = 42), labels)
  expect_warning(
    legacy <- dkge_paint_medoid_map(c(`1` = 7, `2` = 42), labels),
    "deprecated",
    class = "deprecatedWarning"
  )
  expect_equal(neuroim2::values(legacy), neuroim2::values(modern))
})
