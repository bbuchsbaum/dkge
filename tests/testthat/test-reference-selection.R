library(testthat)

make_reference_selection_fixture <- function(seed = 260840L) {
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
  independent_betas <- function(sd, offset) {
    set.seed(seed + offset)
    out <- lapply(latent, function(signal) {
      signal + matrix(
        rnorm(length(signal), sd = sd), nrow(signal), ncol(signal),
        dimnames = dimnames(signal)
      )
    })
    names(out) <- fit$subject_ids
    out
  }
  training <- dkge_alignment_features(
    fit, contrast,
    independent_betas = independent_betas(0.25, 1L),
    independent_data_hash = paste0("reference-training-", seed)
  )
  validation <- dkge_alignment_features(
    fit, contrast,
    independent_betas = independent_betas(0.4, 2L),
    independent_data_hash = paste0("reference-validation-", seed)
  )
  centroids <- lapply(seq_len(S), function(s) {
    cbind(x = seq_len(P), y = sin(seq_len(P) / 2) + s / 20, z = 0)
  })
  names(centroids) <- fit$subject_ids
  sizes <- stats::setNames(lapply(seq_len(S), function(s) {
    seq(1, 2, length.out = P)
  }), fit$subject_ids)
  mapper <- dkge_mapper_spec(
    "sinkhorn", epsilon = 0.1, lambda_emb = 0.8, lambda_spa = 0.2,
    sigma_mm = 15, warm_start = FALSE
  )
  list(
    fit = fit, contrast = contrast, training = training,
    validation = validation, centroids = centroids, sizes = sizes,
    mapper = mapper
  )
}

test_that("explicit references are fixed inputs, not medoids", {
  centroids <- list(
    z = cbind(0:2, 0, 0),
    a = cbind(0:2, 0, 0),
    m = cbind(0:2, 0, 0)
  )
  selected <- dkge_select_reference_subject(
    centroids, method = "explicit", reference_subject = "m"
  )
  expect_s3_class(selected, "dkge_reference_selection")
  expect_identical(selected$reference_subject_id, "m")
  expect_identical(selected$reference_subject, 3L)
  expect_null(selected$medoid)
  expect_identical(selected$eligibility$status, "eligible")

  expect_error(
    dkge_select_reference_subject(
      centroids, method = "explicit", reference_subject = 1.9
    ),
    "valid subject index",
    class = "dkge_reference_selection_error"
  )
  expect_identical(
    dkge_select_reference_subject(
      centroids, method = "explicit", reference_subject = 2
    )$reference_subject,
    2L
  )

  mutated <- selected
  mutated$reference_subject <- 1L
  expect_error(
    print(mutated), "mutated", class = "dkge_reference_selection_error"
  )
})

test_that("geometry-only medoids are order invariant with deterministic ties", {
  centroids <- list(
    z = cbind(0:2, 0, 0),
    a = cbind(0:2, 0, 0),
    m = cbind(0:2, 0, 0)
  )
  base <- dkge_select_reference_subject(
    centroids, method = "geometry_only",
    mapper = dkge_mapper_spec("sinkhorn", epsilon = 0.2,
                              lambda_spa = 1, warm_start = FALSE)
  )
  permutation <- c("m", "z", "a")
  reordered <- dkge_select_reference_subject(
    centroids[permutation], method = "geometry_only",
    mapper = dkge_mapper_spec("sinkhorn", epsilon = 0.2,
                              lambda_spa = 1, warm_start = FALSE)
  )

  expect_identical(base$reference_subject_id, "a")
  expect_identical(reordered$reference_subject_id, "a")
  expect_identical(base$method, "geometry_only")
  expect_identical(base$mapper_spec$params$lambda_emb, 0)
  expect_lt(diff(range(base$scores)), 1e-12)
  expect_false("sinkhorn_divergence" %in% names(base))
})

test_that("reference receipts cannot be reused under different mapper scaling", {
  fx <- make_reference_selection_fixture(seed = 260839L)
  selected <- dkge_select_reference_subject(
    fx$centroids, sizes = fx$sizes, mapper = fx$mapper,
    method = "geometry_only", subject_ids = fx$fit$subject_ids
  )
  changed <- fx$mapper
  changed$params$epsilon <- changed$params$epsilon * 2
  expect_error(
    dkge_prepare_alignment(
      fx$fit, fx$training, fx$centroids, sizes = fx$sizes,
      mapper = changed, reference_selection = selected
    ),
    "different mapper scaling",
    class = "dkge_reference_selection_error"
  )
})

test_that("functional selection fits and scores on distinct independent channels", {
  fx <- make_reference_selection_fixture()
  selected <- dkge_select_reference_subject(
    fx$centroids,
    alignment_features = fx$training,
    validation_features = fx$validation,
    sizes = fx$sizes,
    mapper = fx$mapper,
    method = "functional_heldout"
  )

  expect_identical(selected$method, "functional_heldout")
  expect_identical(selected$criterion,
                   "symmetric_normalized_reconstruction_loss")
  expect_identical(selected$eligibility$status, "eligible")
  expect_identical(selected$reference_subject_id,
                   names(which.min(selected$scores)))
  expect_identical(selected$alignment_features_hash,
                   fx$training$structural_hash)
  expect_identical(selected$validation_features_hash,
                   fx$validation$structural_hash)
  expect_false(selected$fixed_scaling$candidate_specific_rescaling)
  expect_false(selected$fixed_scaling$epsilon_calibration)
  expect_equal(selected$pairwise_loss, t(selected$pairwise_loss))

  expect_error(
    dkge_select_reference_subject(
      fx$centroids,
      alignment_features = fx$training,
      validation_features = fx$training,
      sizes = fx$sizes,
      mapper = fx$mapper,
      method = "functional_heldout"
    ),
    "not independently identified",
    class = "dkge_reference_selection_error"
  )
  expect_error(
    dkge_select_reference_subject(
      fx$centroids,
      alignment_features = fx$training,
      sizes = fx$sizes,
      mapper = fx$mapper,
      method = "functional_heldout"
    ),
    "second, held-out",
    class = "dkge_reference_selection_error"
  )
})

test_that("known functional correspondence can overturn misleading geometry", {
  set.seed(9)
  ids <- c("left", "bridge", "right")
  q <- 3L
  effects <- paste0("e", seq_len(q))
  latent <- list(
    c(-1, -0.4, 0.2),
    seq(-1, 1, length.out = 5),
    c(-0.2, 0.4, 1)
  )
  subjects <- lapply(seq_along(ids), function(s) {
    P <- length(latent[[s]])
    B <- matrix(rnorm(q * P), q, P,
                dimnames = list(effects, paste0("p", seq_len(P))))
    X <- matrix(rnorm(12 * q), 12, q,
                dimnames = list(NULL, effects))
    dkge_subject(B, X, id = ids[[s]])
  })
  K <- diag(q)
  dimnames(K) <- list(effects, effects)
  fit <- suppressWarnings(dkge_fit(
    dkge_data(subjects), K = K, rank = 2,
    w_method = "none", effect_scaling = "none"
  ))
  contrast <- suppressWarnings(dkge_contrast(
    fit, c(1, -1, 0), method = "loso", align = FALSE
  ))
  signature <- function(u, delta) {
    cbind(1 + delta, u + delta, sin(pi * u / 2) + delta * u^3)
  }
  beta_channel <- function(delta) {
    out <- lapply(latent, function(u) {
      B <- t(signature(u, delta))
      dimnames(B) <- list(effects, paste0("p", seq_along(u)))
      B
    })
    names(out) <- ids
    out
  }
  training <- dkge_alignment_features(
    fit, contrast,
    independent_betas = beta_channel(0.03),
    independent_data_hash = "known-correspondence-training"
  )
  validation <- dkge_alignment_features(
    fit, contrast,
    independent_betas = beta_channel(-0.04),
    independent_data_hash = "known-correspondence-validation"
  )
  # Geometry reverses the two sparse supports, while functional signatures
  # identify the denser bridge support shared by both held-out channels.
  centroids <- list(
    left = cbind(rev(seq_along(latent[[1]])), 0, 0),
    bridge = cbind(seq_along(latent[[2]]), 0, 0),
    right = cbind(rev(seq_along(latent[[3]])), 0, 0)
  )
  functional <- dkge_select_reference_subject(
    centroids, training, validation,
    mapper = dkge_mapper_spec(
      "sinkhorn", epsilon = 0.03, lambda_emb = 1, lambda_spa = 0,
      warm_start = FALSE
    ),
    method = "functional_heldout"
  )
  geometry <- dkge_select_reference_subject(
    centroids,
    mapper = dkge_mapper_spec(
      "sinkhorn", epsilon = 0.03, lambda_spa = 1, warm_start = FALSE
    ),
    method = "geometry_only"
  )
  explicit <- dkge_select_reference_subject(
    centroids, method = "explicit", reference_subject = "right"
  )

  expect_identical(functional$reference_subject_id, "bridge")
  expect_identical(geometry$reference_subject_id, "left")
  expect_identical(explicit$reference_subject_id, "right")
  expect_true(functional$scores[["bridge"]] <
                min(functional$scores[c("left", "right")]))
  expect_true(geometry$scores[["left"]] < geometry$scores[["bridge"]])
  expect_true(functional$eligibility$eligible)
  expect_true(geometry$eligibility$eligible)
  expect_null(explicit$medoid)
})

test_that("training-feature selection is labelled descriptive and ineligible", {
  fx <- make_reference_selection_fixture(seed = 260841L)
  selected <- dkge_select_reference_subject(
    fx$centroids,
    alignment_features = fx$training,
    sizes = fx$sizes,
    mapper = fx$mapper,
    method = "descriptive_training"
  )
  expect_identical(selected$method, "descriptive_training")
  expect_identical(selected$eligibility$status, "ineligible")
  expect_false(selected$eligibility$eligible)
  expect_identical(selected$alignment_features_hash,
                   selected$validation_features_hash)
})

test_that("reference selection is structurally bound to fitted transport", {
  fx <- make_reference_selection_fixture(seed = 260842L)
  selected <- dkge_select_reference_subject(
    fx$centroids,
    alignment_features = fx$training,
    validation_features = fx$validation,
    sizes = fx$sizes,
    mapper = fx$mapper,
    method = "functional_heldout"
  )
  alignment <- dkge_prepare_transport(
    fx$fit,
    centroids = fx$centroids,
    sizes = fx$sizes,
    mapper = fx$mapper,
    reference_selection = selected,
    alignment_features = fx$training
  )

  expect_identical(alignment$reference_subject,
                   selected$reference_subject)
  expect_identical(alignment$reference_selection$structural_hash,
                   selected$structural_hash)
  expect_identical(alignment$reference_method, "functional_heldout")
  expect_true(alignment$reference_is_medoid)

  expect_error(
    dkge_prepare_transport(
      fx$fit,
      centroids = fx$centroids,
      sizes = fx$sizes,
      mapper = fx$mapper,
      medoid = if (selected$reference_subject == 1L) 2L else 1L,
      reference_selection = selected,
      alignment_features = fx$training
    ),
    "conflicts",
    class = "dkge_reference_selection_error"
  )
  expect_error(
    dkge_prepare_transport(
      fx$fit,
      centroids = fx$centroids,
      sizes = fx$sizes,
      mapper = dkge_mapper_spec(
        "sinkhorn", epsilon = 0.2, lambda_emb = 0.8, lambda_spa = 0.2
      ),
      reference_selection = selected,
      alignment_features = fx$training
    ),
    "different mapper scaling",
    class = "dkge_reference_selection_error"
  )
})

test_that("diffuse plans expose their point spread instead of masquerading as identity", {
  centroids <- list(
    a = cbind(0:3, 0, 0),
    b = cbind(0:3, 0, 0),
    c = cbind(0:3, 0, 0)
  )
  selected <- dkge_select_reference_subject(
    centroids,
    method = "geometry_only",
    mapper = dkge_mapper_spec(
      "sinkhorn", epsilon = 100, lambda_spa = 1, warm_start = FALSE
    )
  )
  spread <- selected$diagnostics$mean_effective_points
  expect_true(all(spread[!is.na(spread)] > 3.9))
  expect_true(all(selected$diagnostics$converged[!is.na(
    selected$diagnostics$converged
  )]))
  expect_true(all(selected$diagnostics$marginal_error[!is.na(
    selected$diagnostics$marginal_error
  )] <= 1e-4))

  calibrated <- dkge_mapper_spec(
    "sinkhorn", epsilon = 0.1, lambda_spa = 1,
    epsilon_calibration = list(
      target_effective_points = 2,
      epsilon_grid = c(0.01, 0.1),
      calibration_data_hash = "heldout-calibration"
    )
  )
  expect_error(
    dkge_select_reference_subject(
      centroids, method = "geometry_only", mapper = calibrated
    ),
    "fixed mapper scaling",
    class = "dkge_reference_selection_error"
  )
})

test_that("reference scoring fails closed on any non-converged Sinkhorn pair", {
  fx <- make_reference_selection_fixture(seed = 260843L)
  failed <- dkge_mapper_spec(
    "sinkhorn", epsilon = 1e-4, lambda_emb = 1, lambda_spa = 0.2,
    max_iter = 1L, tol = 1e-14, warm_start = FALSE
  )
  expect_error(
    suppressWarnings(dkge_select_reference_subject(
      fx$centroids,
      alignment_features = fx$training,
      validation_features = fx$validation,
      sizes = fx$sizes,
      mapper = failed,
      method = "functional_heldout"
    )),
    "did not converge|numerical",
    class = "dkge_reference_numerical_error"
  )
})
