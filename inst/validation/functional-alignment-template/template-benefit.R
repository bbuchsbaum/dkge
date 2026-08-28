#!/usr/bin/env Rscript

# Known-correspondence efficacy study for iterative functional templates.
#
# Parcel order encodes the latent functional coordinate u. Each subject's
# anatomical x coordinate is cyclically warped with cohort-balanced shifts, so
# geometry and ordinary coordinate-space averaging are knowingly wrong.
# Functional correspondence
# comes from a separate beta-map replicate. A held-out replicate supplies the
# field being predicted. The display support is subject 1; its functional
# features are deliberately noisier than the rest, exposing the raw-medoid
# variance cost that a group template is meant to address.

dktf_config <- function(seed = 911L) {
  list(
    seed = as.integer(seed),
    S = 8L,
    P = 10L,
    q = 5L,
    reference_subject = 1L,
    reference_feature_noise = 1.5,
    other_feature_noise = 0.35,
    heldout_signal_scale = 1,
    heldout_value_noise = 0.2,
    epsilon = 0.015,
    lambda_emb = 0.95,
    lambda_spa = 0.05,
    sigma_mm = 0.2,
    template_max_iter = 25L,
    template_tolerance = 0.002,
    template_objective_tolerance = 0.002,
    update_rate = 0.5
  )
}

dktf_generate <- function(config) {
  set.seed(config$seed)
  S <- config$S
  P <- config$P
  q <- config$q
  ids <- paste0("subject", seq_len(S))
  effects <- paste0("effect", seq_len(q))
  parcels <- paste0("parcel", seq_len(P))
  latent <- seq(0, 1, length.out = P)
  signature <- cbind(
    sin(2 * pi * latent),
    cos(2 * pi * latent),
    sin(4 * pi * latent)
  )
  primary <- alignment <- vector("list", S)
  heldout <- centroids <- vector("list", S)
  for (s in seq_len(S)) {
    primary[[s]] <- matrix(
      rnorm(q * P), q, P,
      dimnames = list(effects, parcels)
    )
    alignment[[s]] <- matrix(
      rnorm(q * P, sd = 0.1), q, P,
      dimnames = list(effects, parcels)
    )
    feature_noise <- if (s == config$reference_subject) {
      config$reference_feature_noise
    } else {
      config$other_feature_noise
    }
    alignment[[s]][2:4, ] <- t(signature) +
      matrix(rnorm(3L * P, sd = feature_noise), 3L, P)
    heldout[[s]] <- config$heldout_signal_scale * sin(2 * pi * latent) +
      rnorm(P, sd = config$heldout_value_noise)
    # Balanced cyclic anatomical warps make geometry ancillary but genuinely
    # non-corresponding. Unlike independent random permutations, their finite-
    # cohort average cannot accidentally correlate strongly with the truth.
    permutation <- ((seq_len(P) - 1L + (s - 1L)) %% P) + 1L
    centroids[[s]] <- cbind(x = latent[permutation], y = 0, z = 0)
  }
  names(primary) <- names(alignment) <- names(heldout) <-
    names(centroids) <- ids
  design <- diag(q)
  colnames(design) <- effects
  K <- diag(q)
  diag(K)[c(1, 5)] <- 0
  dimnames(K) <- list(effects, effects)
  list(
    ids = ids,
    effects = effects,
    parcels = parcels,
    latent = latent,
    truth = config$heldout_signal_scale * sin(2 * pi * latent),
    primary = primary,
    alignment = alignment,
    heldout = heldout,
    centroids = centroids,
    designs = replicate(S, design, simplify = FALSE),
    sizes = stats::setNames(
      replicate(S, rep(1, P), simplify = FALSE), ids
    ),
    K = K
  )
}

dktf_fit_features <- function(data, config) {
  subjects <- lapply(seq_along(data$ids), function(s) {
    dkge_subject(
      data$primary[[s]], data$designs[[s]], id = data$ids[[s]]
    )
  })
  fit <- suppressWarnings(dkge_fit(
    dkge_data(subjects), K = data$K, rank = 3,
    w_method = "none", effect_scaling = "none"
  ))
  contrast <- suppressWarnings(dkge_contrast(
    fit, c(0, 1, 0, 0, 0), method = "loso", align = FALSE
  ))
  features <- dkge_alignment_features(
    fit, contrast,
    independent_betas = data$alignment,
    independent_data_hash = paste0(
      "known-warp-independent-alignment-seed-", config$seed
    )
  )
  list(fit = fit, contrast = contrast, features = features)
}

dktf_mapper <- function(config, geometry_only = FALSE) {
  dkge_mapper_spec(
    "sinkhorn",
    epsilon = config$epsilon,
    max_iter = 8000L,
    tol = 1e-4,
    lambda_emb = if (geometry_only) 0 else config$lambda_emb,
    lambda_spa = if (geometry_only) 1 else config$lambda_spa,
    sigma_mm = config$sigma_mm,
    value_type = "intensive",
    warm_start = FALSE
  )
}

dktf_mni_average <- function(data, support) {
  maps <- operators <- vector("list", length(data$ids))
  for (s in seq_along(data$ids)) {
    fitted <- fit_mapper(
      dkge_mapper("knn", k = 1L, sigx = 1),
      subj_points = data$centroids[[s]],
      anchor_points = support$coordinates
    )
    maps[[s]] <- apply_mapper(fitted, data$heldout[[s]])
    operator <- matrix(0, nrow(data$centroids[[s]]), support$n_locations)
    operator[cbind(seq_len(nrow(operator)), fitted$idx[, 1L])] <-
      fitted$weights[, 1L]
    operators[[s]] <- operator
  }
  subj_values <- do.call(rbind, maps)
  rownames(subj_values) <- data$ids
  list(
    group = colMeans(subj_values),
    subj_values = subj_values,
    operators = operators
  )
}

dktf_operator_metrics <- function(operators, data) {
  latent_error <- mean(vapply(seq_along(operators), function(s) {
    operator <- as.matrix(operators[[s]])
    column_mass <- colSums(operator)
    induced <- as.numeric(crossprod(operator, data$latent)) /
      pmax(column_mass, .Machine$double.xmin)
    mean(abs(induced - data$latent))
  }, numeric(1)))
  point_spread <- mean(vapply(operators, function(operator) {
    operator <- as.matrix(operator)
    probabilities <- sweep(
      abs(operator), 2L,
      pmax(colSums(abs(operator)), .Machine$double.xmin), "/"
    )
    mean(apply(probabilities, 2L, function(p) {
      p <- p[p > 0]
      exp(-sum(p * log(p)))
    }))
  }, numeric(1)))
  list(latent_error = latent_error, point_spread = point_spread)
}

dktf_map_metrics <- function(group, truth, latent) {
  list(
    correlation = stats::cor(group, truth),
    rmse = sqrt(mean((group - truth)^2)),
    amplitude_ratio = stats::sd(group) / stats::sd(truth),
    positive_peak_error = abs(
      latent[[which.max(group)]] - latent[[which.max(truth)]]
    ),
    negative_peak_error = abs(
      latent[[which.min(group)]] - latent[[which.min(truth)]]
    )
  )
}

dktf_run <- function(config = dktf_config()) {
  started <- proc.time()[["elapsed"]]
  data <- dktf_generate(config)
  fitted <- dktf_fit_features(data, config)
  reference <- config$reference_subject
  reference_selection <- dkge_select_reference_subject(
    data$centroids,
    sizes = data$sizes,
    method = "explicit",
    reference_subject = reference,
    subject_ids = data$ids,
    provenance = list(
      fixed_by_protocol = TRUE,
      protocol_role = "known_warp_display_support"
    )
  )
  support <- dkge_reference_support(
    data$centroids[[reference]],
    labels = data$parcels,
    provenance = list(
      kind = "reference_subject_support",
      subject_id = data$ids[[reference]],
      reference_selection_hash = reference_selection$structural_hash,
      known_warp_validation = TRUE
    )
  )
  initializer <- dkge_functional_template(
    support,
    fitted$features$features[[reference]],
    feature_source = "independent",
    provenance = list(
      independent_data_hash =
        fitted$features$provenance$independent_data_hash,
      initialization_channel = "same_as_alignment_features",
      alignment_features_hash = fitted$features$structural_hash
    ),
    reference_selection = reference_selection
  )
  mapper <- dktf_mapper(config)
  template <- dkge_fit_functional_template(
    support,
    fitted$features,
    data$centroids,
    reference_selection = reference_selection,
    sizes = data$sizes,
    mapper = mapper,
    initialization = "supplied_independent",
    initial_template = initializer,
    max_iter = config$template_max_iter,
    tolerance = config$template_tolerance,
    objective_tolerance = config$template_objective_tolerance,
    update_rate = config$update_rate
  )
  aligned <- dkge_align_to_template(
    template, support, data$heldout,
    allow_nonconverged = TRUE
  )
  template_group <- colMeans(aligned$values[[1]])

  raw <- dkge:::.dkge_transport_to_medoid(
    mapper, data$heldout, fitted$features$features,
    data$centroids, data$sizes, reference,
    subject_ids = data$ids,
    preprocessing = list(source = "independent_raw_medoid_comparator")
  )
  geometry_features <- stats::setNames(lapply(seq_along(data$ids), function(s) {
    matrix(0, config$P, 1L)
  }), data$ids)
  geometry <- dkge:::.dkge_transport_to_medoid(
    dktf_mapper(config, geometry_only = TRUE),
    data$heldout, geometry_features,
    data$centroids, data$sizes, reference,
    subject_ids = data$ids,
    preprocessing = list(source = "geometry_only_comparator")
  )
  solver_converged <- c(
    iterative_template = all(vapply(
      template$fitting$diagnostics, function(x) isTRUE(x$converged), logical(1)
    )),
    raw_functional_medoid = isTRUE(
      raw$fitted_alignment$eligibility$solver_converged
    ),
    geometry_only = isTRUE(
      geometry$fitted_alignment$eligibility$solver_converged
    ),
    ordinary_mni_average = TRUE
  )
  if (!all(solver_converged)) {
    stop(
      "A known-warp arm failed its declared numerical solver contract: ",
      paste(names(solver_converged)[!solver_converged], collapse = ", "),
      call. = FALSE
    )
  }
  mni <- dktf_mni_average(data, support)

  groups <- list(
    iterative_template = template_group,
    raw_functional_medoid = colMeans(raw$subj_values),
    geometry_only = colMeans(geometry$subj_values),
    ordinary_mni_average = mni$group
  )
  aligned_rows <- list(
    iterative_template = aligned$values[[1]],
    raw_functional_medoid = raw$subj_values,
    geometry_only = geometry$subj_values,
    ordinary_mni_average = mni$subj_values
  )
  operators <- list(
    iterative_template = template$fitting$operators,
    raw_functional_medoid = raw$operators,
    geometry_only = geometry$operators,
    ordinary_mni_average = mni$operators
  )
  metrics <- do.call(rbind, lapply(names(groups), function(arm) {
    map <- dktf_map_metrics(groups[[arm]], data$truth, data$latent)
    operator <- dktf_operator_metrics(operators[[arm]], data)
    data.frame(
      arm = arm,
      correlation = map$correlation,
      rmse = map$rmse,
      amplitude_ratio = map$amplitude_ratio,
      positive_peak_error = map$positive_peak_error,
      negative_peak_error = map$negative_peak_error,
      latent_error = operator$latent_error,
      point_spread = operator$point_spread
    )
  }))
  gates <- list(
    convergence = isTRUE(template$fitting$converged),
    solver_convergence = all(solver_converged),
    independent_prediction = identical(
      fitted$features$feature_source, "independent"
    ),
    initializer_channel_consistent = identical(
      initializer$provenance$independent_data_hash,
      fitted$features$provenance$independent_data_hash
    ) && identical(
      initializer$provenance$alignment_features_hash,
      fitted$features$structural_hash
    ),
    beats_geometry_correlation =
      metrics$correlation[metrics$arm == "iterative_template"] >=
      metrics$correlation[metrics$arm == "geometry_only"] + 0.10,
    beats_mni_correlation =
      metrics$correlation[metrics$arm == "iterative_template"] >=
      metrics$correlation[metrics$arm == "ordinary_mni_average"] + 0.10,
    amplitude_is_auditable = all(is.finite(metrics$amplitude_ratio)) &&
      metrics$amplitude_ratio[metrics$arm == "iterative_template"] > 0.1 &&
      metrics$amplitude_ratio[metrics$arm == "iterative_template"] < 1.5,
    raw_medoid_comparator_reported =
      "raw_functional_medoid" %in% metrics$arm,
    approximate_not_promoted = identical(
      template$eligibility$status, "approximate"
    )
  )
  list(
    protocol = "dkge-functional-template-known-warp-v1",
    config = config,
    metrics = metrics,
    aligned_rows = aligned_rows,
    solver_converged = solver_converged,
    gates = gates,
    all_gates_pass = all(unlist(gates)),
    template = template,
    initializer = initializer,
    support = support,
    runtime_seconds = proc.time()[["elapsed"]] - started
  )
}
