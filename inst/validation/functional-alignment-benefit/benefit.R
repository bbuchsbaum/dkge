#!/usr/bin/env Rscript

# Functional-alignment EFFICACY study.
#
# The frozen Type-I court (inst/validation/functional-alignment-court) answers
# "does alignment leak?". It cannot answer "does alignment help": in its DGP
# every parcel sits at `seq(0, 1, length.out = P)` and centroids are that same
# coordinate plus 0.01 jitter, so ground-truth correspondence is recoverable
# from geometry alone, and its functional patterns carry a subject-dependent
# phase (`0.17 * s`) that points *away* from the truth.
#
# This module supplies the complementary DGP: a subject-specific monotone
# spatial warp makes anatomy genuinely uninformative about correspondence,
# while a phase-locked functional signature makes function genuinely
# informative. Setting `warp_amount = 0` recovers the anatomically corresponded
# regime, where functional alignment must NOT help.
#
# Ground truth is the latent coordinate `u`; anatomy is `x = phi_s(u)`.

`%fab_or%` <- function(x, y) if (is.null(x)) y else x

fab_base_config <- function() {
  list(
    S = 12L,
    q = 6L,
    n_shared = 3L,
    base_parcels = 24L,
    parcel_heterogeneity = TRUE,
    warp_amount = 0.9,
    signature_snr = 4,
    centroid_noise = 0.004,
    epsilon = 0.05,
    functional_cost_weight = 0.65,
    sigma_spatial = 0.15,
    estimation_rank = 3L,
    reference_subject = 1L
  )
}

# Subject-specific anatomical displacement. amount = 0 is the identity.
#
# The displacement is deliberately NOT monotone. A monotone warp is useless
# here: on a 1-D support with balanced uniform masses, entropic OT recovers the
# monotone rearrangement, which matches rank to rank -- i.e. exactly the true
# latent correspondence -- no matter how severe the warp. Geometry alone then
# scores near-perfectly and no functional signal can add anything. Breaking
# monotonicity is what makes anatomical position genuinely ambiguous, so that
# functional identity is the only thing that can resolve it.
fab_displace <- function(u, amount, coefs) {
  if (amount <= 0) return(u)
  d <- rowSums(vapply(seq_along(coefs), function(j) {
    coefs[[j]] * sin(j * pi * u)
  }, numeric(length(u))))
  u + amount * d
}

# Fraction of adjacent latent pairs whose anatomical order is preserved. At 1
# the map is monotone and geometry can solve the problem by rank alone.
fab_monotone_fraction <- function(u, x) {
  mean(diff(x) > 0)
}

# Functional signature indexed by LATENT coordinate and shared across
# subjects: this is the common signal alignment is supposed to exploit.
fab_signature <- function(u, n_shared) {
  vapply(seq_len(n_shared), function(j) {
    sin(2 * pi * (j * u + 0.13 * j))
  }, numeric(length(u)))
}

fab_generate <- function(cfg, seed) {
  set.seed(seed)
  S <- cfg$S
  q <- cfg$q
  stopifnot(cfg$n_shared + 1L <= q)

  offsets <- if (isTRUE(cfg$parcel_heterogeneity)) {
    rep(c(-4L, -2L, 0L, 2L, 4L), length.out = S)
  } else {
    rep(0L, S)
  }
  P <- pmax(8L, cfg$base_parcels + offsets)
  effects <- paste0("effect", seq_len(q))
  ids <- paste0("subject", seq_len(S))

  latent <- lapply(seq_len(S), function(s) seq(0, 1, length.out = P[[s]]))
  coefs <- lapply(seq_len(S), function(s) stats::rnorm(3L, sd = 1))
  anatomy <- lapply(seq_len(S), function(s) {
    fab_displace(latent[[s]], cfg$warp_amount, coefs[[s]])
  })

  sd_noise <- 1 / cfg$signature_snr
  draw <- function() {
    lapply(seq_len(S), function(s) {
      ps <- P[[s]]
      Sig <- fab_signature(latent[[s]], cfg$n_shared)
      B <- matrix(stats::rnorm(q * ps, sd = sd_noise), q, ps,
                  dimnames = list(effects, paste0("p", seq_len(ps))))
      rows <- seq.int(2L, cfg$n_shared + 1L)
      B[rows, ] <- t(Sig) + B[rows, , drop = FALSE]
      B
    })
  }

  B <- draw()
  B_align <- draw()
  B_eval <- draw()

  centroids <- lapply(seq_len(S), function(s) {
    cbind(x = anatomy[[s]] +
            stats::rnorm(P[[s]], sd = cfg$centroid_noise),
          y = 0, z = 0)
  })
  sizes <- lapply(seq_len(S), function(s) rep(1, P[[s]]))
  designs <- replicate(S, {
    X <- diag(q)
    colnames(X) <- effects
    X
  }, simplify = FALSE)
  K <- diag(q)
  dimnames(K) <- list(effects, effects)

  names(B) <- names(B_align) <- names(B_eval) <- ids
  names(centroids) <- names(sizes) <- names(latent) <- ids

  list(B = B, B_align = B_align, B_eval = B_eval,
       designs = designs, K = K, centroids = centroids, sizes = sizes,
       latent = latent, anatomy = anatomy, P = P, ids = ids,
       monotone_fraction = mean(vapply(seq_len(S), function(s) {
         fab_monotone_fraction(latent[[s]], anatomy[[s]])
       }, numeric(1))),
       n_shared = cfg$n_shared, seed = seed)
}

# Alignment features come from an independent replicate: efficacy is measured
# in the leak-free mode so that "does it help" is not entangled with "does it
# leak".  Geometry-only features are a constant column, as in the court.
fab_features <- function(dat, cfg, geometry_only = FALSE) {
  rows <- seq.int(2L, cfg$n_shared + 1L)
  out <- lapply(seq_along(dat$B_align), function(s) {
    if (geometry_only) {
      matrix(0, dat$P[[s]], 1L)
    } else {
      t(dat$B_align[[s]][rows, , drop = FALSE])
    }
  })
  names(out) <- dat$ids
  out
}

fab_mapper <- function(cfg, geometry_only = FALSE) {
  fw <- if (geometry_only) 0 else cfg$functional_cost_weight
  dkge_mapper_spec(
    "sinkhorn",
    epsilon = cfg$epsilon,
    max_iter = 12000L,
    tol = 1e-4,
    lambda_emb = fw,
    lambda_spa = 1 - fw,
    sigma_mm = cfg$sigma_spatial,
    value_type = "intensive",
    warm_start = FALSE
  )
}

# Transport a per-subject field and return operators alongside values.
fab_transport <- function(dat, cfg, values, geometry_only = FALSE) {
  mapper <- fab_mapper(cfg, geometry_only = geometry_only)
  features <- fab_features(dat, cfg, geometry_only = geometry_only)
  dkge:::.dkge_transport_to_medoid(
    mapper, values, features, dat$centroids, dat$sizes,
    cfg$reference_subject,
    subject_ids = dat$ids,
    preprocessing = list(study = "benefit",
                         arm = if (geometry_only) "geometry" else "functional")
  )
}

# --- metrics ---------------------------------------------------------------

# Column-normalised operator: column j holds the source weights feeding
# target j, so the induced latent position of target j is a weighted mean.
fab_induced_latent <- function(operator, u_source) {
  w <- as.matrix(operator)
  totals <- colSums(w)
  totals[totals <= 0] <- NA_real_
  as.numeric(crossprod(w, u_source)) / totals
}

# Primary recovery metric: mean absolute latent-coordinate error. Robust to
# how much the plan smooths, because it is a position, not an amplitude.
fab_latent_error <- function(transport, dat, cfg) {
  u_ref <- dat$latent[[cfg$reference_subject]]
  errs <- vapply(seq_along(dat$ids), function(s) {
    if (identical(s, as.integer(cfg$reference_subject))) return(NA_real_)
    induced <- fab_induced_latent(transport$operators[[s]], dat$latent[[s]])
    mean(abs(induced - u_ref), na.rm = TRUE)
  }, numeric(1))
  mean(errs, na.rm = TRUE)
}

# Chance baseline: every target draws uniformly from the source.
fab_chance_error <- function(dat, cfg) {
  u_ref <- dat$latent[[cfg$reference_subject]]
  errs <- vapply(seq_along(dat$ids), function(s) {
    if (identical(s, as.integer(cfg$reference_subject))) return(NA_real_)
    mean(abs(mean(dat$latent[[s]]) - u_ref))
  }, numeric(1))
  mean(errs, na.rm = TRUE)
}

# Secondary metric: does better correspondence actually reconstruct the field?
# Uses the held-out replicate B_eval, never the one that built the features.
fab_reconstruction <- function(transport_fun, dat, cfg) {
  u_ref <- dat$latent[[cfg$reference_subject]]
  truth <- fab_signature(u_ref, cfg$n_shared)[, 1]
  row <- 2L
  values <- lapply(seq_along(dat$ids), function(s) {
    as.numeric(dat$B_eval[[s]][row, ])
  })
  names(values) <- dat$ids
  tr <- transport_fun(values)
  Y <- tr$subj_values
  per_subject <- vapply(seq_len(nrow(Y)), function(s) {
    stats::cor(Y[s, ], truth)
  }, numeric(1))
  list(mean_subject_cor = mean(per_subject, na.rm = TRUE),
       group_cor = stats::cor(colMeans(Y), truth),
       operators = tr$operators)
}

# Effective smoothing of each operator, so that any reconstruction advantage
# can be checked not to be a smoothing artefact.
fab_operator_diffusion <- function(operator) {
  w <- as.matrix(operator)
  totals <- colSums(w)
  keep <- totals > 0
  if (!any(keep)) return(NA_real_)
  w <- sweep(w[, keep, drop = FALSE], 2L, totals[keep], "/")
  ent <- apply(w, 2L, function(col) {
    col <- col[col > 0]
    -sum(col * log(col))
  })
  mean(exp(ent))
}

fab_run_arm <- function(dat, cfg, geometry_only) {
  values <- lapply(seq_along(dat$ids), function(s) {
    as.numeric(dat$B[[s]][2L, ])
  })
  names(values) <- dat$ids
  tr <- fab_transport(dat, cfg, values, geometry_only = geometry_only)
  recon <- fab_reconstruction(
    function(vals) fab_transport(dat, cfg, vals, geometry_only = geometry_only),
    dat, cfg
  )
  diffusion <- mean(vapply(tr$operators, fab_operator_diffusion, numeric(1)),
                    na.rm = TRUE)
  list(
    arm = if (geometry_only) "geometry_only" else "functional",
    latent_error = fab_latent_error(tr, dat, cfg),
    mean_subject_cor = recon$mean_subject_cor,
    group_cor = recon$group_cor,
    diffusion = diffusion,
    transport = tr
  )
}

fab_compare <- function(cfg, seed) {
  dat <- fab_generate(cfg, seed)
  list(
    data = dat,
    chance = fab_chance_error(dat, cfg),
    functional = fab_run_arm(dat, cfg, geometry_only = FALSE),
    geometry = fab_run_arm(dat, cfg, geometry_only = TRUE)
  )
}

# --- processing symmetry ---------------------------------------------------

# A reference subject that is special-cased to an identity operator is neither
# smoothed nor attenuated, while every other subject is both. That is a
# categorical processing asymmetry: it biases the group map toward one
# participant's topography and inflates that participant's leverage on the
# across-subject variance at every target location.
fab_operator_symmetry <- function(transport, dat, cfg) {
  rms <- function(z) sqrt(mean(z^2))
  u_ref <- dat$latent[[cfg$reference_subject]]
  f_ref <- fab_signature(u_ref, cfg$n_shared)[, 1]
  diffusion <- vapply(transport$operators, fab_operator_diffusion, numeric(1))
  attenuation <- vapply(seq_along(dat$ids), function(s) {
    f_s <- fab_signature(dat$latent[[s]], cfg$n_shared)[, 1]
    rms(as.numeric(crossprod(as.matrix(transport$operators[[s]]), f_s))) /
      rms(f_ref)
  }, numeric(1))
  ref <- as.integer(cfg$reference_subject)
  others <- setdiff(seq_along(dat$ids), ref)
  z_score <- function(x) {
    (x[[ref]] - mean(x[others])) / (stats::sd(x[others]) + 1e-12)
  }
  list(
    diffusion = diffusion,
    attenuation = attenuation,
    diffusion_ratio = max(diffusion) / min(diffusion),
    attenuation_ratio = max(attenuation) / min(attenuation),
    reference_diffusion_z = z_score(diffusion),
    reference_attenuation_z = z_score(attenuation)
  )
}

# --- predeclared gates -----------------------------------------------------
#
# Calibrated on seeds 201-215 before any assertion was written; see
# data-raw/functional-alignment-benefit/ for the calibration record.
fab_gates <- function() {
  list(
    warp_corresponded = 0,
    warp_idiosyncratic = 0.25,
    calibration_seeds = 201:215,
    # Regime A (idiosyncratic topography). Observed over the calibration
    # seeds: functional wins 15/15; latent-error ratio geometry/functional
    # mean 1.39, min 1.18; correlation gain mean 0.315, min 0.118.
    min_error_ratio = 1.10,
    min_cor_gain = 0.08,
    # Regime B (anatomically corresponded). Observed: correlation gain is
    # negative on every seed, range [-0.023, -0.007]; latent error is worse
    # by ~5x. Functional alignment must not appear to help, and its cost must
    # stay small.
    max_cor_gain_corresponded = 0.02,
    max_cor_loss_corresponded = 0.05,
    # Any regime-A advantage must not be bought with extra smoothing.
    # Observed: functional diffusion 2.14 vs geometry 2.96, i.e. sharper.
    max_diffusion_excess = 0.5,
    # Recovery must beat drawing uniformly from the source parcels.
    # Observed at warp 0.25: chance 0.263 vs functional 0.154.
    min_chance_ratio = 1.5
  )
}
