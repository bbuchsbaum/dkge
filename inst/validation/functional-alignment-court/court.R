#!/usr/bin/env Rscript

# Frozen functional-alignment Type-I court. This is validation code, not a
# package API. See frozen-protocol.json before changing any design constant.

`%fa_or%` <- function(x, y) if (is.null(x)) y else x

fa_court_root <- function() {
  here <- normalizePath(getwd(), mustWork = TRUE)
  repeat {
    desc <- file.path(here, "DESCRIPTION")
    if (file.exists(desc) && any(grepl("^Package: dkge$", readLines(desc)))) {
      return(here)
    }
    parent <- dirname(here)
    if (identical(parent, here)) break
    here <- parent
  }
  stop("Could not locate the DKGE repository root.", call. = FALSE)
}

fa_court_protocol_path <- function() {
  file.path(fa_court_root(), "inst", "validation",
            "functional-alignment-court", "frozen-protocol.json")
}

fa_court_output_dir <- function() {
  file.path(fa_court_root(), "data-raw", "functional-alignment-court")
}

fa_court_hash <- function(x) {
  digest::digest(x, algo = "sha256", serialize = TRUE)
}

fa_court_hash_file <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

fa_court_source_tree_hash <- function(root = fa_court_root()) {
  roots <- c("DESCRIPTION", "NAMESPACE", "R", "src")
  paths <- unlist(lapply(roots, function(path) {
    absolute <- file.path(root, path)
    if (dir.exists(absolute)) {
      list.files(absolute, recursive = TRUE, full.names = TRUE,
                 all.files = TRUE, no.. = TRUE)
    } else if (file.exists(absolute)) {
      absolute
    } else {
      character()
    }
  }), use.names = FALSE)
  paths <- paths[file.info(paths)$isdir %in% FALSE]
  paths <- paths[!grepl("\\.(o|so|dll|dylib)$", paths, ignore.case = TRUE)]
  paths <- sort(paths, method = "radix")
  relative <- substring(paths, nchar(root) + 2L)
  records <- paste(relative, vapply(paths, fa_court_hash_file, character(1)),
                   sep = "=")
  digest::digest(paste(records, collapse = "\n"), algo = "sha256",
                 serialize = FALSE)
}

fa_court_wilson <- function(x, n, level = 0.95) {
  if (!is.finite(n) || n <= 0) return(c(lower = NA_real_, upper = NA_real_))
  z <- stats::qnorm(1 - (1 - level) / 2)
  p <- x / n
  den <- 1 + z^2 / n
  centre <- (p + z^2 / (2 * n)) / den
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / den
  c(lower = max(0, centre - half), upper = min(1, centre + half))
}

fa_court_base_config <- function() {
  list(
    dgp = "gaussian_separable",
    S = 8L,
    q = 6L,
    estimation_rank = 3L,
    kernel_rank = 6L,
    contrast_family_dimension = 1L,
    eigengap = "large",
    nuisance_magnitude = 1,
    spatial_correlation = 0.35,
    spatial_heteroscedasticity = 1.5,
    parcel_count_heterogeneity = TRUE,
    base_parcels = 9L,
    epsilon = 0.08,
    functional_cost_weight = 0.65,
    sigma_spatial = 0.18,
    w_method = "mfa_sigma1",
    w_tau = 0.3,
    reference_subject = 1L
  )
}

fa_court_grid <- function() {
  b <- fa_court_base_config()
  variants <- list(
    baseline = list(),
    small_S = list(S = 8L),
    medium_S = list(S = 16L),
    large_S = list(S = 32L),
    rank_2 = list(estimation_rank = 2L),
    rank_4 = list(estimation_rank = 4L),
    kernel_rank_4 = list(kernel_rank = 4L),
    family_2 = list(contrast_family_dimension = 2L),
    small_gap = list(eigengap = "small"),
    weak_nuisance = list(nuisance_magnitude = 0.5),
    strong_nuisance = list(nuisance_magnitude = 2),
    independent_space = list(spatial_correlation = 0),
    correlated_space = list(spatial_correlation = 0.65),
    homoscedastic = list(spatial_heteroscedasticity = 1),
    heteroscedastic = list(spatial_heteroscedasticity = 3),
    equal_parcels = list(parcel_count_heterogeneity = FALSE),
    low_epsilon = list(epsilon = 0.03),
    high_epsilon = list(epsilon = 0.15),
    geometry_cost = list(functional_cost_weight = 0),
    functional_cost = list(functional_cost_weight = 1),
    heavy_tail = list(dgp = "student_t_spatial_heteroscedastic")
  )
  rows <- lapply(seq_along(variants), function(i) {
    cfg <- utils::modifyList(b, variants[[i]])
    cfg$cell <- names(variants)[[i]]
    cfg$cell_index <- i
    cfg$exact_oracle <- i %in% c(1L, 3L, 5L, 9L, 15L, 21L)
    as.data.frame(cfg, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

fa_court_as_config <- function(row) {
  x <- as.list(row)
  ints <- c("S", "q", "estimation_rank", "kernel_rank",
            "contrast_family_dimension", "base_parcels",
            "reference_subject", "cell_index")
  x[ints] <- lapply(x[ints], as.integer)
  nums <- c("nuisance_magnitude", "spatial_correlation",
            "spatial_heteroscedasticity", "epsilon",
            "functional_cost_weight", "sigma_spatial", "w_tau")
  x[nums] <- lapply(x[nums], as.numeric)
  x$parcel_count_heterogeneity <- as.logical(x$parcel_count_heterogeneity)
  x$exact_oracle <- as.logical(x$exact_oracle)
  x
}

fa_court_noise <- function(n, dgp) {
  if (identical(dgp, "student_t_spatial_heteroscedastic")) {
    stats::rt(n, df = 5) / sqrt(5 / 3)
  } else {
    stats::rnorm(n)
  }
}

fa_court_spatial_draw <- function(P, rho, hetero, dgp) {
  x <- seq(0, 1, length.out = P)
  Sigma <- if (rho <= 0) diag(P) else exp(-abs(outer(x, x, "-")) / rho)
  L <- chol(Sigma + diag(1e-8, P))
  scale <- seq(1 / sqrt(hetero), sqrt(hetero), length.out = P)
  as.numeric(fa_court_noise(P, dgp) %*% L) * scale
}

fa_court_generate <- function(cfg, seed) {
  set.seed(seed)
  S <- cfg$S
  q <- cfg$q
  m <- cfg$contrast_family_dimension
  if (m >= q || cfg$kernel_rank < m || cfg$estimation_rank > cfg$kernel_rank) {
    stop("Court configuration violates rank/contrast constraints.", call. = FALSE)
  }
  offsets <- if (isTRUE(cfg$parcel_count_heterogeneity)) {
    rep(c(-2L, -1L, 0L, 1L, 2L), length.out = S)
  } else {
    rep(0L, S)
  }
  P <- pmax(5L, cfg$base_parcels + offsets)
  effects <- paste0("effect", seq_len(q))
  designs <- replicate(S, {
    X <- diag(q)
    colnames(X) <- effects
    X
  }, simplify = FALSE)
  K <- diag(c(rep(1, cfg$kernel_rank), rep(0, q - cfg$kernel_rank)))
  dimnames(K) <- list(effects, effects)

  build_replicate <- function(independent = FALSE) {
    lapply(seq_len(S), function(s) {
      ps <- P[[s]]
      x <- seq(0, 1, length.out = ps)
      B <- matrix(0, q, ps, dimnames = list(effects, paste0("p", seq_len(ps))))
      for (j in seq_len(m)) {
        B[j, ] <- fa_court_spatial_draw(
          ps, cfg$spatial_correlation,
          cfg$spatial_heteroscedasticity, cfg$dgp
        )
      }
      nuisance_ids <- seq.int(m + 1L, q)
      nuisance_scale <- if (identical(cfg$eigengap, "small")) {
        seq(1, 0.94, length.out = length(nuisance_ids))
      } else {
        seq(1.8, 0.45, length.out = length(nuisance_ids))
      }
      for (jj in seq_along(nuisance_ids)) {
        j <- nuisance_ids[[jj]]
        phase <- 0.17 * s + 0.11 * jj
        shared <- sin(2 * pi * (jj * x + phase))
        amp <- cfg$nuisance_magnitude * nuisance_scale[[jj]]
        B[j, ] <- amp * shared + 0.35 * fa_court_spatial_draw(
          ps, cfg$spatial_correlation,
          cfg$spatial_heteroscedasticity, cfg$dgp
        )
      }
      if (independent) {
        # Independent replicate: same geometry/topological nuisance law, fresh
        # tested and nuisance noise. It shares no tested beta realization.
        B <- B + matrix(stats::rnorm(q * ps, sd = 0.05), q, ps)
      }
      B
    })
  }
  B <- build_replicate(FALSE)
  B_independent <- build_replicate(TRUE)
  centroids <- lapply(seq_len(S), function(s) {
    x <- seq(0, 1, length.out = P[[s]])
    cbind(x = x + stats::rnorm(P[[s]], sd = 0.01),
          y = 0.03 * sin(2 * pi * x + 0.2 * s),
          z = 0)
  })
  sizes <- lapply(seq_len(S), function(s) {
    exp(stats::rnorm(P[[s]], sd = 0.2))
  })
  contrasts <- diag(q)[, seq_len(m), drop = FALSE]
  colnames(contrasts) <- paste0("c", seq_len(m))
  list(B = B, B_independent = B_independent, designs = designs, K = K,
       contrasts = contrasts, centroids = centroids, sizes = sizes,
       P = P, tested_rows = seq_len(m), seed = seed)
}

fa_court_fit <- function(dat, cfg, independent = FALSE) {
  dkge_fit(
    if (independent) dat$B_independent else dat$B,
    designs = dat$designs,
    K = dat$K,
    rank = cfg$estimation_rank,
    w_method = cfg$w_method,
    w_tau = cfg$w_tau
  )
}

fa_court_legacy_values <- function(fit, contrasts) {
  S <- length(fit$Btil)
  r <- ncol(fit$U)
  out <- lapply(seq_len(ncol(contrasts)), function(j) vector("list", S))
  names(out) <- colnames(contrasts)
  for (s in seq_len(S)) {
    train <- setdiff(seq_len(S), s)
    Chat <- Reduce(`+`, Map(function(w, M) w * M,
                            fit$weights[train], fit$contribs[train]))
    eig <- eigen((Chat + t(Chat)) / 2, symmetric = TRUE)
    U <- fit$Kihalf %*% eig$vectors[, seq_len(r), drop = FALSE]
    U <- dkge_k_orthonormalize(U, fit$K)
    Bmodel <- dkge:::.dkge_fit_model_btil(fit, s)
    A <- t(Bmodel) %*% fit$K %*% U
    for (j in seq_len(ncol(contrasts))) {
      ct <- backsolve(fit$R, contrasts[, j], transpose = FALSE)
      alpha <- as.numeric(crossprod(U, fit$K %*% ct))
      out[[j]][[s]] <- as.numeric(A %*% alpha)
    }
  }
  ids <- fit$subject_ids %fa_or% paste0("subject", seq_len(S))
  out <- lapply(out, function(x) stats::setNames(x, ids))
  out
}

fa_court_kernel_factor <- function(K) {
  ee <- eigen((K + t(K)) / 2, symmetric = TRUE)
  keep <- ee$values > max(ee$values, 0) * 1e-10
  sweep(ee$vectors[, keep, drop = FALSE], 2L,
        sqrt(ee$values[keep]), "*")
}

fa_court_residual_features <- function(fit, contrast_obj) {
  L <- fa_court_kernel_factor(fit$K)
  receipts <- contrast_obj$metadata$alignment_receipts
  out <- vector("list", length(receipts))
  diagnostics <- vector("list", length(receipts))
  for (s in seq_along(receipts)) {
    receipt <- receipts[[s]]
    Bmodel <- dkge:::.dkge_fit_model_btil(
      fit, s, voxel_weights = receipt$preprocessing$voxel_weights
    )
    G <- t(Bmodel) %*% L
    Gamma <- do.call(cbind, lapply(receipt$alphas, function(alpha) {
      as.numeric(crossprod(L, receipt$basis %*% alpha))
    }))
    qr_gamma <- qr(Gamma, tol = 1e-10)
    if (qr_gamma$rank < ncol(Gamma)) {
      stop("Prototype residual contrast span is rank deficient.", call. = FALSE)
    }
    Qg <- qr.Q(qr_gamma)[, seq_len(qr_gamma$rank), drop = FALSE]
    F <- G %*% (diag(ncol(G)) - tcrossprod(Qg))
    out[[s]] <- F
    reconstructed <- vapply(seq_len(ncol(Gamma)), function(j) {
      max(abs(as.numeric(G %*% Gamma[, j]) -
                contrast_obj$values[[j]][[s]]))
    }, numeric(1))
    diagnostics[[s]] <- list(
      available_dimension = ncol(G) - qr_gamma$rank,
      reconstruction_error = max(reconstructed),
      orthogonality_error = norm(F %*% Gamma, "F")
    )
  }
  names(out) <- names(receipts)
  attr(out, "diagnostics") <- diagnostics
  out
}

fa_court_mapper <- function(cfg, geometry_only = FALSE) {
  fw <- if (geometry_only) 0 else cfg$functional_cost_weight
  dkge_mapper_spec(
    "sinkhorn",
    epsilon = cfg$epsilon,
    max_iter = 5000L,
    tol = 1e-4,
    lambda_emb = fw,
    lambda_spa = 1 - fw,
    sigma_mm = cfg$sigma_spatial,
    value_type = "intensive",
    warm_start = FALSE
  )
}

fa_court_map <- function(values, features, dat, cfg, arm,
                         geometry_only = FALSE) {
  mapper <- fa_court_mapper(cfg, geometry_only = geometry_only)
  cache <- NULL
  Y <- vector("list", length(values))
  first <- NULL
  for (j in seq_along(values)) {
    tr <- dkge:::.dkge_transport_to_medoid(
      mapper, values[[j]], features, dat$centroids, dat$sizes,
      cfg$reference_subject, transport_cache = cache,
      subject_ids = names(features) %fa_or% paste0("subject", seq_along(features)),
      preprocessing = list(court_arm = arm)
    )
    if (is.null(first)) first <- tr
    Y[[j]] <- tr$subj_values
    cache <- tr$fitted_alignment
  }
  names(Y) <- names(values)
  list(Y = Y, alignment = cache, diagnostics = first$diagnostics,
       operators = first$operators)
}

fa_court_t_vector <- function(Y) {
  unlist(lapply(Y, function(M) {
    colMeans(M) / (apply(M, 2, stats::sd) / sqrt(nrow(M)) + 1e-12)
  }), use.names = FALSE)
}

fa_court_mean_vector <- function(Y) {
  unlist(lapply(Y, colMeans), use.names = FALSE)
}

fa_court_signs <- function(S, B, seed) {
  set.seed(seed)
  max_unique <- 2^(S - 1L) - 1L
  if (S <= 16L && B >= max_unique) {
    idx <- 0:(max_unique - 1L)
    bits <- vapply(idx, function(i) as.integer(intToBits(i))[seq_len(S - 1L)],
                   integer(S - 1L))
    return(rbind(1, ifelse(bits == 1L, 1, -1)))
  }
  signs <- matrix(sample(c(-1, 1), S * B, replace = TRUE), S, B)
  signs[1, ] <- 1
  signs
}

fa_court_signed_Y <- function(Y, sign) {
  lapply(Y, function(M) sign * M)
}

fa_court_operator_distance <- function(a, b) {
  if (is.null(a) || is.null(b) || length(a) != length(b)) return(NA_real_)
  mean(vapply(seq_along(a), function(i) {
    A <- as.matrix(a[[i]])
    B <- as.matrix(b[[i]])
    if (!identical(dim(A), dim(B))) return(NA_real_)
    norm(A - B, "F") / max(norm(B, "F"), 1e-12)
  }, numeric(1)), na.rm = TRUE)
}

fa_court_convergence <- function(diagnostics) {
  ok <- vapply(diagnostics, function(x) isTRUE(x$converged %fa_or% TRUE),
               logical(1))
  mean(ok)
}

fa_court_frozen_test <- function(mapped, signs, alpha = 0.05) {
  observed_t <- fa_court_t_vector(mapped$Y)
  observed_side <- abs(observed_t)
  B <- ncol(signs)
  null <- vapply(seq_len(B), function(b) {
    max(abs(fa_court_t_vector(fa_court_signed_Y(mapped$Y, signs[, b]))))
  }, numeric(1))
  coordinate_null <- vapply(seq_len(B), function(b) {
    abs(fa_court_t_vector(fa_court_signed_Y(mapped$Y, signs[, b])))
  }, numeric(length(observed_t)))
  p_fwer <- (1 + sum(null >= max(observed_side))) / (B + 1)
  p_unadjusted <- vapply(seq_along(observed_side), function(j) {
    (1 + sum(coordinate_null[j, ] >= observed_side[[j]])) / (B + 1)
  }, numeric(1))
  means <- fa_court_mean_vector(mapped$Y)
  list(
    p_fwer = p_fwer,
    reject_fwer = p_fwer <= alpha,
    reject_unadjusted = any(p_unadjusted <= alpha),
    min_p_unadjusted = min(p_unadjusted),
    effect_bias = mean(means),
    effect_abs_bias = mean(abs(means)),
    observed_max_t = max(observed_side),
    convergence = fa_court_convergence(mapped$diagnostics)
  )
}

fa_court_apply_raw_sign <- function(B_list, tested_rows, sign) {
  Map(function(B, z) {
    out <- B
    out[tested_rows, ] <- z * out[tested_rows, , drop = FALSE]
    out
  }, B_list, sign)
}

fa_court_observed_objects <- function(dat, cfg) {
  fit <- fa_court_fit(dat, cfg)
  contrast <- dkge_contrast(fit, dat$contrasts, method = "loso", align = FALSE)
  legacy_values <- fa_court_legacy_values(fit, dat$contrasts)
  full_features <- dkge:::.dkge_fit_subject_loadings(fit)
  fold_features <- dkge:::.dkge_fold_receipt_loadings(
    fit, contrast, cfg$reference_subject
  )$loadings
  residual_features <- fa_court_residual_features(fit, contrast)
  independent_fit <- fa_court_fit(dat, cfg, independent = TRUE)
  independent_features <- dkge:::.dkge_fit_subject_loadings(independent_fit)
  geometry_features <- lapply(dat$P, function(P) matrix(0, P, 1L))
  ids <- fit$subject_ids %fa_or% paste0("subject", seq_len(cfg$S))
  names(geometry_features) <- ids
  list(
    fit = fit,
    contrast = contrast,
    values = contrast$values,
    legacy_values = legacy_values,
    features = list(
      legacy_fullfit_functional = full_features,
      l1_fixed_fullfit_functional = full_features,
      l2_fold_loading = fold_features,
      geometry_only = geometry_features,
      kernel_image_residual_prototype = residual_features,
      independent_alignment = independent_features
    )
  )
}

fa_court_run_frozen_arms <- function(dat, cfg, observed, signs) {
  arm_names <- names(observed$features)
  mapped <- vector("list", length(arm_names))
  names(mapped) <- arm_names
  for (arm in arm_names) {
    vals <- if (identical(arm, "legacy_fullfit_functional")) {
      observed$legacy_values
    } else {
      observed$values
    }
    mapped[[arm]] <- fa_court_map(
      vals, observed$features[[arm]], dat, cfg, arm,
      geometry_only = identical(arm, "geometry_only")
    )
  }
  geometry_ops <- mapped$geometry_only$operators
  rows <- lapply(arm_names, function(arm) {
    res <- fa_court_frozen_test(mapped[[arm]], signs)
    data.frame(
      arm = arm,
      p_fwer = res$p_fwer,
      reject_fwer = res$reject_fwer,
      reject_unadjusted = res$reject_unadjusted,
      min_p_unadjusted = res$min_p_unadjusted,
      effect_bias = res$effect_bias,
      effect_abs_bias = res$effect_abs_bias,
      observed_max_t = res$observed_max_t,
      operator_sensitivity = fa_court_operator_distance(
        mapped[[arm]]$operators, geometry_ops
      ),
      convergence = res$convergence,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

fa_court_full_reestimation <- function(dat, cfg, signs, observed = NULL) {
  observed <- observed %fa_or% fa_court_observed_objects(dat, cfg)
  observed_map <- fa_court_map(
    observed$values, observed$features$l2_fold_loading, dat, cfg,
    "full_permutation_reestimation"
  )
  obs_t <- abs(fa_court_t_vector(observed_map$Y))
  B <- ncol(signs)
  null_max <- numeric(B)
  null_coordinate <- matrix(NA_real_, nrow = length(obs_t), ncol = B)
  sensitivity <- numeric(B)
  convergence <- numeric(B)
  for (b in seq_len(B)) {
    perm_dat <- dat
    perm_dat$B <- fa_court_apply_raw_sign(
      dat$B, dat$tested_rows, signs[, b]
    )
    perm_fit <- fa_court_fit(perm_dat, cfg)
    perm_contrast <- dkge_contrast(
      perm_fit, dat$contrasts, method = "loso", align = FALSE
    )
    perm_features <- dkge:::.dkge_fold_receipt_loadings(
      perm_fit, perm_contrast, cfg$reference_subject
    )$loadings
    perm_map <- fa_court_map(
      perm_contrast$values, perm_features, dat, cfg,
      "full_permutation_reestimation"
    )
    stat <- abs(fa_court_t_vector(perm_map$Y))
    null_coordinate[, b] <- stat
    null_max[[b]] <- max(stat)
    sensitivity[[b]] <- fa_court_operator_distance(
      perm_map$operators, observed_map$operators
    )
    convergence[[b]] <- fa_court_convergence(perm_map$diagnostics)
  }
  p_fwer <- (1 + sum(null_max >= max(obs_t))) / (B + 1)
  p_unadjusted <- vapply(seq_along(obs_t), function(j) {
    (1 + sum(null_coordinate[j, ] >= obs_t[[j]])) / (B + 1)
  }, numeric(1))
  means <- fa_court_mean_vector(observed_map$Y)
  data.frame(
    arm = "full_permutation_reestimation",
    p_fwer = p_fwer,
    reject_fwer = p_fwer <= 0.05,
    reject_unadjusted = any(p_unadjusted <= 0.05),
    min_p_unadjusted = min(p_unadjusted),
    effect_bias = mean(means),
    effect_abs_bias = mean(abs(means)),
    observed_max_t = max(obs_t),
    operator_sensitivity = mean(sensitivity),
    convergence = mean(convergence),
    stringsAsFactors = FALSE
  )
}

fa_court_gap_metrics <- function(observed, cfg) {
  gaps <- vapply(observed$contrast$metadata$alignment_receipts, function(x) {
    x$basis_diagnostics$relative_eigengap
  }, numeric(1))
  gap <- mean(gaps[is.finite(gaps) & gaps > 0])
  if (!is.finite(gap)) gap <- NA_real_
  c(relative_eigengap = gap,
    latent_span_predictor = if (is.finite(gap) && gap > 0) {
      1 / (cfg$S * gap)
    } else {
      NA_real_
    })
}

fa_court_run_one <- function(cfg, seed, n_perm, include_full = FALSE,
                             full_first = FALSE) {
  dat <- fa_court_generate(cfg, seed)
  signs <- fa_court_signs(cfg$S, n_perm, seed + 900000L)
  observed <- fa_court_observed_objects(dat, cfg)
  exact <- NULL
  exact_elapsed <- NA_real_
  if (include_full && full_first) {
    started <- proc.time()[["elapsed"]]
    exact <- fa_court_full_reestimation(dat, cfg, signs, observed)
    exact_elapsed <- proc.time()[["elapsed"]] - started
  }
  frozen <- fa_court_run_frozen_arms(dat, cfg, observed, signs)
  if (include_full && is.null(exact)) {
    started <- proc.time()[["elapsed"]]
    exact <- fa_court_full_reestimation(dat, cfg, signs, observed)
    exact_elapsed <- proc.time()[["elapsed"]] - started
  }
  result <- if (is.null(exact)) frozen else rbind(frozen, exact)
  gaps <- fa_court_gap_metrics(observed, cfg)
  result$seed <- seed
  result$cell <- cfg$cell %fa_or% "pilot"
  result$cell_index <- cfg$cell_index %fa_or% 0L
  result$S <- cfg$S
  result$estimation_rank <- cfg$estimation_rank
  result$kernel_rank <- cfg$kernel_rank
  result$contrast_family_dimension <- cfg$contrast_family_dimension
  result$eigengap_setting <- cfg$eigengap
  result$relative_eigengap <- unname(gaps[["relative_eigengap"]])
  result$latent_span_predictor <- unname(gaps[["latent_span_predictor"]])
  result$nuisance_magnitude <- cfg$nuisance_magnitude
  result$spatial_correlation <- cfg$spatial_correlation
  result$spatial_heteroscedasticity <- cfg$spatial_heteroscedasticity
  result$parcel_count_heterogeneity <- cfg$parcel_count_heterogeneity
  result$epsilon <- cfg$epsilon
  result$functional_cost_weight <- cfg$functional_cost_weight
  result$dgp <- cfg$dgp
  result$n_perm <- ncol(signs)
  result$exact_elapsed_seconds <- exact_elapsed
  result
}

fa_court_select_budget <- function(seconds_per_exact_perm) {
  formal_cores <- 4L
  candidates <- expand.grid(
    n_sim = c(40L, 64L, 96L),
    n_perm = c(31L, 63L, 127L)
  )
  candidates$predicted_seconds <- seconds_per_exact_perm *
    candidates$n_sim * candidates$n_perm * 6 * 2.5 / formal_cores
  eligible <- candidates$predicted_seconds <= 1200
  if (!any(eligible)) {
    pick <- which(candidates$n_sim == 40L & candidates$n_perm == 31L)[[1]]
    resource_limited <- TRUE
  } else {
    score <- candidates$n_sim * candidates$n_perm
    score[!eligible] <- -Inf
    best <- which(score == max(score))
    pick <- best[which.max(candidates$n_sim[best])]
    resource_limited <- FALSE
  }
  out <- candidates[pick, , drop = FALSE]
  out$seconds_per_exact_perm <- seconds_per_exact_perm
  out$resource_limited <- resource_limited
  out$exact_cells <- 6L
  out$formal_cores <- formal_cores
  out
}

fa_court_summarize <- function(raw) {
  keys <- interaction(raw$cell, raw$arm, drop = TRUE, lex.order = TRUE)
  pieces <- split(raw, keys)
  out <- lapply(pieces, function(x) {
    n <- nrow(x)
    fwer <- sum(x$reject_fwer)
    unadj <- sum(x$reject_unadjusted)
    ci <- fa_court_wilson(fwer, n)
    data.frame(
      cell = x$cell[[1]],
      arm = x$arm[[1]],
      n_sim = n,
      fwer_rejections = fwer,
      fwer_rate = fwer / n,
      fwer_ci_lower = ci[["lower"]],
      fwer_ci_upper = ci[["upper"]],
      unadjusted_rate = unadj / n,
      mean_effect_bias = mean(x$effect_bias),
      mean_effect_abs_bias = mean(x$effect_abs_bias),
      mean_operator_sensitivity = mean(x$operator_sensitivity, na.rm = TRUE),
      convergence = mean(x$convergence),
      mean_relative_eigengap = mean(x$relative_eigengap, na.rm = TRUE),
      mean_latent_span_predictor = mean(x$latent_span_predictor, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

fa_court_verdict <- function(summary) {
  exact <- summary[summary$arm == "full_permutation_reestimation", , drop = FALSE]
  negative <- summary[summary$arm %in% c("geometry_only", "independent_alignment"),
                      , drop = FALSE]
  exact_ok <- nrow(exact) > 0 && all(exact$fwer_ci_lower <= 0.05 &
                                      exact$fwer_ci_upper >= 0.05 &
                                      exact$convergence >= 0.95)
  negative_ok <- nrow(negative) > 0 && all(
    negative$fwer_ci_lower <= 0.05 & negative$fwer_ci_upper >= 0.05 &
      negative$fwer_ci_upper <= 0.12 & negative$convergence >= 0.95
  )
  inflated <- summary[
    summary$arm %in% c("legacy_fullfit_functional",
                       "l1_fixed_fullfit_functional", "l2_fold_loading",
                       "kernel_image_residual_prototype") &
      summary$fwer_ci_lower > 0.05 & summary$fwer_rate > 0.10,
    c("cell", "arm", "fwer_rate", "fwer_ci_lower", "fwer_ci_upper"),
    drop = FALSE
  ]
  list(
    candidate_valid = exact_ok && negative_ok,
    exact_oracle_gate = exact_ok,
    negative_control_gate = negative_ok,
    materially_inflated = inflated,
    claim = if (exact_ok && negative_ok) {
      "Court controls passed for the frozen scope; same-data modes remain approximate."
    } else {
      "Court controls failed or are incomplete; no inferential promotion is permitted."
    }
  )
}

fa_court_manifest <- function(mode, elapsed, extra = character()) {
  root <- fa_court_root()
  script <- file.path(root, "inst", "validation",
                      "functional-alignment-court", "court.R")
  sha <- tryCatch(system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
                  error = function(e) NA_character_)
  c(
    paste0("mode=", mode),
    paste0("generated_at=", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    paste0("elapsed_seconds=", format(elapsed, digits = 12)),
    paste0("source_sha=", sha[[1]] %fa_or% NA_character_),
    paste0("source_tree_sha256=", fa_court_source_tree_hash(root)),
    paste0("script_sha256=", fa_court_hash_file(script)),
    paste0("protocol_sha256=", fa_court_hash_file(fa_court_protocol_path())),
    extra,
    paste0("R=", R.version.string),
    paste0("platform=", R.version$platform)
  )
}

fa_court_write_manifest <- function(lines, path) {
  writeLines(lines, con = path, useBytes = TRUE)
}

fa_court_run_pilot <- function() {
  outdir <- fa_court_output_dir()
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  cfg <- fa_court_base_config()
  cfg$cell <- "pilot"
  cfg$cell_index <- 0L
  started <- proc.time()[["elapsed"]]
  rows <- lapply(seq_len(2L), function(i) {
    # full_first is the contractual point: exact re-estimation runs before any
    # frozen-plan arm is timed or inspected.
    fa_court_run_one(cfg, 7202601L + i - 1L, 15L,
                     include_full = TRUE, full_first = TRUE)
  })
  raw <- do.call(rbind, rows)
  elapsed <- proc.time()[["elapsed"]] - started
  exact <- raw[raw$arm == "full_permutation_reestimation", , drop = FALSE]
  seconds_per_perm <- sum(exact$exact_elapsed_seconds) / sum(exact$n_perm)
  budget <- fa_court_select_budget(seconds_per_perm)
  utils::write.csv(raw, file.path(outdir, "pilot-raw.csv"), row.names = FALSE)
  utils::write.csv(budget, file.path(outdir, "formal-budget.csv"), row.names = FALSE)
  fa_court_write_manifest(
    fa_court_manifest(
      "pilot", elapsed,
      c(paste0("pilot_n_sim=2"), paste0("pilot_n_perm=15"),
        paste0("seconds_per_exact_perm=", seconds_per_perm),
        paste0("selected_n_sim=", budget$n_sim),
        paste0("selected_n_perm=", budget$n_perm),
        paste0("formal_cores=", budget$formal_cores),
        paste0("predicted_formal_seconds=", budget$predicted_seconds),
        paste0("resource_limited=", budget$resource_limited))
    ),
    file.path(outdir, "pilot-manifest.txt")
  )
  invisible(list(raw = raw, budget = budget))
}

fa_court_run_formal <- function() {
  outdir <- fa_court_output_dir()
  budget_path <- file.path(outdir, "formal-budget.csv")
  if (!file.exists(budget_path)) {
    stop("Run the exact-oracle pilot before the formal court.", call. = FALSE)
  }
  budget <- utils::read.csv(budget_path, stringsAsFactors = FALSE)
  n_sim <- as.integer(budget$n_sim[[1]])
  n_perm <- as.integer(budget$n_perm[[1]])
  requested_cores <- if ("formal_cores" %in% names(budget)) {
    as.integer(budget$formal_cores[[1]])
  } else {
    4L
  }
  declared_cores <- suppressWarnings(as.integer(
    Sys.getenv("DKGE_COURT_CORES", unset = NA_character_)
  ))
  available_cores <- if (is.finite(declared_cores) && declared_cores > 0L) {
    declared_cores
  } else {
    parallel::detectCores(logical = FALSE)
  }
  if (!is.finite(available_cores)) available_cores <- 1L
  formal_cores <- if (.Platform$OS.type == "windows") {
    1L
  } else {
    max(1L, min(requested_cores, as.integer(available_cores)))
  }
  grid <- fa_court_grid()
  started <- proc.time()[["elapsed"]]
  tasks <- expand.grid(
    sim = seq_len(n_sim),
    cell_i = seq_len(nrow(grid)),
    KEEP.OUT.ATTRS = FALSE
  )
  run_task <- function(task_i) {
    cell_i <- tasks$cell_i[[task_i]]
    sim <- tasks$sim[[task_i]]
    cfg <- fa_court_as_config(grid[cell_i, , drop = FALSE])
    seed <- 7300000L + 10000L * cell_i + sim
    tryCatch(
      list(
        ok = TRUE,
        value = fa_court_run_one(
          cfg, seed, n_perm,
          include_full = isTRUE(cfg$exact_oracle), full_first = FALSE
        )
      ),
      error = function(e) list(
        ok = FALSE,
        cell = cfg$cell,
        sim = sim,
        seed = seed,
        message = conditionMessage(e)
      )
    )
  }
  task_ids <- seq_len(nrow(tasks))
  task_results <- if (formal_cores > 1L) {
    parallel::mclapply(
      task_ids, run_task, mc.cores = formal_cores,
      mc.preschedule = TRUE, mc.set.seed = FALSE
    )
  } else {
    lapply(task_ids, run_task)
  }
  failed <- which(!vapply(task_results, `[[`, logical(1), "ok"))
  if (length(failed)) {
    detail <- vapply(task_results[failed], function(x) {
      sprintf("cell=%s sim=%d seed=%d: %s",
              x$cell, x$sim, x$seed, x$message)
    }, character(1))
    stop("Formal court task failure(s):\n", paste(detail, collapse = "\n"),
         call. = FALSE)
  }
  rows <- lapply(task_results, `[[`, "value")
  raw <- do.call(rbind, rows)
  summary <- fa_court_summarize(raw)
  verdict <- fa_court_verdict(summary)
  elapsed <- proc.time()[["elapsed"]] - started
  utils::write.csv(grid, file.path(outdir, "formal-grid.csv"), row.names = FALSE)
  utils::write.csv(raw, file.path(outdir, "formal-raw.csv"), row.names = FALSE)
  utils::write.csv(summary, file.path(outdir, "formal-summary.csv"), row.names = FALSE)
  utils::write.csv(verdict$materially_inflated,
                   file.path(outdir, "formal-inflation-flags.csv"),
                   row.names = FALSE)
  latent_rows <- raw[raw$arm %in% c("l2_fold_loading",
                                    "kernel_image_residual_prototype"), ]
  latent_model <- stats::lm(
    operator_sensitivity ~ latent_span_predictor + arm,
    data = latent_rows
  )
  utils::write.csv(as.data.frame(summary(latent_model)$coefficients),
                   file.path(outdir, "latent-span-model.csv"))
  fa_court_write_manifest(
    fa_court_manifest(
      "formal", elapsed,
      c(paste0("n_sim=", n_sim), paste0("n_perm=", n_perm),
        paste0("formal_cores=", formal_cores),
        paste0("requested_formal_cores=", requested_cores),
        paste0("candidate_valid=", verdict$candidate_valid),
        paste0("exact_oracle_gate=", verdict$exact_oracle_gate),
        paste0("negative_control_gate=", verdict$negative_control_gate),
        paste0("claim=", verdict$claim))
    ),
    file.path(outdir, "formal-manifest.txt")
  )
  invisible(list(raw = raw, summary = summary, verdict = verdict))
}

fa_court_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (!length(args) || !args[[1]] %in% c("pilot", "formal")) {
    stop("Usage: court.R pilot|formal", call. = FALSE)
  }
  loader_script <- file.path(
    fa_court_root(), "inst", "validation",
    "functional-alignment-certification", "evidence-utils.R"
  )
  if (!file.exists(loader_script)) {
    stop("Exact-source package loader is unavailable.", call. = FALSE)
  }
  source(loader_script, local = TRUE)
  dkfa_load_source_package(fa_court_root(), quiet = TRUE)
  if (identical(args[[1]], "pilot")) {
    fa_court_run_pilot()
  } else {
    fa_court_run_formal()
  }
}

court_args <- commandArgs(trailingOnly = TRUE)
if (length(court_args) && court_args[[1]] %in% c("pilot", "formal")) {
  fa_court_main(court_args)
}
