# Optional integration with the `neuralign` package
#
# This file defines a minimal adapter that lets dkge register an aligner into
# neuralign's registry at load time, without creating a dependency cycle
# (neuralign should not Import dkge).

.dkge_neuralign_capabilities <- function() {
  list(
    supports_cv = TRUE,
    cv_axes = c("subject"),
    needs_geometry = FALSE,
    needs_design = TRUE,
    returns_invertible = TRUE,
    transform_type = "orthogonal",
    mass_preserving = FALSE,
    returns = "operator",
    supports_new_subject = TRUE,
    supports_new_data = TRUE,
    reference_types = c("subject", "consensus", "template")
  )
}

.dkge_neuralign_is_matrixish <- function(x) {
  is.matrix(x) || inherits(x, "Matrix")
}

.dkge_neuralign_default_effects <- function(q) {
  paste0("effect", seq_len(q))
}

.dkge_neuralign_extract_design <- function(design, q) {
  if (is.null(design) || length(design) == 0L) {
    stop("dkge(neuralign): AlignmentData@design must provide at least K and effects", call. = FALSE)
  }

  if (is.matrix(design) || inherits(design, "Matrix")) {
    K <- as.matrix(design)
    effects <- colnames(K) %||% rownames(K) %||% .dkge_neuralign_default_effects(nrow(K))
    return(list(K = K, effects = as.character(effects)))
  }

  if (!is.list(design)) {
    stop("dkge(neuralign): AlignmentData@design must be a list (or a K matrix)", call. = FALSE)
  }

  K <- design$K %||% design$design_kernel %||% design$kernel
  if (is.list(K) && !is.null(K$K)) {
    K <- K$K
  }
  if (!.dkge_neuralign_is_matrixish(K)) {
    stop("dkge(neuralign): design$K must be a matrix (or list with $K)", call. = FALSE)
  }
  K <- as.matrix(K)
  if (nrow(K) != ncol(K)) {
    stop("dkge(neuralign): design$K must be square", call. = FALSE)
  }
  if (ncol(K) != q) {
    stop(sprintf(
      "dkge(neuralign): design kernel dimension mismatch: K is %d x %d, data has q=%d effects",
      nrow(K), ncol(K), q
    ), call. = FALSE)
  }

  effects <- design$effects %||% colnames(K) %||% rownames(K)
  if (is.null(effects)) {
    effects <- .dkge_neuralign_default_effects(q)
  }
  effects <- as.character(effects)
  if (length(effects) != q) {
    stop("dkge(neuralign): design$effects must have length q (number of effects)", call. = FALSE)
  }

  # Ensure K has dimnames consistent with effects.
  if (is.null(rownames(K)) || is.null(colnames(K))) {
    dimnames(K) <- list(effects, effects)
  } else if (!identical(rownames(K), effects) || !identical(colnames(K), effects)) {
    if (!setequal(rownames(K), effects) || !setequal(colnames(K), effects)) {
      stop("dkge(neuralign): K dimnames must match design$effects", call. = FALSE)
    }
    K <- K[effects, effects, drop = FALSE]
  }

  list(K = K, effects = effects)
}

.dkge_neuralign_validate_data_list <- function(data_list, effects, context) {
  if (!is.list(data_list) || !length(data_list) || is.null(names(data_list))) {
    stop(context, ": expected a named list of matrices", call. = FALSE)
  }

  bad <- names(data_list)[!vapply(data_list, .dkge_neuralign_is_matrixish, logical(1))]
  if (length(bad)) {
    stop(context, ": all subjects must be matrices; invalid: ", paste(bad, collapse = ", "), call. = FALSE)
  }

  dims <- vapply(data_list, function(m) paste(dim(m), collapse = "x"), character(1))
  if (length(unique(dims)) != 1L) {
    stop(context, ": all subjects must have the same matrix dimensions (r x q). Got: ",
      paste(sprintf("%s=%s", names(dims), dims), collapse = ", "),
      call. = FALSE
    )
  }

  q <- length(effects)
  one <- as.matrix(data_list[[1]])
  if (ncol(one) != q) {
    if (nrow(one) == q) {
      stop(context, ": matrices appear transposed (q x r). Expected r x q (components x effects). ",
        "If you are passing dkge bases U (q x r), use t(U) instead.", call. = FALSE
      )
    }
    stop(context, ": expected r x q matrices with q = ", q, " effects", call. = FALSE)
  }

  # If columns have names, enforce/reorder to match effects.
  for (subj in names(data_list)) {
    m <- as.matrix(data_list[[subj]])
    cn <- colnames(m)
    if (!is.null(cn)) {
      cn <- as.character(cn)
      if (!identical(cn, effects)) {
        if (!setequal(cn, effects)) {
          stop(context, ": colnames for '", subj, "' do not match design$effects", call. = FALSE)
        }
        m <- m[, effects, drop = FALSE]
      }
    } else {
      colnames(m) <- effects
    }
    data_list[[subj]] <- m
  }

  data_list
}

.dkge_neuralign_reference_to_U <- function(reference,
                                          U_list,
                                          train_subjects,
                                          K,
                                          allow_reflection) {
  if (is.character(reference) && length(reference) == 1L) {
    if (reference == "consensus") {
      U_train <- U_list[train_subjects]
      cons <- dkge_consensus_basis_K(U_train, K, allow_reflection = allow_reflection)
      return(cons$U)
    }
    if (reference %in% names(U_list)) {
      return(U_list[[reference]])
    }
    stop(
      "dkge(neuralign): unknown reference subject id (not found in data): ",
      reference,
      call. = FALSE
    )
  }

  if (.dkge_neuralign_is_matrixish(reference)) {
    ref <- as.matrix(reference)
    q <- nrow(K)
    if (nrow(ref) == q) {
      # q x r basis
      return(ref)
    }
    if (ncol(ref) == q) {
      # r x q "X" representation
      return(t(ref))
    }
    stop("dkge(neuralign): template reference has incompatible dims; expected qxr (U) or rxq (t(U))", call. = FALSE)
  }

  stop("dkge(neuralign): unsupported reference; expected 'consensus', a subject id, or a template matrix", call. = FALSE)
}

.dkge_neuralign_fit <- function(data,
                                reference = "consensus",
                                train_idx = NULL,
                                allow_reflection = TRUE,
                                ...) {
  if (is.null(train_idx)) train_idx <- seq_along(data@subjects)
  subjects <- data@subjects
  train_subjects <- subjects[train_idx]

  data_list <- neuralign::get_data_list(data)
  q <- ncol(as.matrix(data_list[[1]]))

  design <- .dkge_neuralign_extract_design(data@design, q = q)
  data_list <- .dkge_neuralign_validate_data_list(data_list, design$effects, "dkge(neuralign) fit")

  U_list <- lapply(data_list, function(X) t(X)) # q x r
  U_list <- lapply(U_list, function(U) dkge_k_orthonormalize(U, design$K))

  U_ref <- .dkge_neuralign_reference_to_U(
    reference = reference,
    U_list = U_list,
    train_subjects = train_subjects,
    K = design$K,
    allow_reflection = allow_reflection
  )
  U_ref <- dkge_k_orthonormalize(U_ref, design$K)

  transforms <- lapply(U_list, function(U) {
    pr <- dkge_procrustes_K(U_ref, U, design$K, allow_reflection = allow_reflection)
    t(pr$R) # (r x r) left-multiply operator on X = t(U)
  })
  names(transforms) <- names(data_list)

  space_spec <- data@space
  if (is.null(space_spec)) {
    space_spec <- "dkge_effect_basis"
  }

  list(
    transforms = transforms,
    reference_data = t(U_ref), # r x q
    space_from = space_spec,
    space_to = space_spec,
    method_state = list(
      U_ref = U_ref,
      K = design$K,
      effects = design$effects,
      allow_reflection = allow_reflection
    )
  )
}

.dkge_neuralign_apply <- function(fit_result, new_data, ...) {
  data_list <- neuralign::get_data_list(new_data)
  if (length(data_list) != 1L) {
    stop("dkge(neuralign): apply expects exactly one new subject", call. = FALSE)
  }
  subj_id <- names(data_list)[[1]]

  st <- fit_result$method_state %||% list()
  K <- st$K
  effects <- st$effects
  U_ref <- st$U_ref
  allow_reflection <- isTRUE(st$allow_reflection)

  if (is.null(K) || is.null(effects) || is.null(U_ref)) {
    stop("dkge(neuralign): missing method_state (expected K/effects/U_ref)", call. = FALSE)
  }

  data_list <- .dkge_neuralign_validate_data_list(data_list, effects, "dkge(neuralign) apply")
  X <- as.matrix(data_list[[1]])
  U <- dkge_k_orthonormalize(t(X), K)

  pr <- dkge_procrustes_K(U_ref, U, K, allow_reflection = allow_reflection)
  op <- t(pr$R)

  list(
    transforms = stats::setNames(list(op), subj_id),
    reference_data = fit_result$reference_data,
    space_from = fit_result$space_from,
    space_to = fit_result$space_to,
    method_state = fit_result$method_state
  )
}

.dkge_register_neuralign_aligner <- function() {
  if (!requireNamespace("neuralign", quietly = TRUE)) return(invisible(FALSE))

  if (isTRUE(neuralign::is_aligner_registered("dkge"))) return(invisible(TRUE))

  neuralign::register_aligner(
    name = "dkge",
    fit_fn = .dkge_neuralign_fit,
    apply_fn = .dkge_neuralign_apply,
    capabilities = .dkge_neuralign_capabilities(),
    package = "dkge",
    description = "DKGE effect-space alignment (K-Procrustes; operators are rxr)",
    version = as.character(utils::packageVersion("dkge"))
  )

  invisible(TRUE)
}

