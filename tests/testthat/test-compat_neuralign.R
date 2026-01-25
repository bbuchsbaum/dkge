test_that("dkge registers a neuralign aligner and fits K-Procrustes operators", {
  skip_if_not_installed("neuralign")

  # Ensure aligner is registered (even if .onLoad didn't run).
  dkge:::.dkge_register_neuralign_aligner()

  set.seed(1)
  q <- 6
  r <- 3
  effects <- paste0("eff", seq_len(q))
  K <- diag(q)
  dimnames(K) <- list(effects, effects)

  # Build K-orthonormal bases U_s (q x r), then store X_s = t(U_s) (r x q).
  make_subject <- function() {
    W <- matrix(rnorm(q * r), q, r, dimnames = list(effects, NULL))
    U <- dkge_k_orthonormalize(W, K)
    X <- t(U)
    colnames(X) <- effects
    X
  }
  X_list <- list(
    "sub-01" = make_subject(),
    "sub-02" = make_subject(),
    "sub-03" = make_subject()
  )

  adat <- neuralign::AlignmentData(
    data = X_list,
    design = list(K = K, effects = effects),
    space = "dkge_effect_basis"
  )

  res <- neuralign::fit_alignment(
    adat,
    method = "dkge",
    reference = "consensus",
    compute_quality = FALSE
  )

  model <- neuralign::get_model(res)
  expect_identical(names(model@transforms), names(X_list))
  dims_ok <- vapply(model@transforms, function(A) all(dim(A) == c(r, r)), logical(1))
  expect_true(all(dims_ok))
})

