# test-info-haufe.R

library(testthat)


test_that("Haufe: synthetic end-to-end recovers encoding pattern", {
  skip_on_cran()
  set.seed(1)

  ## ---------------- Synthetic DKGE fit ----------------
  S <- 18L                  # subjects
  q <- 5L                   # design effects
  r <- 3L                   # latent rank
  U <- qr.Q(qr(matrix(rnorm(q * r), q, r)))  # q×r orthonormal
  K <- diag(q); Khalf <- diag(q); Kihalf <- diag(q); R <- diag(q)

  P_list <- sample(150:200, S, replace = TRUE) # clusters per subject
  y <- factor(rep(c("A","B"), length.out = S), levels = c("A","B"))

  a_true <- c(1, 0, 0)      # planted encoding axis in latent space
  delta  <- 0.6             # subject-level mean shift along comp 1
  sigma  <- 0.5             # cluster noise

  # A_s (P_s×r) ~ N(mu_class, sigma^2 I)  -- also equals latent features Z_s
  A_list <- lapply(seq_len(S), function(s) {
    P <- P_list[s]
    mu <- if (y[s] == "A") c(+delta, 0, 0) else c(-delta, 0, 0)
    matrix(rep(mu, each = P), nrow = P, byrow = FALSE) +
      matrix(rnorm(P * r, sd = sigma), P, r)
  })

  # Btil_s = U %*% t(A_s)  so that  t(U) %*% Btil_s = A_s^T  (latent projection recovers A_s)
  Btil <- lapply(A_list, function(A) U %*% t(A))
  fit <- structure(list(
    Btil = Btil, U = U, K = K, Khalf = Khalf, Kihalf = Kihalf, R = R,
    weights = rep(1, S)
  ), class = "dkge")

  ## ---------------- Simple renderer (cluster -> anchors) ----------------
  Q <- 60L  # anchors
  mapper_fits <- lapply(seq_len(S), function(s) {
    P <- P_list[s]
    grp <- sample.int(Q, P, replace = TRUE)
    M <- matrix(0, Q, P)
    for (j in 1:Q) {
      idx <- which(grp == j)
      if (length(idx) > 0) M[j, idx] <- 1 / length(idx)  # row-normalized average
    }
    list(M = M)  # mapper 'fit'
  })
  renderer <- list(
    anchors = matrix(0, Q, 3),
    mapper_fits = mapper_fits,
    graph = list(L = NULL, deg = rep(1, Q)),  # no smoothing in this test
    H = diag(Q),                               # identity decode
    weights = rep(1, S)
  )

  ## ---------------- Mock small rendering helpers ----------------
  testthat::local_mocked_bindings(
    dkge_cluster_loadings = function(fit) {
      # A_s = t(Btil_s) K U  -> with our construction returns A_list
      lapply(fit$Btil, function(Bs) t(Bs) %*% fit$K %*% fit$U)
    },
    apply_mapper = function(mapper_fit, x) {
      as.numeric(mapper_fit$M %*% x)
    },
    dkge_anchor_aggregate = function(anchor_list, subj_weights, L = NULL, lambda = 0) {
      if (is.null(subj_weights)) subj_weights <- rep(1, length(anchor_list))
      Y <- do.call(rbind, anchor_list)    # S × Q
      w <- subj_weights / sum(subj_weights)
      list(y = colSums(Y * w))            # weighted mean across subjects
    },
    dkge_anchor_to_voxel_apply = function(H, y) as.numeric(H %*% y),
    .package = "dkge"
  )

  ## ---------------- Train cross-fitted latent classifier ----------------
  set.seed(42)
  clf <- dkge_cv_train_latent_classifier(
    fit, y, folds = 5, model = "lda", level = "subject", standardize = TRUE
  )

  # 1) Classifier orientation: avg cosine(beta, e1) must be strongly positive
  beta_mat <- do.call(rbind, clf$beta_by_subject)   # S × r
  e1 <- c(1, 0, 0)
  cosines <- apply(beta_mat, 1, function(b) {
    sum(b * e1) / (sqrt(sum(b^2)) * 1.0)
  })
  expect_true(all(cosines > 0))          # oriented by construction
  expect_gt(mean(cosines), 0.60)         # robust alignment

  ## ---------------- Haufe encoding map ----------------
  info <- dkge_info_map_haufe(
    fit, clf, renderer, Z_by_subject = NULL, inference = "none", to_vox = FALSE
  )

  # Ground-truth anchor map from planted encoding a_true
  subj_cluster_true <- lapply(seq_len(S), function(s) as.numeric(A_list[[s]] %*% a_true))
  subj_anchor_true  <- lapply(seq_len(S), function(s) as.numeric(mapper_fits[[s]]$M %*% subj_cluster_true[[s]]))
  y_true <- colMeans(do.call(rbind, subj_anchor_true))

  # 2) Haufe mean anchor correlates with planted truth
  cor_enc <- suppressWarnings(cor(info$mean_anchor, y_true))
  expect_gt(cor_enc, 0.70)

  ## ---------------- Haufe standardization invariance ----------------
  set.seed(42)
  clf_no_std <- dkge_cv_train_latent_classifier(
    fit, y, folds = 5, model = "lda", level = "subject", standardize = FALSE
  )
  info_no_std <- dkge_info_map_haufe(fit, clf_no_std, renderer, inference = "none", to_vox = FALSE)

  # 3) Same encoding up to scale -> correlation ~ 1
  cor_haufe_std <- suppressWarnings(cor(info$mean_anchor, info_no_std$mean_anchor))
  expect_gt(cor_haufe_std, 0.95)

  ## ---------------- Sanity: latent projection equals A_s ----------------
  Z_list <- dkge_project_clusters_to_latent(fit) # should be A_list under our construction
  expect_equal(Z_list[[1]], A_list[[1]], tolerance = 1e-8)
})
