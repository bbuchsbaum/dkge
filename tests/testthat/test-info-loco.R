# test-info-loco.R

library(testthat)


test_that("LOCO proxy ranks planted anchor hotspots highest", {
  skip_on_cran()
  set.seed(7)

  ## ---- synthetic DKGE like in the Haufe test, but tailor mapper to concentrate signal
  S <- 16L; q <- 5L; r <- 3L
  U <- qr.Q(qr(matrix(rnorm(q * r), q, r)))
  K <- diag(q); Khalf <- diag(q); Kihalf <- diag(q); R <- diag(q)
  y <- factor(rep(c("A","B"), length.out = S), levels = c("A","B"))
  P_list <- sample(120:160, S, replace = TRUE)
  delta <- 0.7; sigma <- 0.45

  A_list <- lapply(seq_len(S), function(s) {
    P <- P_list[s]
    mu <- if (y[s] == "A") c(+delta, 0, 0) else c(-delta, 0, 0)
    matrix(rep(mu, each = P), nrow = P, byrow = FALSE) +
      matrix(rnorm(P * r, sd = sigma), P, r)
  })
  Btil <- lapply(A_list, function(A) U %*% t(A))
  fit <- structure(list(Btil = Btil, U = U, K = K, Khalf = Khalf, Kihalf = Kihalf, R = R), class = "dkge")

  # Mock helpers (A_s and mapper; aggregator not needed here)
  testthat::local_mocked_bindings(
    dkge_cluster_loadings = function(fit) {
      lapply(fit$Btil, function(Bs) t(Bs) %*% fit$K %*% fit$U)
    },
    apply_mapper = function(mapper_fit, x) as.numeric(mapper_fit$M %*% x),
    .package = "dkge"
  )

  # Build mapper so that a small set of "ROI" anchors capture the biggest |A_s[,1]| clusters
  Q <- 80L
  n_roi <- 6L
  roi <- 1:n_roi  # the first anchors are ROI by design
  mapper_fits <- vector("list", S)
  for (s in seq_len(S)) {
    P <- P_list[s]
    # score clusters by absolute signal along latent dim 1
    score <- abs(A_list[[s]][, 1])
    ord <- order(score, decreasing = TRUE)
    top <- ord[1:round(0.15 * P)]      # top 15% assigned to ROI anchors evenly
    rest <- setdiff(seq_len(P), top)

    M <- matrix(0, Q, P)
    # distribute 'top' clusters over ROI anchors
    chunks <- split(top, cut(seq_along(top), n_roi, labels = FALSE))
    for (j in seq_along(chunks)) {
      idx <- chunks[[j]]
      if (length(idx) > 0) M[roi[j], idx] <- 1 / length(idx)
    }
    # distribute remaining clusters uniformly over non-ROI anchors
    others <- setdiff(1:Q, roi)
    grp <- sample(others, length(rest), replace = TRUE)
    for (j in others) {
      idx <- which(grp == j)
      if (length(idx) > 0) M[j, idx] <- 1 / length(idx)
    }
    mapper_fits[[s]] <- list(M = M)
  }

  renderer <- list(
    anchors = matrix(0, Q, 3),
    mapper_fits = mapper_fits,
    graph = list(L = NULL, deg = rep(1, Q)),  # neighborhoods passed explicitly
    H = NULL,
    weights = rep(1, S)
  )

  # Cross-fitted classifier (subject-level, LDA)
  set.seed(101)
  clf <- dkge_cv_train_latent_classifier(fit, y, folds = 4, model = "lda", level = "subject", standardize = TRUE)

  # Neighborhoods: singletons so we get per-anchor LOCO
  neighborhoods <- lapply(seq_len(Q), function(j) j)

  loco <- dkge_info_map_loco(fit, clf, renderer, neighborhoods = neighborhoods, aggregate = "mean")

  # ROI anchors should dominate the top-10 LOCO scores
  top10 <- order(loco$loco_anchor, decreasing = TRUE)[1:10]
  hits <- sum(top10 %in% roi)
  expect_gte(hits, 5)   # >=5 of the top-10 belongs to ROI

  # Median importance in ROI > non-ROI
  med_roi    <- stats::median(loco$loco_anchor[roi])
  med_nonroi <- stats::median(loco$loco_anchor[-roi])
  expect_gt(med_roi, med_nonroi)
})
