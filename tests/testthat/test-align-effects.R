test_that("dkge_align_effects completes union via Nyström", {
  set.seed(1)
  make_psd <- function(eff) {
    p <- length(eff)
    L <- matrix(rnorm(p * p), p)
    K <- tcrossprod(L)
    dimnames(K) <- list(eff, eff)
    K
  }

  E1 <- c("A", "B", "C")
  E2 <- c("B", "C", "D")
  K1 <- make_psd(E1)
  K2 <- make_psd(E2)

  ref <- dkge_align_effects(
    K_list = list(s1 = K1, s2 = K2),
    effects = list(s1 = E1, s2 = E2),
    mode = "nystrom",
    ridge = 1e-6
  )

  expect_equal(names(ref$K_aligned), c("s1", "s2"))
  expect_equal(length(ref$effect_ids), 4L)
  expect_true(all(ref$effect_ids %in% c("A", "B", "C", "D")))

  for (Khat in ref$K_aligned) {
    expect_equal(nrow(Khat), length(ref$effect_ids))
    eig <- eigen(0.5 * (Khat + t(Khat)), symmetric = TRUE)$values
    expect_true(min(eig) > -1e-8)
  }

  expect_true(min(eigen(ref$G, symmetric = TRUE)$values) > -1e-8)

  mask1 <- ref$obs_mask[["s1"]]
  mask2 <- ref$obs_mask[["s2"]]
  expect_equal(mask1, ref$effect_ids %in% E1)
  expect_equal(mask2, ref$effect_ids %in% E2)

  diag_counts <- diag(ref$pair_counts)
  names(diag_counts) <- ref$effect_ids
  expect_equal(
    diag_counts[c("A", "B", "C", "D")],
    setNames(c(1L, 2L, 2L, 1L), c("A", "B", "C", "D"))
  )
})

test_that("shrinkage and intersection modes behave as expected", {
  set.seed(2)
  make_psd <- function(eff) {
    p <- length(eff)
    M <- matrix(rnorm(p * p), p)
    K <- tcrossprod(M)
    dimnames(K) <- list(eff, eff)
    K
  }

  E1 <- c("A", "B", "C")
  E2 <- c("B", "C", "D")
  K1 <- make_psd(E1)
  K2 <- make_psd(E2)

  shrink <- dkge_align_effects(
    list(K1, K2),
    list(E1, E2),
    mode = "shrinkage",
    alpha = 0.4,
    ensure_psd = FALSE
  )
  union_ids <- shrink$effect_ids
  G <- shrink$G

  obs1 <- union_ids %in% E1
  idx1 <- match(union_ids[obs1], E1)
  expected1 <- G
  expected1[obs1, obs1] <- 0.6 * K1[idx1, idx1, drop = FALSE] + 0.4 * G[obs1, obs1, drop = FALSE]
  expect_equal(shrink$K_aligned[[1]], expected1)

  obs2 <- union_ids %in% E2
  idx2 <- match(union_ids[obs2], E2)
  expected2 <- G
  expected2[obs2, obs2] <- 0.6 * K2[idx2, idx2, drop = FALSE] + 0.4 * G[obs2, obs2, drop = FALSE]
  expect_equal(shrink$K_aligned[[2]], expected2)

  inter <- dkge_align_effects(
    list(K1, K2),
    list(E1, E2),
    mode = "intersection"
  )
  expect_equal(inter$effect_ids, intersect(E1, E2))
  expect_equal(nrow(inter$K_aligned[[1]]), length(inter$effect_ids))
  expect_true(all(inter$obs_mask[[1]]))
  expect_true(all(inter$obs_mask[[2]]))
})

test_that("fold-aware alignment keeps training/test separation", {
  make_mat <- function(eff, seed) {
    set.seed(seed)
    p <- length(eff)
    L <- matrix(rnorm(p * p), p)
    K <- tcrossprod(L)
    dimnames(K) <- list(eff, eff)
    K
  }

  E1 <- c("A", "B")
  E2 <- c("B", "C", "D")
  E3 <- c("A")
  K1 <- make_mat(E1, 11)
  K2 <- make_mat(E2, 22)
  K3 <- make_mat(E3, 33)

  folds <- list(c(2))
  res <- dkge_align_effects(
    K_list = list(s1 = K1, s2 = K2, s3 = K3),
    effects = list(s1 = E1, s2 = E2, s3 = E3),
    folds = folds,
    mode = "nystrom"
  )

  fold1 <- res$folds$fold_1
  expect_equal(fold1$train_idx, c(1L, 3L))
  expect_equal(fold1$test_idx, c(2L))

  expect_equal(sort(fold1$effect_ids), sort(unique(c(E1, E3))))
  expect_equal(fold1$obs_mask[["s2"]], fold1$effect_ids %in% E2)

  diag_counts <- diag(fold1$pair_counts)
  names(diag_counts) <- fold1$effect_ids
  expect_equal(
    diag_counts[c("A", "B")],
    setNames(c(2L, 1L), c("A", "B"))
  )
  expect_true(all(diag_counts <= length(fold1$train_idx)))
})
