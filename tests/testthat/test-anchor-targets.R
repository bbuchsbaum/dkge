test_that("dkge_anchor_targets_from_prototypes builds class matrix", {
  set.seed(1)
  anchors <- matrix(rnorm(15), nrow = 5)
  protos <- list(
    class1 = matrix(rnorm(6), nrow = 2),
    class2 = rnorm(3)
  )
  negs <- list(class1 = rnorm(3))
  W <- dkge_anchor_targets_from_prototypes(anchors, protos, negatives = negs)
  expect_equal(dim(W), c(2, nrow(anchors)))
  expect_equal(rownames(W), c("class1", "class2"))
})

test_that("dkge_anchor_targets_from_directions accepts list or matrix", {
  set.seed(2)
  anchors <- matrix(rnorm(12), nrow = 4)
  dirs <- matrix(rnorm(6), nrow = 2, dimnames = list(c("up", "down"), NULL))
  W1 <- dkge_anchor_targets_from_directions(anchors, dirs)
  expect_equal(dim(W1), c(2, nrow(anchors)))

  dir_list <- list(A = rnorm(ncol(anchors)), B = rnorm(ncol(anchors)))
  W2 <- dkge_anchor_targets_from_directions(anchors, dir_list)
  expect_equal(rownames(W2), c("A", "B"))
})
