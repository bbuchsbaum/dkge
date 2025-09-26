library(testthat)

context("design kernel")

test_that("cell-basis kernel is PSD and sized correctly", {
  factors <- list(A = list(L = 3, type = "nominal"),
                  B = list(L = 4, type = "ordinal", l = 1.5))
  terms <- list("A", "B", c("A", "B"))
  Kres <- design_kernel(factors, terms = terms,
                        rho = c("A" = 1, "B" = 1, "A:B" = 0.4),
                        basis = "cell", normalize = "unit_trace")
  expect_equal(dim(Kres$K), c(12, 12))
  eigvals <- eigen((Kres$K + t(Kres$K))/2, symmetric = TRUE)$values
  expect_true(all(eigvals >= -1e-8))
})

test_that("effect-basis map matches block metadata", {
  factors <- list(A = list(L = 3, type = "nominal"),
                  B = list(L = 4, type = "ordinal", l = 1.0))
  Ctrs <- sum_contrasts(c(A = 3, B = 4))
  Kres <- design_kernel(factors, terms = list("A", "B", c("A", "B")),
                        basis = "effect", contrasts = Ctrs)
  T <- Kres$info$map
  blocks <- Kres$info$blocks
  expect_true(is.matrix(T))
  expect_equal(sum(lengths(blocks)), ncol(T))
  expect_equal(sort(names(blocks)), sort(Kres$info$term_names))
  expect_true(all(unlist(blocks) %in% seq_len(ncol(T))))
})

test_that("normalisation removes dependence on rho scaling", {
  factors <- list(A = list(L = 3, type = "nominal"))
  K1 <- design_kernel(factors, basis = "cell", normalize = "unit_trace")$K
  K2 <- design_kernel(factors, rho = c("A" = 5), basis = "cell", normalize = "unit_trace")$K
  expect_equal(sum(diag(K1)), sum(diag(K2)))
})

test_that("contrasts helpers have expected shapes", {
  sc <- sum_contrasts(c(A = 4))
  hc <- helmert_contrasts(c(A = 4))
  expect_equal(dim(sc$A), c(4, 3))
  expect_equal(dim(hc$A), c(4, 3))
  expect_equal(unname(crossprod(hc$A)), diag(3), tolerance = 1e-8)
})

test_that("kernel roots and alignment behave", {
  factors <- list(A = list(L = 3, type = "nominal"))
  K <- design_kernel(factors, basis = "cell", normalize = "none")$K
  roots <- kernel_roots(K)
  expect_equal(dim(roots$Khalf), dim(K))
  expect_equal(dim(roots$Kihalf), dim(K))
  reconstructed <- roots$Khalf %*% roots$Khalf
  expect_equal(reconstructed, K, tolerance = 1e-6)
  expect_equal(kernel_alignment(K, K), 1, tolerance = 1e-8)
})

test_that("kernel_roots reports clamped eigenvalues", {
  K <- diag(c(1, 1e-12, 0))
  roots <- kernel_roots(K, jitter = NULL)
  expect_equal(roots$n_clamped, 1L)
  expect_equal(roots$rank, 3L)
  expect_true(all(roots$evals >= 0))
})
