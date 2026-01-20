library(testthat)

# test-design-kernel.R
# Tests for design kernel construction, including edge cases

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

# ---------------------------------------------------------------------------
# Edge case tests: Small kernels
# ---------------------------------------------------------------------------

test_that("design_kernel handles 1-level factor (1x1 kernel)", {
  factors <- list(A = list(L = 1, type = "nominal"))
  K_res <- design_kernel(factors, basis = "cell", normalize = "none")

  expect_equal(dim(K_res$K), c(1, 1))
  expect_true(K_res$K[1, 1] > 0)  # Should be positive
})

test_that("design_kernel handles 2-level factor (2x2 kernel)", {
  factors <- list(A = list(L = 2, type = "nominal"))
  K_res <- design_kernel(factors, basis = "cell", normalize = "none")

  expect_equal(dim(K_res$K), c(2, 2))
  expect_true(isSymmetric(K_res$K, tol = 1e-10))

  # For nominal, diagonal should be equal (both 1 + jitter + rho0)
  expect_equal(K_res$K[1, 1], K_res$K[2, 2], tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# Edge case tests: Input validation
# ---------------------------------------------------------------------------

test_that("design_kernel rejects negative rho values", {
  factors <- list(A = list(L = 3, type = "nominal"))
  expect_error(design_kernel(factors, rho = c("A" = -1)),
               "non-negative|rho")
})

test_that("design_kernel rejects missing L for discrete factors", {
  factors <- list(A = list(type = "nominal"))  # Missing L
  expect_error(design_kernel(factors), "L")
})

test_that("design_kernel rejects missing values for continuous factors", {
  factors <- list(A = list(type = "continuous"))  # Missing values
  expect_error(design_kernel(factors), "values")
})

test_that("design_kernel rejects unnamed factors", {
  factors <- list(list(L = 3, type = "nominal"))  # Unnamed
  expect_error(design_kernel(factors), "named")
})

# ---------------------------------------------------------------------------
# Edge case tests: kernel_roots edge cases
# ---------------------------------------------------------------------------

test_that("kernel_roots handles diagonal matrices", {
  K <- diag(c(2, 1, 0.5))
  roots <- kernel_roots(K)

  expect_equal(diag(roots$Khalf), sqrt(c(2, 1, 0.5)), tolerance = 1e-10)
  reconstructed <- roots$Khalf %*% roots$Khalf
  expect_equal(reconstructed, K, tolerance = 1e-10)
})

test_that("kernel_roots handles near-zero eigenvalues with jitter", {
  K <- diag(c(1, 1e-12, 1e-15))
  roots <- kernel_roots(K, jitter = 1e-10)

  # All eigenvalues should be clamped to at least jitter
  expect_true(all(roots$evals >= 1e-10))
  expect_equal(roots$n_clamped, 2L)
})

test_that("kernel_roots warns for asymmetric input", {
  K <- matrix(c(1, 0.5, 0.6, 1), 2, 2)  # Asymmetric
  expect_warning(kernel_roots(K), "symmetric")
})

# ---------------------------------------------------------------------------
# Edge case tests: Circular factor
# ---------------------------------------------------------------------------

test_that("circular factor with L=2 is well-defined", {
  factors <- list(A = list(L = 2, type = "circular", l = 1.0))
  K_res <- design_kernel(factors, basis = "cell", normalize = "none")

  expect_equal(dim(K_res$K), c(2, 2))
  expect_true(isSymmetric(K_res$K, tol = 1e-10))

  # Circular distance between positions 1 and 2 with L=2: min(1, 1) = 1
  # So off-diagonal should be exp(-1/(2*1^2)) = exp(-0.5)
  # Note: diagonal has jitter and rho0 added, off-diagonal only has jitter
  # The actual kernel is constructed by Kronecker products, so check relative value
  expect_equal(K_res$K[1, 2], K_res$K[2, 1], tolerance = 1e-10)
})
