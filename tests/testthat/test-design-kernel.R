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

test_that("effect grid records scopes and cell labels", {
  grid <- dkge_effect_grid(
    factors = list(
      group = c("control", "sdam"),
      task = c("recog", "nback"),
      measure = c("low_mid", "high_sem")
    ),
    scope = c(group = "between", task = "within", measure = "within"),
    block_factors = "group"
  )

  expect_s3_class(grid, "dkge_effect_grid")
  expect_equal(grid$scope["group"], c(group = "between"))
  expect_equal(
    grid$cell_labels,
    c("control:recog:low_mid", "control:recog:high_sem",
      "control:nback:low_mid", "control:nback:high_sem",
      "sdam:recog:low_mid", "sdam:recog:high_sem",
      "sdam:nback:low_mid", "sdam:nback:high_sem")
  )
})

test_that("design_kernel carries scope metadata and cell dimnames", {
  grid <- dkge_effect_grid(
    factors = list(
      group = c("control", "sdam"),
      task = c("recog", "nback"),
      measure = c("low_mid", "high_sem")
    ),
    scope = c(group = "between", task = "within", measure = "within")
  )

  Kres <- design_kernel(
    grid,
    terms = list("group", "task", c("group", "task")),
    basis = "cell",
    normalize = "none",
    include_intercept = FALSE,
    jitter = 0
  )

  expect_equal(Kres$info$term_scope["group"], c(group = "between"))
  expect_equal(Kres$info$term_scope["task"], c(task = "within"))
  expect_equal(Kres$info$term_scope["group:task"], c("group:task" = "mixed"))
  expect_equal(rownames(Kres$K), grid$cell_labels)
  expect_equal(Kres$info$blocks$cells, seq_along(grid$cell_labels))
})

test_that("block factors prevent accidental cross-group coupling", {
  grid <- dkge_effect_grid(
    factors = list(
      group = c("control", "sdam"),
      task = c("recog", "nback"),
      measure = c("low_mid", "high_sem")
    ),
    scope = c(group = "between", task = "within", measure = "within"),
    block_factors = "group"
  )

  terms <- list("task", "measure", c("task", "measure"))
  K_block <- design_kernel(
    grid,
    terms = terms,
    basis = "cell",
    normalize = "none",
    include_intercept = FALSE,
    jitter = 0
  )$K
  K_coupled <- design_kernel(
    grid$factors,
    terms = terms,
    basis = "cell",
    normalize = "none",
    include_intercept = FALSE,
    jitter = 0
  )$K

  expect_true(all(K_block[1:4, 5:8] == 0))
  expect_true(any(K_coupled[1:4, 5:8] > 0))
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

# -------------------------------------------------------------------------
# Specification validation for design_kernel() and dkge_effect_grid()
# -------------------------------------------------------------------------

test_that("design_kernel rejects terms naming unknown factors", {
  factors <- list(A = list(L = 2), B = list(L = 3))
  expect_error(
    design_kernel(factors, terms = list("A", "Bee")),
    "unknown factor"
  )
  expect_error(
    design_kernel(factors, terms = list("A", c("A", "A"))),
    "repeats a factor"
  )
  expect_error(
    design_kernel(factors, terms = list("A", "B", "A")),
    "Duplicate kernel terms"
  )
})

test_that("design_kernel rejects rho names that are not term names", {
  factors <- list(A = list(L = 2), B = list(L = 3))
  expect_error(
    design_kernel(factors, terms = list("A", "B"), rho = c(A = 1, Bee = 2)),
    "Unknown: Bee"
  )
})

test_that("design_kernel validates level labels against L", {
  expect_error(
    design_kernel(list(A = list(L = 3, levels = c("a", "b")))),
    "`levels` has 2 entries but `L` is 3"
  )
  # levels alone is enough to determine L
  kern <- design_kernel(list(A = list(levels = c("lo", "hi"))), basis = "cell")
  expect_equal(rownames(kern$K), c("lo", "hi"))
})

test_that("dkge_effect_grid handles list specs without L", {
  grid <- dkge_effect_grid(list(
    time = list(type = "continuous", values = c(0, 1, 2))
  ))
  expect_equal(grid$factors$time$L, 3L)
  expect_equal(grid$factors$time$type, "continuous")
  expect_equal(grid$cell_labels, c("0", "1", "2"))

  grid2 <- dkge_effect_grid(list(A = list(levels = c("x", "y"))))
  expect_equal(grid2$factors$A$L, 2L)
  expect_equal(grid2$cell_labels, c("x", "y"))

  expect_error(
    dkge_effect_grid(list(A = list(type = "nominal"))),
    "must provide `L`"
  )
  expect_error(
    dkge_effect_grid(list(A = list(type = "continuous"))),
    "must provide `values`"
  )
  expect_error(
    dkge_effect_grid(list(A = list(L = 2, levels = c("a", "b", "c")))),
    "`levels` has 3 entries but `L` is 2"
  )
})

test_that("dkge_effect_grid and design_kernel produce the same cell order", {
  grid <- dkge_effect_grid(list(A = c("a1", "a2"), B = c("b1", "b2", "b3")))
  kern <- design_kernel(grid, basis = "cell")
  expect_equal(rownames(kern$K), grid$cell_labels)
  expect_equal(
    grid$cell_labels,
    dkge:::.dkge_effect_grid_cell_labels(grid$factors)$labels
  )
  # First factor varies slowest.
  expect_equal(grid$cell_labels,
               c("a1:b1", "a1:b2", "a1:b3", "a2:b1", "a2:b2", "a2:b3"))
})

test_that("effect grid order can be pinned onto dkge_data", {
  grid <- dkge_effect_grid(list(A = c("a1", "a2"), B = c("b1", "b2")))
  labels <- grid$cell_labels
  design <- diag(length(labels))
  colnames(design) <- labels

  # Deliberately present subjects in an order whose first-appearance union
  # disagrees with the grid order.
  subj <- function(id, keep) {
    set.seed(nchar(id) + which(c("s1", "s2") == id))
    B <- matrix(rnorm(length(labels) * 6), length(labels), 6,
                dimnames = list(labels, NULL))
    B <- B[keep, , drop = FALSE]
    dkge_subject(B, design = design[, keep, drop = FALSE], id = id)
  }
  data_default <- dkge_data(list(subj("s1", c(4L, 3L)), subj("s2", 1:4)))
  expect_false(identical(data_default$effects, labels))

  data_grid <- dkge_data(list(subj("s1", c(4L, 3L)), subj("s2", 1:4)),
                         effects = labels)
  expect_equal(data_grid$effects, labels)

  kern <- design_kernel(grid, basis = "cell")
  fit <- dkge_fit(data_grid, K = kern$K, rank = 2, w_method = "none")
  expect_equal(fit$effects, labels)
  expect_equal(colnames(fit$K), labels)
})

test_that("ordinal length-scale is not shadowed by the `levels` element", {
  # `f$l` partially matches `f$levels`, and every dkge_effect_grid() spec
  # carries `levels`, so the RBF length-scale must be read exactly.
  grid <- dkge_effect_grid(list(load = list(L = 4, type = "ordinal", l = 1)))
  kern <- design_kernel(grid, basis = "cell", normalize = "none")
  Kc <- kern$K_cell
  # Similarity must decay with distance between ordinal levels.
  expect_gt(Kc[1, 2], Kc[1, 3])
  expect_gt(Kc[1, 3], Kc[1, 4])

  # A longer length-scale couples distant levels more strongly.
  grid_wide <- dkge_effect_grid(list(load = list(L = 4, type = "ordinal", l = 4)))
  Kw <- design_kernel(grid_wide, basis = "cell", normalize = "none")$K_cell
  expect_gt(Kw[1, 4], Kc[1, 4])

  # Same kernel whether the spec arrives with or without level labels.
  bare <- design_kernel(list(load = list(L = 4, type = "ordinal", l = 1)),
                        basis = "cell", normalize = "none")
  expect_equal(unname(Kc), unname(bare$K_cell), tolerance = 1e-12)

  expect_error(
    design_kernel(list(a = list(L = 3, type = "ordinal", l = "wide")),
                  basis = "cell"),
    "positive finite scalar"
  )
})

test_that("partially named rho scales only the terms it names", {
  factors <- list(A = list(L = 2), B = list(L = 3))
  terms <- list("A", "B", c("A", "B"))
  cell <- function(rho) {
    design_kernel(factors, terms = terms, rho = rho, basis = "cell",
                  normalize = "none", include_intercept = FALSE, jitter = 0)$K
  }

  # A partially named `rho` used to abort with "subscript out of bounds";
  # unnamed terms must simply keep their default weight of 1.
  K_full <- cell(c(A = 1, B = 1, "A:B" = 0.4))
  expect_equal(cell(c("A:B" = 0.4)), K_full, tolerance = 1e-12)
  # Order of the names is irrelevant.
  expect_equal(cell(c("A:B" = 0.4, B = 1, A = 1)), K_full, tolerance = 1e-12)
  # Unnamed rho stays positional.
  expect_equal(cell(c(1, 1, 0.4)), K_full, tolerance = 1e-12)

  expect_error(design_kernel(factors, terms = terms, rho = c(A = 1, A = 2)),
               "repeats term name")
  expect_error(design_kernel(factors, terms = terms, rho = c(1, 2)),
               "one entry per term")
  expect_error(design_kernel(factors, terms = terms, rho = c(A = NA_real_)),
               "finite")
  expect_error(design_kernel(factors, terms = terms, rho = list(A = 1)),
               "numeric")
})

test_that("factor specs reject vector-valued L and type cleanly", {
  # These used to reach `if (is.na(f$L) || f$L < 1L)` with a length-2 condition
  # and abort with R's "the condition has length > 1".
  expect_error(design_kernel(list(A = list(L = c(2, 3)))),
               "single level count")
  expect_error(dkge_effect_grid(list(A = list(L = c(2, 3)))),
               "single level count")
  expect_error(design_kernel(list(A = list(L = 2, type = c("nominal", "ordinal")))),
               "single character string")
  expect_error(dkge_effect_grid(list(A = list(L = 2, type = c("nominal", "ordinal")))),
               "single character string")
  # A length-2 `levels` against L = 3 still gets the specific message.
  expect_error(design_kernel(list(A = list(L = 3, levels = c("x", "y")))),
               "`levels` has 2 entries but `L` is 3")
  # A plain level vector remains valid.
  kern <- design_kernel(list(A = list(levels = c("x", "y"))), basis = "cell")
  expect_equal(rownames(kern$K), c("x", "y"))
})

test_that("duplicated level labels are rejected in both constructors", {
  expect_error(design_kernel(list(A = list(levels = c("x", "x")))),
               "duplicated level label")
  expect_error(design_kernel(list(A = list(L = 3, levels = c("x", "y", "x")))),
               "duplicated level label")
  expect_error(dkge_effect_grid(list(A = c("x", "x"))),
               "duplicated level label")
  expect_error(dkge_effect_grid(list(A = list(L = 2, levels = c("a", "a")))),
               "duplicated level label")
  expect_error(
    dkge_effect_grid(list(t = list(type = "continuous", values = c(0, 1, 0)))),
    "duplicated level label"
  )

  # Cell labels stay unique, so the kernel can still be indexed by name.
  grid <- dkge_effect_grid(list(A = c("a1", "a2"), B = c("b1", "b2")))
  expect_false(anyDuplicated(grid$cell_labels) > 0)

  # The historical single-factor default (terms = list("A", "A")) still repeats
  # *effect* names; that is deliberate and unchanged.
  K_eff <- design_kernel(list(A = list(L = 5)), basis = "effect")$K
  expect_equal(colnames(K_eff), rep(paste0("A", 1:4), 2))
})

test_that("continuous specs coming from an effect grid keep their numeric values", {
  # dkge_effect_grid() attaches `levels` to every spec, so factor-spec access
  # must use `[[` -- `$` would partially match `l` onto `levels`.
  spec <- list(t = list(type = "continuous", values = c(0, 1, 4)))
  grid <- dkge_effect_grid(spec)
  from_grid <- design_kernel(grid, basis = "cell", normalize = "none",
                             include_intercept = FALSE, jitter = 0)$K
  bare <- design_kernel(spec, basis = "cell", normalize = "none",
                        include_intercept = FALSE, jitter = 0)$K

  expect_equal(unname(from_grid), unname(bare), tolerance = 1e-12)
  expect_equal(rownames(from_grid), c("0", "1", "4"))
  # Similarity follows the numeric spacing, not the level index.
  expect_gt(from_grid[1, 2], from_grid[1, 3])
})

test_that("design_kernel and dkge_effect_grid reject unknown factor-spec fields", {
  # `typ` (typo for `type`) and `lev` (typo for `levels`) used to be silently
  # ignored, turning an ordinal factor into a nominal one.
  expect_error(
    design_kernel(factors = list(t = list(L = 4, typ = "ordinal", l = 1.5)),
                  terms = list("t"), basis = "cell"),
    "unrecognised field"
  )
  expect_error(
    design_kernel(factors = list(t = list(L = 3, lev = c("a", "b", "c"))),
                  terms = list("t"), basis = "cell"),
    "unrecognised field"
  )
  expect_error(
    design_kernel(factors = list(t = list(L = 3, type = "ordnal")),
                  terms = list("t"), basis = "cell"),
    "unsupported type"
  )
  expect_error(dkge_effect_grid(list(t = list(L = 3, lev = c("a", "b", "c")))),
               "unrecognised field")
  expect_error(
    dkge_effect_grid(list(t = list(L = 4, type = "ordinal", l = -1))),
    "length-scale"
  )
  expect_error(
    dkge_effect_grid(list(t = list(L = 4, type = "ordinal", l = c(1, 2)))),
    "length-scale"
  )
  # the correct spelling still works and is genuinely ordinal
  K <- design_kernel(factors = list(t = list(L = 4, type = "ordinal", l = 1.5)),
                     terms = list("t"), basis = "cell", normalize = "none",
                     include_intercept = FALSE, jitter = 0)$K
  expect_gt(K[1, 2], K[1, 3])
})

test_that("kernel_roots handles a 1x1 kernel without the diag() scalar footgun", {
  # diag(scalar) builds an identity of that dimension, so a 1x1 kernel silently
  # produced Khalf = 1 instead of sqrt(K). Passing the length restores correctness.
  K <- matrix(2, 1, 1)
  kr <- kernel_roots(K, jitter = 1e-12)
  expect_equal(dim(kr$Khalf), c(1L, 1L))
  expect_equal(as.numeric(kr$Khalf), sqrt(2), tolerance = 1e-8)
  expect_equal(as.numeric(kr$Kihalf), 1 / sqrt(2), tolerance = 1e-8)
  expect_equal(kr$Khalf %*% kr$Khalf, K, tolerance = 1e-8)
})
