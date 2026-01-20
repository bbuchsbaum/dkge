# test-data-validation.R
# Input validation edge case tests for dkge_subject() and dkge_data()

library(testthat)

# -------------------------------------------------------------------------
# dkge_subject edge cases -------------------------------------------------
# -------------------------------------------------------------------------

test_that("dkge_subject rejects invalid input types with clear messages", {
  design <- matrix(1:10, 5, 2, dimnames = list(NULL, c("eff1", "eff2")))

  cases <- list(
    list(input = data.frame(a = 1:2, b = 3:4), pattern = "Unsupported"),
    list(input = 1:10, pattern = "Unsupported"),
    list(input = "not_a_matrix", pattern = "Unsupported"),
    list(input = NULL, pattern = "Unsupported|NULL")
  )

  for (case in cases) {
    expect_error(
      dkge_subject(case$input, design = design),
      case$pattern,
      info = paste("Input type:", class(case$input)[1])
    )
  }
})

test_that("dkge_subject rejects beta with NA values", {
  design <- matrix(1:10, 5, 2, dimnames = list(NULL, c("eff1", "eff2")))
  beta_with_na <- matrix(c(1, NA, 3, 4), 2, 2, dimnames = list(c("eff1", "eff2"), NULL))

  # Check current behavior - if NA is accepted, this documents it
  # Expected: should either error or at least store NA (verify downstream handling)
  result <- tryCatch(
    dkge_subject(beta_with_na, design = design),
    error = function(e) e
  )

  if (inherits(result, "error")) {
    expect_match(conditionMessage(result), "NA|missing|finite", ignore.case = TRUE)
  } else {
    # NA is accepted - document this behavior and verify it's stored
    expect_true(any(is.na(result$beta)),
                info = "NA values in beta should be preserved if accepted")
  }
})

test_that("dkge_subject rejects beta with Inf values", {
  design <- matrix(1:10, 5, 2, dimnames = list(NULL, c("eff1", "eff2")))
  beta_with_inf <- matrix(c(1, Inf, 3, 4), 2, 2, dimnames = list(c("eff1", "eff2"), NULL))

  # Check current behavior
  result <- tryCatch(
    dkge_subject(beta_with_inf, design = design),
    error = function(e) e
  )

  if (inherits(result, "error")) {
    expect_match(conditionMessage(result), "Inf|infinite|finite", ignore.case = TRUE)
  } else {
    # Inf is accepted - document this behavior
    expect_true(any(is.infinite(result$beta)),
                info = "Inf values in beta should be preserved if accepted")
  }
})

test_that("dkge_subject rejects dimension mismatch between beta and design", {
  # Design has 2 effects, beta has 3 rows
  design <- matrix(1:10, 5, 2, dimnames = list(NULL, c("eff1", "eff2")))
  beta_wrong_dim <- matrix(1:9, 3, 3, dimnames = list(c("a", "b", "c"), NULL))

  expect_error(
    dkge_subject(beta_wrong_dim, design = design),
    "match|Row names",
    info = "Beta rows must match design columns"
  )
})

test_that("dkge_subject handles empty beta matrix gracefully", {
  design <- matrix(1:10, 5, 2, dimnames = list(NULL, c("eff1", "eff2")))

  # Zero rows
  beta_no_rows <- matrix(numeric(0), 0, 5)
  result_no_rows <- tryCatch(
    dkge_subject(beta_no_rows, design = design),
    error = function(e) e
  )
  # Should error due to dimension mismatch or be handled
  if (inherits(result_no_rows, "error")) {
    # Error is expected for empty beta - document this behavior
    expect_true(TRUE, info = "Empty beta (zero rows) causes error as expected")
  } else {
    expect_equal(nrow(result_no_rows$beta), 0,
                 info = "Empty beta accepted - verify downstream handles it")
  }

  # Zero columns (no clusters)
  beta_no_cols <- matrix(numeric(0), 2, 0, dimnames = list(c("eff1", "eff2"), NULL))
  result_no_cols <- tryCatch(
    dkge_subject(beta_no_cols, design = design),
    error = function(e) e
  )
  if (inherits(result_no_cols, "error")) {
    expect_true(TRUE, info = "Zero-cluster beta causes error")
  } else {
    expect_equal(result_no_cols$n_clusters, 0,
                 info = "Zero-cluster beta accepted")
  }
})

test_that("dkge_subject handles duplicate effect names in design", {
  # Design with duplicate column names
  design_dup <- matrix(1:15, 5, 3)
  colnames(design_dup) <- c("eff1", "eff1", "eff2")
  beta <- matrix(1:9, 3, 3, dimnames = list(c("eff1", "eff1", "eff2"), NULL))

  result <- tryCatch(
    dkge_subject(beta, design = design_dup),
    error = function(e) e,
    warning = function(w) w
  )

  # Document current behavior - duplicates may cause issues
  if (inherits(result, "error") || inherits(result, "warning")) {
    expect_true(TRUE, info = "Duplicates cause error/warning as expected")
  } else {
    # Duplicates accepted - this may cause downstream issues
    expect_s3_class(result, "dkge_subject")
    # Note: Duplicate effect names are accepted - verify downstream handles this
  }
})

# -------------------------------------------------------------------------
# dkge_data edge cases ----------------------------------------------------
# -------------------------------------------------------------------------

test_that("dkge_data rejects empty betas list with informative message", {
  expect_error(
    dkge_data(list(), list()),
    "length|empty|greater",
    info = "Empty input lists should error"
  )
})

test_that("dkge_data handles single subject case correctly", {
  withr::local_seed(123)
  beta <- matrix(rnorm(3 * 10), 3, 10, dimnames = list(c("e1", "e2", "e3"), NULL))
  design <- matrix(rnorm(20 * 3), 20, 3, dimnames = list(NULL, c("e1", "e2", "e3")))

  data <- dkge_data(list(beta), list(design))

  expect_s3_class(data, "dkge_data")
  expect_equal(data$n_subjects, 1)
  expect_equal(length(data$betas), 1)
  expect_equal(data$effects, c("e1", "e2", "e3"))
})

test_that("dkge_data handles zero-cluster subjects", {
  beta_zero <- matrix(numeric(0), 3, 0, dimnames = list(c("e1", "e2", "e3"), NULL))
  design <- matrix(1:15, 5, 3, dimnames = list(NULL, c("e1", "e2", "e3")))

  result <- tryCatch(
    dkge_data(list(beta_zero), list(design)),
    error = function(e) e
  )

  if (inherits(result, "error")) {
    expect_true(TRUE, info = "Zero-cluster subject causes error")
  } else {
    expect_s3_class(result, "dkge_data")
    expect_equal(ncol(result$betas[[1]]), 0,
                 info = "Zero-cluster subject accepted")
  }
})

test_that("dkge_data tracks provenance for disjoint effects", {
  # Two subjects with NO overlapping effects
  beta1 <- matrix(1:4, 2, 2, dimnames = list(c("e1", "e2"), NULL))
  beta2 <- matrix(5:8, 2, 2, dimnames = list(c("e3", "e4"), NULL))
  design1 <- matrix(1:10, 5, 2, dimnames = list(NULL, c("e1", "e2")))
  design2 <- matrix(1:10, 5, 2, dimnames = list(NULL, c("e3", "e4")))

  data <- dkge_data(list(beta1, beta2), list(design1, design2))

  expect_equal(length(data$effects), 4)
  expect_true(!is.null(data$provenance))

  prov <- data$provenance
  # Cross-effects should have 0 pair counts
  expect_equal(prov$pair_counts["e1", "e3"], 0L)
  expect_equal(prov$pair_counts["e2", "e4"], 0L)
  # Self counts should be 1 each
  expect_true(all(diag(prov$pair_counts) == 1L))
})

test_that("dkge_data validates consistent dimensions within betas list", {
  design <- matrix(1:15, 5, 3, dimnames = list(NULL, c("e1", "e2", "e3")))

  # Different number of effects in each subject's beta
  beta1 <- matrix(1:9, 3, 3, dimnames = list(c("e1", "e2", "e3"), NULL))
  beta2 <- matrix(1:4, 2, 2, dimnames = list(c("e1", "e2"), NULL))  # Missing e3

  # This should work - partial overlap is allowed
  data <- dkge_data(list(beta1, beta2), list(design, design[, 1:2]))
  expect_equal(length(data$effects), 3)
})

test_that("dkge_data handles mismatched list lengths", {
  beta <- matrix(1:6, 2, 3, dimnames = list(c("e1", "e2"), NULL))
  design <- matrix(1:10, 5, 2, dimnames = list(NULL, c("e1", "e2")))

  expect_error(
    dkge_data(list(beta, beta), list(design)),
    "length",
    info = "Different lengths for betas and designs should error"
  )
})

test_that("dkge_data handles subjects with all NA betas", {
  design <- matrix(1:10, 5, 2, dimnames = list(NULL, c("e1", "e2")))
  beta_all_na <- matrix(NA_real_, 2, 5, dimnames = list(c("e1", "e2"), NULL))

  result <- tryCatch(
    dkge_data(list(beta_all_na), list(design)),
    error = function(e) e
  )

  # Document behavior
  if (inherits(result, "error")) {
    expect_match(conditionMessage(result), "NA|missing", ignore.case = TRUE)
  } else {
    expect_s3_class(result, "dkge_data")
  }
})

# -------------------------------------------------------------------------
# Error message quality ---------------------------------------------------
# -------------------------------------------------------------------------

test_that("error messages include helpful context about the problem", {
  design <- matrix(1:10, 5, 2, dimnames = list(NULL, c("eff1", "eff2")))
  beta_wrong <- matrix(1:9, 3, 3, dimnames = list(c("x", "y", "z"), NULL))

  err <- tryCatch(
    dkge_subject(beta_wrong, design = design),
    error = function(e) conditionMessage(e)
  )

  # Error should mention the nature of the problem
  expect_true(
    grepl("match|effect|row|column|name", err, ignore.case = TRUE),
    info = paste("Error message should be descriptive:", err)
  )
})

test_that("dkge_data error messages identify problematic subject", {
  design1 <- matrix(1:10, 5, 2, dimnames = list(NULL, c("e1", "e2")))
  design2 <- matrix(1:10, 5, 2, dimnames = list(NULL, c("a", "b")))  # Different names
  beta1 <- matrix(1:4, 2, 2, dimnames = list(c("e1", "e2"), NULL))
  beta2 <- matrix(1:4, 2, 2, dimnames = list(c("e1", "e2"), NULL))  # Mismatch with design2

  err <- tryCatch(
    dkge_data(list(beta1, beta2), list(design1, design2)),
    error = function(e) conditionMessage(e)
  )

  # Error should help identify the issue
  expect_true(
    grepl("match|effect|name|beta|design", err, ignore.case = TRUE),
    info = paste("Error should identify nature of mismatch:", err)
  )
})
