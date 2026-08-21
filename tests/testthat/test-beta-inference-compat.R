as.matrix.beta_inference_contrasts <- function(x, contrast, ...) {
  x$values[[contrast]]
}

beta_inference_contrasts <- function(n_targets = 1L, q = 1L) {
  values <- lapply(seq_len(n_targets), function(i) {
    matrix(seq_len(6L * q) + i, nrow = 6L, ncol = q,
           dimnames = list(paste0("s", 1:6), paste0("feature", 1:q)))
  })
  names(values) <- paste0("target", seq_len(n_targets))
  structure(
    list(
      contrasts = setNames(as.list(seq_len(n_targets)), names(values)),
      method = "loso",
      values = values
    ),
    class = c("beta_inference_contrasts", "list")
  )
}

test_that("the inference compatibility table fails closed", {
  expect_error(
    .dkge_validate_inference_compatibility(
      inference = "parametric", correction = "maxT",
      has_transport = FALSE, center = "mean", n_targets = 1L,
      parallel = FALSE
    ),
    "[Pp]arametric.*maxT",
    class = "dkge_inference_compatibility_error"
  )
  expect_error(
    .dkge_validate_inference_compatibility(
      inference = "parametric", correction = "none",
      has_transport = TRUE, center = "mean", n_targets = 1L,
      parallel = FALSE
    ),
    "[Pp]arametric.*transport",
    class = "dkge_inference_compatibility_error"
  )
  expect_error(
    .dkge_validate_inference_compatibility(
      inference = "freedman-lane", correction = "none",
      has_transport = FALSE, center = "mean", n_targets = 1L,
      parallel = FALSE
    ),
    "Freedman-Lane.*not implemented",
    class = "dkge_inference_compatibility_error"
  )
  expect_error(
    .dkge_validate_inference_compatibility(
      inference = "signflip", correction = "none",
      has_transport = FALSE, center = "median", n_targets = 1L,
      parallel = FALSE
    ),
    "center.*mean",
    class = "dkge_inference_compatibility_error"
  )
})

test_that("dkge_infer rejects unsupported combinations before fitting", {
  expect_error(
    dkge_infer(
      list(), 1, inference = "signflip", correction = "none", n_perm = 0
    ),
    "strictly positive integer",
    class = "dkge_validation_error"
  )
  expect_error(
    dkge_infer(list(), 1, inference = "parametric", correction = "maxT"),
    "[Pp]arametric.*maxT",
    class = "dkge_inference_compatibility_error"
  )
  expect_error(
    dkge_infer(
      list(), 1, inference = "parametric", correction = "none",
      transport = list(provenance = list(class = "prespecified_frozen"))
    ),
    "[Pp]arametric.*transport",
    class = "dkge_inference_compatibility_error"
  )
  expect_error(
    dkge_infer(list(), 1, inference = "freedman-lane", correction = "none"),
    "Freedman-Lane.*not implemented",
    class = "dkge_inference_compatibility_error"
  )
})

test_that("the full inference compatibility Cartesian table is explicit", {
  table <- expand.grid(
    inference = c("signflip", "parametric", "freedman-lane"),
    correction = c("maxT", "fdr", "bonferroni", "none"),
    has_transport = c(FALSE, TRUE),
    center = c("mean", "median"),
    n_targets = c(1L, 2L),
    parallel = c(FALSE, TRUE),
    stringsAsFactors = FALSE
  )
  table$allowed <- with(
    table,
    center == "mean" &
      (inference == "signflip" |
       (inference == "parametric" & correction != "maxT" & !has_transport))
  )

  for (i in seq_len(nrow(table))) {
    row <- table[i, ]
    call <- function() {
      .dkge_validate_inference_compatibility(
        inference = row$inference,
        correction = row$correction,
        has_transport = row$has_transport,
        transport_provenance = if (row$has_transport) {
          list(class = "prespecified_frozen")
        } else {
          NULL
        },
        center = row$center,
        n_targets = row$n_targets,
        parallel = row$parallel
      )
    }
    if (row$allowed) {
      expect_silent(call())
    } else {
      expect_error(
        call(), class = "dkge_inference_compatibility_error",
        info = paste(row, collapse = "/")
      )
    }
  }
})

test_that("one-target sign-flip inference keeps a named matrix schema", {
  contrasts <- beta_inference_contrasts(n_targets = 1L, q = 1L)
  set.seed(102)
  result <- .infer_signflip(
    contrasts, n_perm = 100, correction = "none",
    mapped_values = contrasts$values,
    center = "mean", parallel = FALSE
  )

  expect_named(result$statistics, "target1")
  expect_named(result$p_values, "target1")
  expect_named(result$p_adjusted, "target1")
  expect_named(result$statistics$target1, "feature1")
  expect_length(result$p_values$target1, 1L)
  expect_length(result$p_adjusted$target1, 1L)
})

test_that("serial and parallel inference have the same stable schema", {
  skip_if_not_installed("future")
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  future::plan(future::sequential)

  contrasts <- beta_inference_contrasts(n_targets = 2L, q = 2L)
  set.seed(913)
  serial <- .infer_signflip(
    contrasts, n_perm = 100, correction = "bonferroni",
    mapped_values = contrasts$values,
    center = "mean", parallel = FALSE
  )
  set.seed(913)
  parallel <- suppressWarnings(.infer_signflip(
    contrasts, n_perm = 100, correction = "bonferroni",
    mapped_values = contrasts$values,
    center = "mean", parallel = TRUE
  ))

  expect_identical(names(serial$statistics), names(parallel$statistics))
  expect_equal(serial$statistics, parallel$statistics, tolerance = 0)
  expect_equal(serial$p_values, parallel$p_values, tolerance = 0)
  expect_equal(serial$p_adjusted, parallel$p_adjusted, tolerance = 0)
})

test_that("public sign-flip center semantics reject unsupported statistics", {
  Y <- matrix(c(-4, -1, 0, 1, 9, 2, 3, 4, 5, 6), ncol = 2)
  expect_error(
    dkge_signflip_maxT(Y, B = 100, center = "median"),
    "center.*mean",
    class = "dkge_inference_compatibility_error"
  )
})
