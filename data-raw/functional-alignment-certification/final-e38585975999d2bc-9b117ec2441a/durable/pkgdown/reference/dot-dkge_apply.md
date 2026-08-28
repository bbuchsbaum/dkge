# Apply helper with optional parallelism

Wraps [`lapply()`](https://rdrr.io/r/base/lapply.html) with an optional
future.apply backend so callers can enable `parallel = TRUE` without
repeating boilerplate dependency checks.

## Usage

``` r
.dkge_apply(X, FUN, parallel = FALSE, ...)
```

## Arguments

- X:

  Vector or list to iterate over.

- FUN:

  Function to apply.

- parallel:

  Logical; if `TRUE`, uses
  [`future.apply::future_lapply()`](https://future.apply.futureverse.org/reference/future_lapply.html).

- ...:

  Additional arguments passed to the apply backend.

## Value

List of results matching
[`lapply()`](https://rdrr.io/r/base/lapply.html) semantics.
