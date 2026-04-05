# Apply helper with optional parallelism

Wraps \`lapply()\` with an optional future.apply backend so callers can
enable \`parallel = TRUE\` without repeating boilerplate dependency
checks.

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

  Logical; if \`TRUE\`, uses \`future.apply::future_lapply()\`.

- ...:

  Additional arguments passed to the apply backend.

## Value

List of results matching \`lapply()\` semantics.
