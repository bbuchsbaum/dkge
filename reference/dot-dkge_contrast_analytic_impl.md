# Analytic cross-fitting for multiple contrasts

Internal implementation called by dkge_contrast() for method="analytic".
Uses first-order perturbation theory to approximate LOSO contrasts.

## Usage

``` r
.dkge_contrast_analytic_impl(
  fit,
  contrast_list,
  ridge,
  parallel,
  verbose,
  tol = 1e-06,
  fallback = TRUE,
  align = TRUE,
  ...
)
```

## Arguments

- fit:

  dkge object

- contrast_list:

  List of normalized contrasts

- ridge:

  Ridge parameter (unused in analytic, kept for consistency)

- parallel:

  Logical; enables future.apply-based parallelism for per-subject
  computations (requires future.apply)

- verbose:

  Print progress

- tol:

  Tolerance for perturbation stability

- fallback:

  Allow fallback to full eigen when unstable

- ...:

  Additional arguments

## Value

List with values, metadata, etc.
