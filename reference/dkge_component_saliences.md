# DKGE component saliences in effect space

Extracts the design/effect-space saliences \\K U\\. These are the raw
coordinates most users should inspect first: each row is an effect or
design cell, each column is a DKGE latent variable.

## Usage

``` r
dkge_component_saliences(
  fit,
  comps = NULL,
  scale = c("raw", "unit", "zscore"),
  long = TRUE
)
```

## Arguments

- fit:

  Fitted \`dkge\` object.

- comps:

  Components to include (defaults to first min(rank, 6)).

- scale:

  Optional within-component display scaling. \`"raw"\` leaves saliences
  on their original scale, \`"unit"\` rescales each component to unit
  K-norm of the underlying latent vector (see Details), and \`"zscore"\`
  z-scores each component across effects.

- long:

  Logical; return a tidy long data frame when \`TRUE\`, otherwise a
  numeric effects-by-components matrix.

## Value

A data frame or matrix of component saliences.

## Details

The salience matrix is exactly \\K U\[, comps\]\\. Because \`fit\$U\` is
K-orthonormal (\\U^\top K U = I\\), each salience column already has
unit norm in the metric dual to \`K\`, i.e. \\s_j^\top K^{-1} s_j =
u_j^\top K u_j = 1\\. \`scale = "unit"\` therefore divides column \`j\`
by \\\sqrt{u_j^\top K u_j}\\ and is a no-op for a well-formed fit; it is
retained so that hand-assembled or perturbed bases are put back on a
comparable footing. It deliberately does \*\*not\*\* divide by the
Euclidean norm, which would distort the K geometry.

## Examples

``` r
toy <- dkge_sim_toy(
  factors = list(A = list(L = 2), B = list(L = 3)),
  active_terms = c("A", "B"), S = 3, P = 15, snr = 5
)
fit <- dkge(toy$B_list, toy$X_list, K = toy$K, rank = 2)
head(dkge_component_saliences(fit, comps = 1:2))
#>   effect component component_id      salience
#> 1      A       LV1            1  0.6022807855
#> 2     B1       LV1            1 -0.0105409164
#> 3     B2       LV1            1  0.3896172538
#> 4   A:B1       LV1            1 -0.0028486683
#> 5   A:B2       LV1            1 -0.0003776993
#> 6      A       LV2            2 -0.5506129871
dkge_component_saliences(fit, comps = 1, long = FALSE)
#>                LV1
#> A     0.6022807855
#> B1   -0.0105409164
#> B2    0.3896172538
#> A:B1 -0.0028486683
#> A:B2 -0.0003776993
```
