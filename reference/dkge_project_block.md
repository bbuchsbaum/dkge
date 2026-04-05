# Project DKGE data into component space

Functions for projecting subject-standardised betas or new blocks into
DKGE component coordinates.

Convenience wrapper for projecting row-standardised betas onto DKGE
components.

## Usage

``` r
dkge_project_btil(fit, Btil)

dkge_project_block(
  fit,
  s,
  B_s,
  Omega_s = NULL,
  w_s = NULL,
  least_squares = TRUE
)
```

## Arguments

- fit:

  A \`dkge\` object.

- Btil:

  Either a qxP matrix or a list of such matrices (e.g. \`fit\$Btil\`).

- s:

  Block index (subject) to project against.

- B_s:

  Beta matrix for the new data block.

- Omega_s:

  Optional weights (vector length P or PxP matrix) matching the columns
  of \`B_s\`.

- w_s:

  Optional subject-level weight (defaults to 1 when omitted).

- least_squares:

  Logical; pass to \[multivarious::project_block()\].

## Value

List of Pxrank matrices; returns a single matrix when \`Btil\` is a
matrix.

Projection scores (qxrank) restricted to block \`s\`.

## Functions

- `dkge_project_btil()`: Project subject-standardised betas into
  component space

## Examples

``` r
toy <- dkge_sim_toy(
  factors = list(A = list(L = 2), B = list(L = 3)),
  active_terms = c("A", "B"), S = 3, P = 20, snr = 5
)
fit <- dkge_fit(toy$B_list, toy$X_list, toy$K, rank = 2)
A <- dkge_project_btil(fit, fit$Btil[[1]])
dim(A)
#> [1] 20  2
```
