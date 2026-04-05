# Fit DKGE with CPCA filtering

Convenience wrapper around \[dkge()\] that enables CPCA filtering with a
slightly terser interface. Either \`cpca_blocks\` or \`cpca_T\` must be
provided. Additional arguments are forwarded to \[dkge()\].

## Usage

``` r
dkge_cpca_fit(
  ...,
  cpca_blocks = NULL,
  cpca_T = NULL,
  cpca_part = c("design", "resid", "both"),
  cpca_ridge = 0
)
```

## Arguments

- ...:

  Additional arguments passed to \[dkge()\].

- cpca_blocks:

  Optional integer vector specifying the effect rows that span a CPCA
  design subspace. Ignored when \`cpca_part = "none"\` or when
  \`cpca_T\` is provided.

- cpca_T:

  Optional qxq0 matrix giving the CPCA design basis explicitly.
  Overrides \`cpca_blocks\` when supplied.

- cpca_part:

  Which CPCA-filtered component to fit: \`"none"\` (default) performs
  the standard DKGE fit; \`"design"\` uses only the CPCA design part;
  \`"resid"\` uses the residual part; \`"both"\` fits the design part
  but also stores the residual basis.

- cpca_ridge:

  Optional ridge applied to the CPCA-filtered matrices before
  eigen-decomposition.
