# Build a flexible design-similarity kernel

Constructs a PSD kernel that captures factorial similarity across design
effects, optionally mapping from cell space to effect space using
contrasts.

## Usage

``` r
design_kernel(
  factors,
  terms = NULL,
  rho = NULL,
  include_intercept = TRUE,
  rho0 = 1e-08,
  basis = c("cell", "effect"),
  contrasts = NULL,
  block_structure = NULL,
  normalize = c("unit_trace", "none", "unit_fro", "max_diag"),
  jitter = 1e-08
)

dkge_design_kernel(
  factors,
  terms = NULL,
  rho = NULL,
  include_intercept = TRUE,
  rho0 = 1e-08,
  basis = c("cell", "effect"),
  contrasts = NULL,
  block_structure = NULL,
  normalize = c("unit_trace", "none", "unit_fro", "max_diag"),
  jitter = 1e-08
)
```

## Arguments

- factors:

  Named list of factor specifications. Each factor is described by a
  list containing:

  type

  :   "nominal" \| "ordinal" \| "circular" \| "continuous" (default
      "nominal").

  L

  :   Number of levels (for discrete types).

  values

  :   Numeric coordinates for continuous factors.

  l

  :   Optional length-scale for ordinal/circular/continuous factors.

- terms:

  List of character vectors describing which factors appear in each
  kernel term (e.g. list("A","B", c("A","B"))). Defaults to all main
  effects plus the full interaction.

- rho:

  Named numeric weights per term (names like "A", "A:B"). Defaults to 1
  for each term if omitted. Must be non-negative.

- include_intercept:

  Logical; if TRUE adds a small identity ridge (controlled by \`rho0\`)
  to keep the kernel full rank (default TRUE).

- rho0:

  Non-negative scalar ridge weight added when \`include_intercept\` is
  TRUE (default 1e-8).

- basis:

  Either "cell" (kernel over all design cells) or "effect" (kernel over
  regressors/effects). Default "cell".

- contrasts:

  Optional named list of per-factor contrast matrices used when
  \`basis="effect"\`. Defaults to sum-to-zero contrasts for discrete
  factors and a single column of ones for continuous factors.

- block_structure:

  Optional ordering of effect blocks (names matching terms). If NULL,
  uses the order of \`terms\`.

- normalize:

  One of "unit_trace", "none", "unit_fro", "max_diag". Controls how the
  kernel is scaled after construction (default "unit_trace").

- jitter:

  Small diagonal jitter added to \`K_cell\` for numerical stability
  (default 1e-8).

## Value

A list with elements \`K\` (kernel in requested basis), \`K_cell\`
(always returned), and \`info\` containing metadata such as factor/term
names, mapping matrix, and block indices.

## Examples

``` r
# Simple 2x3 factorial design
kern <- design_kernel(
  factors = list(A = list(L = 2), B = list(L = 3)),
  basis = "effect"
)
dim(kern$K)  # 5x5 effect-space kernel
#> [1] 5 5

# Ordinal factor with RBF smoothing
kern_ord <- design_kernel(
  factors = list(time = list(L = 5, type = "ordinal", l = 1.5)),
  basis = "effect"
)
```
