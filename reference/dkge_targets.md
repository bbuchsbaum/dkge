# Build classification targets from a DKGE fit

Generates a list of target specifications that map cell-level design
structure to the DKGE effect space. Each target provides a weight matrix
that can be multiplied with subject effect coefficients to obtain
class-specific patterns for decoding.

## Usage

``` r
dkge_targets(
  fit,
  spec = NULL,
  residualize = TRUE,
  collapse = NULL,
  scope = "within_subject",
  restrict_factors = NULL
)
```

## Arguments

- fit:

  dkge object containing \`kernel_info\$map\` metadata.

- spec:

  Target specification. Accepts a formula (e.g. \`~ A + B + A:B\`), a
  character vector of term labels, the string "fullcell", or an existing
  list of \`dkge_target\` objects (returned unchanged).

- residualize:

  Logical; if \`TRUE\` (default) residualise higher-order targets
  against previously constructed lower-order targets.

- collapse:

  Optional named list describing how to collapse factors that do not
  appear in a target. Each entry may be \`"mean"\`, \`list(method =
  "mean", window = 3:8)\`, or a numeric vector of length equal to the
  number of levels providing custom weights (automatically normalised).

- scope:

  Permutation/exchangeability scope stored with each target (default
  "within_subject").

- restrict_factors:

  Optional character vector restricting the factors used when \`spec =
  "fullcell"\`. When \`NULL\`, all factors are used.

## Value

List of objects with class \`dkge_target\`. Each target contains
\`name\`, \`factors\`, \`labels\`, \`weight_matrix\`, \`indicator\`,
\`residualized\`, \`collapse\`, and \`scope\` fields.

## Examples

``` r
# \donttest{
toy <- dkge_sim_toy(
  factors = list(A = list(L = 2), B = list(L = 3)),
  active_terms = c("A", "B"), S = 3, P = 20, snr = 5
)
kern <- design_kernel(factors = list(A = list(L = 2), B = list(L = 3)), basis = "effect")
fit <- dkge(toy$B_list, toy$X_list, kernel = kern, rank = 2)
#> Warning: Argument 'kernel' is deprecated; use 'K' instead.
tg <- dkge_targets(fit, ~ A + B + A:B)
length(tg)
#> [1] 3
# }
```
