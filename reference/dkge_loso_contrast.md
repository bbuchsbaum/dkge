# Leave-one-subject-out DKGE contrast

Leave-one-subject-out DKGE contrast

## Usage

``` r
dkge_loso_contrast(fit, s, contrasts, ridge = 0)
```

## Arguments

- fit:

  \`dkge\` object

- s:

  Subject index (1-based)

- contrasts:

  Contrast vector in the original design basis

- ridge:

  Optional ridge when recomputing the held-out compressed matrix

## Value

List with fields \`v\`, \`alpha\`, and \`basis\`

## Examples

``` r
# \donttest{
toy <- dkge_sim_toy(
  factors = list(A = list(L = 2), B = list(L = 3)),
  active_terms = c("A", "B"), S = 4, P = 20, snr = 5
)
fit <- dkge_fit(toy$B_list, toy$X_list, toy$K, rank = 2)
c_vec <- c(1, -1, rep(0, 3))
result <- dkge_loso_contrast(fit, s = 1, contrasts = c_vec)
# }
```
