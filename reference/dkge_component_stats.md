# Component-level consensus statistics

Transports each subject's component loadings to a reference
parcellation, computes inference statistics, and returns tidy summaries
ready for visualisation.

## Usage

``` r
dkge_component_stats(
  fit,
  mapper = "sinkhorn",
  centroids = NULL,
  sizes = NULL,
  inference = "signflip",
  medoid = 1L,
  components = NULL,
  adjust = "fdr",
  ...
)

dkge_write_component_stats(fit, file, ...)
```

## Arguments

- fit:

  A fitted \`dkge\` object.

- mapper:

  Mapper strategy (string or \[dkge_mapper_spec()\]). Defaults to
  "sinkhorn".

- centroids:

  Optional list of subject centroid matrices; defaults to centroids
  stored in \`fit\` if available.

- sizes:

  Optional list of cluster masses (one vector per subject).

- inference:

  One of "signflip" or "parametric", or a list providing \`type\`,
  \`B\`, \`tail\`, and \`alpha\`.

- medoid:

  Reference subject index (defaults to 1).

- components:

  Optional vector of component indices or names; default is all
  components.

- adjust:

  Method supplied to \[stats::p.adjust()\] for multiple testing
  correction in the tidy summary.

- ...:

  Additional mapper-specific parameters (e.g. \`epsilon\`).

- file:

  Path to the CSV file where component statistics will be written.

## Value

A list with fields: - \`summary\`: tidy data frame of statistics and
p-values. - \`statistics\`: per-component statistic vectors. -
\`transport\`: per-component transported subject matrices.

## Examples

``` r
# \donttest{
toy <- dkge_sim_toy(
  factors = list(A = list(L = 2), B = list(L = 3)),
  active_terms = c("A", "B"), S = 3, P = 15, snr = 5
)
fit <- dkge(toy$B_list, toy$X_list, kernel = toy$K, rank = 2)
#> Warning: Argument 'kernel' is deprecated; use 'K' instead.
centroids <- lapply(toy$B_list, function(B) matrix(rnorm(ncol(B) * 3), ncol(B), 3))
res <- dkge_component_stats(fit,
                            centroids = centroids,
                            mapper = "ridge",
                            inference = "parametric",
                            components = 1)
head(res$summary)
#>   component cluster       stat           p      p_adj significant
#> 1         1       1   2.475140 0.131733897 0.14114346       FALSE
#> 2         1       2   2.061469 0.175389568 0.17538957       FALSE
#> 3         1       3   4.964777 0.038256630 0.05738495       FALSE
#> 4         1       4   6.560291 0.022455895 0.03750143        TRUE
#> 5         1       5 -15.980768 0.003892808 0.03750143        TRUE
#> 6         1       6  11.387447 0.007623576 0.03750143        TRUE
# }
```
