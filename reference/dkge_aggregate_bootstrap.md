# Bootstrap aggregate targets over subjects

Resamples subjects with replacement, recomputes aggregate rows, refits
the aggregate decomposition, aligns components, and evaluates the
statistic.

## Usage

``` r
dkge_aggregate_bootstrap(
  target,
  K = NULL,
  statistic = c("singular_value", "salience_contrast", "component_score_contrast",
    "contrast_score", "between_group_contrast"),
  B = 999L,
  strata = NULL,
  rank = NULL,
  center = c("none", "grand", "row", "column"),
  component = 1L,
  conf = 0.95,
  interval = c("percentile", "basic"),
  component_contrasts = NULL,
  component_scale = c("score", "salience"),
  return_features = FALSE,
  parallel = FALSE,
  seed = NULL,
  ...
)

# S3 method for class 'dkge_aggregate_bootstrap'
print(x, ...)
```

## Arguments

- target:

  A \`dkge_aggregate_target\`.

- K:

  Optional aggregate-row design kernel. See \[dkge_aggregate_fit()\] for
  the row-space convention (\`q\` = number of aggregate rows).

- statistic:

  Built-in statistic name or function.

- B:

  Number of bootstrap resamples.

- strata:

  Optional subject-level strata for within-stratum bootstrap: either a
  vector with one entry per source subject, or the name of a column in
  the subject data. Defaults to the interaction of
  \`target\$group_vars\`. Explicit strata \*\*must\*\* nest within
  \`group_vars\` (each stratum lies inside a single group level) and are
  rejected up front otherwise: a coarser or cross-cutting stratification
  can draw a resample in which an entire group level is absent, which
  changes the aggregate row set and fails against the fixed observed
  kernel. Every draw is additionally checked against the observed row
  set.

- rank, center:

  Passed to \[dkge_aggregate_fit()\].

- component:

  Component index passed to \[dkge_aggregate_stat()\] for the built-in
  statistics. It is a formal argument rather than part of \`...\`
  because R's partial matching would otherwise bind a \`component =\`
  argument to \`component_scale\`/\`component_contrasts\`. User-supplied
  statistic functions do not receive it.

- conf:

  Confidence level for the interval. Quantiles use
  \`stats::quantile(..., type = 6)\` (Weibull plotting positions), which
  is slightly wider than the default type 7 in small samples.

- interval:

  \`"percentile"\` uses the bootstrap quantiles directly. \`"basic"\` is
  the reflection interval \\(2\hat\theta - q\_{1-\alpha/2},\\
  2\hat\theta - q\_{\alpha/2})\\, which recentres a biased top singular
  value around the observed statistic.

- component_contrasts:

  Optional aggregate-row contrast matrix used to summarize each aligned
  component in every bootstrap draw.

- component_scale:

  Whether component contrasts are applied to singular-value-scaled cell
  scores (\`"score"\`) or raw saliences (\`"salience"\`).

- return_features:

  Logical; if \`TRUE\`, accumulate aligned feature-space component maps
  across bootstrap draws and return streaming mean, SD, and
  bootstrap-ratio z maps (\`observed / bootstrap SD\`, the standard PLS
  bootstrap ratio).

- parallel:

  Logical; if \`TRUE\`, evaluate draws with
  \`future.apply::future_lapply()\`. Bootstrap indices are drawn up
  front in the calling RNG stream, so results are identical for either
  setting of \`parallel\` given the same \`seed\`.

- seed:

  Optional random seed.

- ...:

  Additional arguments passed to \[dkge_aggregate_stat()\]. Every other
  setting is a named formal, so a misspelled argument (\`sed = 1\` for
  \`seed = 1\`) lands here and is reported as an error whenever
  \`statistic\` is one of the built-ins; \`...\` reaches only
  user-supplied statistic functions.

- x:

  A \`dkge_aggregate_bootstrap\` object to print.

## Value

Object of class \`dkge_aggregate_bootstrap\`. \`observed\` and
\`statistics\` are signed. \`interval\` is the requested bootstrap
interval of the signed statistic. \`excludes_zero\` reports whether that
interval excludes zero, or \`NA\` for the non-negative
\`"singular_value"\` statistic (the comparison is uninformative there).

## Details

Resampling is \*\*stratified by default\*\*: when \`strata\` is \`NULL\`
and the target has subject-level \`group_vars\`, subjects are resampled
with replacement \*within\*
\`interaction(subject_data\[target\$group_vars\])\`, so every draw
preserves the observed group sizes and the aggregate row set. Pass
\`strata\` explicitly (a vector, or the name of a column in the subject
data) to stratify differently; a stratum with a single subject
contributes that subject to every draw, which is reported with a
warning.

## Examples

``` r
set.seed(1)
values <- lapply(1:4, function(i) {
  M <- matrix(rnorm(2 * 3), 2, 3)
  dimnames(M) <- list(c("c1", "c2"), paste0("f", 1:3))
  M
})
names(values) <- paste0("s", 1:4)
subject_data <- data.frame(subject_id = names(values),
                           grp = c("A", "A", "B", "B"))
target <- dkge_aggregate_target(values, subject_data, group_vars = "grp")
# \donttest{
dkge_aggregate_bootstrap(target, statistic = "singular_value",
                         B = 49, rank = 1, seed = 1, interval = "basic")
#> <dkge_aggregate_bootstrap>
#>   statistic: 1.367
#>   interval:  0.0515, 1.367
#>   B:         49
# }
```
