# Evaluate aggregate-fit statistics

Aggregate statistics operate on a cell-mean decomposition, not on a
subject-level second-level model. In particular,
\`"between_group_contrast"\` is a convenience bridge statistic: it
applies a supplied group-by-cell row contrast to the selected
component's singular-value-scaled aggregate saliences. It should not be
interpreted as a replacement for \[dkge_between_rrr()\] or
\[dkge_between_permute()\] on subject-level targets.

## Usage

``` r
dkge_aggregate_stat(
  fit,
  statistic = c("singular_value", "salience_contrast", "component_score_contrast",
    "contrast_score", "between_group_contrast"),
  component = 1L,
  contrast = NULL,
  ...
)
```

## Arguments

- fit:

  A \`dkge_aggregate_fit\`.

- statistic:

  Built-in statistic name or a function accepting \`fit\`.
  \`"singular_value"\` returns the selected component singular value.
  \`"salience_contrast"\` projects the selected component salience
  vector onto \`contrast\`. \`"component_score_contrast"\` applies the
  same contrast to the component's singular-value-scaled salience
  vector. \`"contrast_score"\` is a legacy alias for
  \`"salience_contrast"\`. \`"between_group_contrast"\` is a
  bridge-analysis alias for \`"component_score_contrast"\` and should be
  supplied a group-by-within row contrast.

- component:

  Component index.

- contrast:

  Numeric row contrast. Required for every statistic except
  \`"singular_value"\`, i.e. for \`"salience_contrast"\`,
  \`"component_score_contrast"\`, \`"contrast_score"\`, and
  \`"between_group_contrast"\`. May be named, in which case it is
  reordered to the aggregate row IDs.

- ...:

  Additional arguments passed to a user-supplied statistic function.
  Built-in statistics take only \`component\` and \`contrast\`, so any
  other argument reaching \`...\` is reported as an error rather than
  being silently ignored.

## Value

Numeric scalar statistic. Contrast statistics are returned
\*\*signed\*\*; the sign is only interpretable relative to the component
orientation of the fit the statistic is evaluated on.
\[dkge_aggregate_bootstrap()\] orients resampled fits with
\[dkge_aggregate_align()\] before taking the statistic;
\[dkge_aggregate_permute()\] does not (null draws stay unaligned;
two-sided \`abs()\` absorbs the arbitrary component sign).
\[dkge_aggregate_bootstrap()\] reports a percentile interval on the
signed statistic so that it can legitimately contain zero.

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
fit <- dkge_aggregate_fit(target, rank = 1)
dkge_aggregate_stat(fit, "singular_value")
#> [1] 1.367119
```
