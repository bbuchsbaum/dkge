# One standard-error rule selection helper

Aggregates cross-validation scores by parameter setting and returns both
the best-performing parameter and the one within one standard error of
the best.

## Usage

``` r
dkge_one_se(scores, param_col = "param", metric_col = "score")
```

## Arguments

- scores:

  Data frame containing per-fold scores.

- param_col:

  Column name identifying the tuning parameter.

- metric_col:

  Column name holding the metric (larger is better).

## Value

List with \`best\`, \`pick\`, and \`summary\` table of mean/se by
parameter.

## Examples

``` r
scores <- data.frame(param = rep(1:3, each = 2), score = c(0.2, 0.1, 0.25, 0.2, 0.24, 0.23))
dkge_one_se(scores, param_col = "param", metric_col = "score")$pick
#> [1] 3
```
