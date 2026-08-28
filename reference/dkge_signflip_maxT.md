# One-sample sign-flip max-T inference on transported subject maps

Computes cluster-wise one-sample t-statistics across subjects on
transported values (SxQ matrix), and calibrates p-values by the
max-\|t\| distribution under random subject-wise sign flips (symmetric
null). The helper conditions on the supplied matrix and is exact only
when that complete matrix is jointly row-sign invariant. It does not
establish that an upstream adaptive DKGE estimator has this property.

## Usage

``` r
dkge_signflip_maxT(
  Y,
  B = 2000,
  center = c("mean", "median"),
  tail = c("two.sided", "greater", "less")
)
```

## Arguments

- Y:

  SxQ matrix of aligned subject values on one identified reference
  support (rows = subjects, columns = support locations).

- B:

  number of sign-flip permutations

- center:

  "mean" or "median" for the location statistic (t uses mean)

- tail:

  "two.sided" \| "greater" \| "less"

## Value

A list with fields: `stat` (Q-vector of observed t-statistics), `p`
(Q-vector of max-T family-wise-error adjusted p-values), `p_unadj`
(Q-vector of per-column unadjusted permutation p-values), `maxnull`
(B-vector of permutation maximum statistics), and `flips` (S-by-B sign
matrix). Statistic and p-value names follow `colnames(Y)` (or stable
`feature*` defaults); flip rows follow `rownames(Y)` (or `subject*`
defaults).
