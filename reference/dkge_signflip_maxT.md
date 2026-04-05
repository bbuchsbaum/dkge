# One-sample sign-flip max-T inference on transported subject maps

Computes cluster-wise one-sample t-statistics across subjects on
transported values (SxQ matrix), and calibrates p-values by the
max-\|t\| distribution under random subject-wise sign flips (symmetric
null). This does not re-estimate DKGE, leveraging LOSO independence of
each subject's value.

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

  SxQ matrix of subject values on the medoid parcellation
  (rows=subjects, cols=clusters)

- B:

  number of sign-flip permutations

- center:

  "mean" or "median" for the location statistic (t uses mean)

- tail:

  "two.sided" \| "greater" \| "less"

## Value

list with fields: stat (Q-vector), p (Q-vector), maxnull (B-vector),
flips (SxB signs)
