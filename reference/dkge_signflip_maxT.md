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
  center = "mean",
  tail = c("two.sided", "greater", "less"),
  flips = NULL
)
```

## Arguments

- Y:

  SxQ matrix of subject values on the medoid parcellation
  (rows=subjects, cols=clusters)

- B:

  number of sign-flip permutations

- center:

  Location statistic. The beta API supports only \`"mean"\`, matching
  the one-sample t statistic used for observed and randomized data.

- tail:

  "two.sided" \| "greater" \| "less"

- flips:

  Optional precomputed S-by-B matrix of -1/+1 signs. This is an advanced
  reproducibility hook used to make serial and parallel execution
  consume exactly the same randomization descriptors.

## Value

A list with fields: \`stat\` (Q-vector of observed t-statistics), \`p\`
(Q-vector of max-T family-wise-error adjusted p-values), \`p_unadj\`
(Q-vector of per-column unadjusted permutation p-values), \`maxnull\`
(B-vector of permutation maximum statistics), and \`flips\` (S-by-B sign
matrix). Statistic and p-value names follow \`colnames(Y)\` (or stable
\`feature\*\` defaults); flip rows follow \`rownames(Y)\` (or
\`subject\*\` defaults).
