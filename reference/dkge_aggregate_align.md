# Align aggregate components to a reference fit

Aligns a resampled aggregate fit to a reference using sign alignment for
rank one and K-Procrustes alignment for rank greater than one.

## Usage

``` r
dkge_aggregate_align(reference, fit, rank = NULL, degeneracy_tol = 0.001)
```

## Arguments

- reference:

  Reference \`dkge_aggregate_fit\`.

- fit:

  Fit to align.

- rank:

  Optional rank truncation before alignment.

- degeneracy_tol:

  Relative singular-value gap below which components are flagged as
  near-tied. The gap is computed on the \*unrotated\* (sorted) spectrum
  of \`fit\`, which is the quantity that governs how identifiable the
  rotation is. Treat \`alignment_summary\$min_cosine\` as the primary
  diagnostic; \`near_tie\` is a coarse flag on this relative gap.

## Value

Aligned \`dkge_aggregate_fit\` with \`alignment\` metadata. \`U\`,
\`saliences\`, \`scores_feature\`, and \`singular_values\` are all
carried through the rotation together: \`singular_values\` is recomputed
as \`sqrt(colSums(scores_feature^2))\`, which reproduces the unaligned
singular values exactly when the rotation is the identity. \`Chat\` and
\`eig_values\` are rotation-invariant properties of the data and are
left untouched, so after alignment \`singular_values\` no longer equals
\`sqrt(eig_values\[seq_len(rank)\])\`.

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
dkge_aggregate_align(fit, fit)
#> <dkge_aggregate_fit>
#>   estimand: aggregate_cell_mean
#>   rows:     4
#>   rank:     1
#>   singular: 1.367
```
