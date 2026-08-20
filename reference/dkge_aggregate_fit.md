# Fit a DKGE-regularized aggregate decomposition

Fits a row-space decomposition of an aggregate target using only q x q
linear algebra, where q is the number of aggregate rows. This is
intended for cell-mean/PLS-bridge analyses, not subject-level DKGE
inference.

## Usage

``` r
dkge_aggregate_fit(
  target,
  K = NULL,
  rank = NULL,
  center = c("none", "grand", "row", "column")
)

# S3 method for class 'dkge_aggregate_fit'
print(x, ...)
```

## Arguments

- target:

  A \`dkge_aggregate_target\` or row-by-feature numeric matrix.

- K:

  Optional design kernel, or an object returned by \[design_kernel()\],
  for the \*aggregate row\* space. Note that \`q\` here is the number of
  aggregate rows (for example \`group:task:measure\` cells), not the
  number of subject-level GLM effects used elsewhere in the package:
  \`K\` must therefore be \`nrow(target\$Y)\` by \`nrow(target\$Y)\`.
  When \`K\` carries dimnames they are matched against the aggregate row
  IDs; a kernel with only row names (or only column names) is validated
  and reordered on whichever labels are present. Rank-deficient PSD
  kernels keep a true square root: null directions stay at zero rather
  than receiving the jitter that \[\`.dkge_kernel_roots()\`\] uses for
  invertibility elsewhere.

- rank:

  Number of components to retain. Requests larger than \`min(nrow(Y),
  ncol(Y))\` are capped with a message.

- center:

  Centering applied to the aggregate matrix before fitting. The default
  \`"none"\` keeps the grand mean in the decomposition, which makes the
  leading component largely a mean-level effect; under
  \[dkge_aggregate_permute()\] with \`statistic = "singular_value"\`
  that mean survives subject-label permutation, so the resulting test
  has little power. Use \`"grand"\`/\`"column"\` centering, or a
  contrast statistic, when the question is about differences between
  aggregate rows.

- x:

  A \`dkge_aggregate_fit\` object to print.

- ...:

  Unused; present for S3 method compatibility.

## Value

Object of class \`dkge_aggregate_fit\`. Notable fields:

- U:

  q x rank K-orthonormal row basis.

- saliences:

  `K %*% U`.

- scores_feature:

  p x rank feature-space scores, `t(Yc) %*% saliences`.

- singular_values:

  Per-component energies, `sqrt(colSums(scores_feature^2))`. For an
  unrotated fit these equal `sqrt(eig_values[seq_len(rank)])` exactly;
  after \[dkge_aggregate_align()\] they are the energies of the
  \*rotated\* components and are no longer sorted.

- eig_values:

  Full length-q eigenvalue spectrum of `Chat`.

- Chat:

  \\K^{1/2} Y_c Y_c^\top K^{1/2}\\. Both `Chat` and `eig_values`
  describe the data in the kernel metric and are invariant to the
  component rotation applied by \[dkge_aggregate_align()\]; only `U`,
  `saliences`, `scores_feature`, and `singular_values` are
  basis-dependent.

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
fit$singular_values
#>      LV1 
#> 1.367119 
```
