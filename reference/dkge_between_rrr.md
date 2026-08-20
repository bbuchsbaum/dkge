# Reduced-rank regression on a DKGE subject target

Fits a directional multivariate second-level model \$\$Y = X B + E,\quad
rank(B) \<= r\$\$ using a compact SVD of the fitted design-space
response. This is intended for subject-level factors, traits, and
covariates after DKGE has produced a subject-by-feature target.

## Usage

``` r
dkge_between_rrr(
  target,
  design,
  rank = NULL,
  weights = c("none", "target"),
  feature_mask = NULL,
  tol = 1e-08
)
```

## Arguments

- target:

  A \`dkge_target\` from \[dkge_make_target()\] or a numeric matrix.

- design:

  A \`dkge_subject_model\` from \[dkge_subject_model()\] or a numeric
  model matrix.

- rank:

  Reduced rank. Defaults to the maximum estimable rank.

- weights:

  \`"none"\`, \`"target"\`, or a list with optional \`subject\` and
  \`feature\` numeric weights.

- feature_mask:

  Optional logical vector selecting target features.

- tol:

  Numerical tolerance for rank checks; a single positive finite number.
  It is stored on the fit and reused for every refit performed by
  \[dkge_between_permute()\], so a design accepted here is also accepted
  by the permutation machinery.

## Value

Object of class \`dkge_between_rrr\`.

## Examples

``` r
set.seed(1)
n <- 12
dat <- data.frame(
  subject_id = paste0("s", seq_len(n)),
  group = factor(rep(c("A", "B"), length.out = n)),
  trait = scale(rnorm(n), center = TRUE, scale = FALSE)[, 1]
)
design <- dkge_subject_model(~ group * trait, dat)
Y <- design$X %*% matrix(c(0, 1, 0.5, -0.2), ncol(design$X), 1)
target <- dkge_make_target(Y = Y, subject_ids = dat$subject_id)
dkge_between_rrr(target, design, rank = 1)
#> <dkge_between_rrr>
#>   subjects : 12 
#>   features : 1 
#>   rank     : 1 
#>   terms    : (Intercept), group, trait, group:trait 
```
