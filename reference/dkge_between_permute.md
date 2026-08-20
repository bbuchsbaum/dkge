# Resampling tests for between-subject DKGE RRR terms

Performs term-specific global tests for a fitted \[dkge_between_rrr()\]
model. Residual-space rotation fixes the nuisance-design column space
and rotates its orthogonal complement before comparing the reduced and
full reduced-rank residual sums of squares. Freedman-Lane residual
permutation remains available as a legacy approximate method. The beta
API requires an explicit method because the frozen audit both found
material size inflation for the legacy method and failed its predeclared
promotion rule for rotation. Rotation is exact only under its documented
matrix-normal, homoscedastic, unweighted, single-block scope; other
error laws require cautious interpretation.

## Usage

``` r
dkge_between_permute(
  object,
  terms = NULL,
  method = NULL,
  B = 999L,
  blocks = NULL,
  seed = NULL,
  adjust = "none",
  statistic = c("frob"),
  scope = c("global", "features", "both"),
  feature_adjust = c("none", "fdr", "maxT"),
  parallel = FALSE
)
```

## Arguments

- object:

  A \`dkge_between_rrr\` object.

- terms:

  Character vector of terms or model-matrix columns to test. When
  \`NULL\`, all non-intercept formula terms except those listed in
  \`object\$design\$nuisance\` are tested.

- method:

  Required resampling method. \`"rotation"\` uses Haar rotations in the
  orthogonal complement of the reduced design; \`"freedman_lane"\`
  permutes reduced-model residual rows. There is deliberately no default
  in the beta API because neither method is universally certified.

- B:

  Number of rotations or permutations.

- blocks:

  Optional exchangeability blocks of length \`n_subjects\`.

- seed:

  Optional random seed.

- adjust:

  P-value adjustment method passed to \[stats::p.adjust()\].

- statistic:

  Test statistic. Currently \`"frob"\` for reduced-minus-full weighted
  residual sum of squares.

- scope:

  \`"global"\` for the term-level test only, \`"features"\` for
  featurewise tests only, or \`"both"\`.

- feature_adjust:

  Multiplicity adjustment for featurewise p-values.

- parallel:

  Logical; evaluate resampling replicates with
  \[future.apply::future_lapply()\]. Randomization descriptors are
  always drawn serially up front, so results are identical to the serial
  path for a given \`seed\`. The serial path streams replicates and is
  preferred when \`scope != "global"\` with many features.

## Value

Object of class \`dkge_between_permutation\`.

## Details

For \`method = "rotation"\`, let \`Q0\` span the reduced design and
\`Qe\` span its orthogonal complement. Each null response is equivalent
to \`Q0 Q0' Y + Qe G Qe' Y\`, where \`G\` is Haar-uniform on the
residual-space orthogonal group. This fixes every reduced-model mean and
preserves \`crossprod(Y)\` exactly. The test is finite-sample exact
under a matrix-normal null with independent homoscedastic subject rows
and arbitrary feature covariance. Version 1 supports unweighted fits and
one exchangeability block only; unsupported weights or blocks are
rejected. With non-Gaussian or heteroscedastic subject errors, rotation
is approximate rather than exact.

Both methods refit the reduced model at the same fitted rank
(\`object\$rank\`) as the full model. For \`method = "freedman_lane"\`,
when the tested term removes enough design columns that the reduced
model's estimable rank drops below that value, the reduced fit is
clipped to its own maximum rank while the full fit is not, and the
reduced-minus-full residual statistic then also absorbs that rank
difference. This affects the observed statistic and every permutation
replicate alike, but the null is not exactly a test of the omitted term
alone. Rotation does not have this caveat: the reduced design is held
fixed in \`Q0\` and the Haar draw lives in its orthogonal complement.

The legacy Freedman-Lane method's frozen package audit
(\`null-calibration-v1\`) found empirical type-I error of 0.074, 0.076,
and 0.077 at nominal 0.05 for three terms with 20 subjects (1,000 null
data sets per term and 399 permutations). Distortion attenuated but was
mixed at 40 subjects; at 80 subjects the three estimates were 0.030,
0.048, and 0.052. These results establish a material small-sample
limitation, not merely uncertainty from the number of permutations.
Treat small-sample p-values close to alpha as exploratory, report the
limitation, and use an independently justified exchangeability design.
Increasing \`B\` improves Monte Carlo resolution but does not remove
this size distortion. The frozen plan, replicate p-values, and summary
are shipped under \`data-raw/dkge-null-calibration-\*\` and
\`inst/extdata/dkge-null-calibration-\*\`.

A separate predeclared audit of residual rotation found global-null
sizes of 0.049, 0.055, and 0.058 and nonzero-nuisance null sizes of
0.040, 0.042, and 0.066 for the same three terms at 20 subjects. One
nuisance-interaction result missed the frozen promotion ceiling by one
rejection, and rotation's strong-interaction power was 0.430 versus
0.580 for Freedman-Lane. That raw gap is mostly Freedman-Lane size
inflation (matched-condition sizes 0.050 versus 0.086; size-matched
Freedman-Lane power about 0.49). The predeclared gates are not relaxed:
the missed nuisance-interaction ceiling and the power trade-off prevent
default promotion. The beta API therefore requires an explicit choice
instead of silently selecting either the demonstrably size-inflated
legacy method or a rotation method that missed its frozen promotion
rule. Its frozen plan, 8,100 results, and decision report are shipped as
\`data-raw/dkge-between-rotation-\*\` and
\`inst/extdata/dkge-between-rotation-\*\`.

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
fit <- dkge_between_rrr(target, design, rank = 1)
# \donttest{
dkge_between_permute(fit, terms = "group", method = "rotation",
                     B = 49, seed = 1)
#> <dkge_between_permutation>
#>   method       : rotation 
#>   resamples    : 49 
#>   scope        : global 
#>   term statistic    p  B   method statistic_name p_adjusted
#>  group  2.996028 0.02 49 rotation           frob       0.02
# }
```
