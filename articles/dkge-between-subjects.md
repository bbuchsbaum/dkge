# Between-Subject Multivariate Modeling

Suppose you have one brain-derived feature vector per subject and want
to ask whether group, trait, or covariate effects are distributed across
that vector. The between-subject DKGE layer fits one supervised
multivariate model and returns a coefficient map and a resampling test
for each model term. The features may be component scores, transported
parcels, or any common-space subject-by-feature matrix.

The workflow has three steps:

1.  build a subject-by-feature target with
    [`dkge_make_target()`](https://bbuchsbaum.github.io/dkge/reference/dkge_make_target.md)
2.  build a second-level design with
    [`dkge_subject_model()`](https://bbuchsbaum.github.io/dkge/reference/dkge_subject_model.md)
3.  fit and test a reduced-rank model with
    [`dkge_between_rrr()`](https://bbuchsbaum.github.io/dkge/reference/dkge_between_rrr.md)
    and
    [`dkge_between_permute()`](https://bbuchsbaum.github.io/dkge/reference/dkge_between_permute.md)

## What does the second-level model estimate?

The between-subject model is a reduced-rank regression:

``` math
Y = X B + E, \qquad \mathrm{rank}(B) \le r
```

Here `Y` is a subject-by-feature matrix, `X` is the subject-level design
matrix, and `B` contains distributed coefficient maps for your terms of
interest.

## How do you fit a direct subject-by-feature matrix?

In the simplest case you already have a common-space matrix with one row
per subject.

``` r

subject_data <- data.frame(
  subject_id = paste0("s", 1:18),
  group = factor(rep(c("control", "patient"), each = 9)),
  trait = seq(-1, 1, length.out = 18)
)

design <- dkge_subject_model(~ group * trait, subject_data)
X <- design$X

beta <- matrix(0, nrow = ncol(X), ncol = 6,
               dimnames = list(colnames(X), paste0("feature", 1:6)))
beta["grouppatient:trait", ] <- c(2, 1, 0, -1, -2, -1)

Y <- X %*% beta + matrix(rnorm(nrow(X) * ncol(beta), sd = 0.08), nrow(X), ncol(beta))
rownames(Y) <- subject_data$subject_id
colnames(Y) <- colnames(beta)

target <- dkge_make_target(Y = Y)
fit_between <- dkge_between_rrr(target, design, rank = 1)
```

You can inspect the fitted term map directly.

``` r

round(dkge_term_map(fit_between, "grouppatient:trait"), 2)
#> feature1 feature2 feature3 feature4 feature5 feature6 
#>     1.94     1.02    -0.01    -0.99    -1.89    -0.96
```

For an unweighted, single-block model with a scientifically defensible
Gaussian row-sphericity assumption, request a residual-space rotation
test. The `terms` argument controls which formula terms are tested; it
does not refit a different full model for each term.

``` r

perm_between <- dkge_between_permute(
  fit_between,
  terms = c("group", "trait", "group:trait"),
  method = "rotation",
  B = 99,
  seed = 42,
  scope = "both",
  feature_adjust = "maxT"
)

perm_between$summary[, c("term", "statistic", "p", "p_adjusted")]
#>                    term    statistic    p p_adjusted
#> group             group 5.040060e-07 0.99       0.99
#> trait             trait 2.163909e-03 0.63       0.63
#> group:trait group:trait 4.263363e+00 0.01       0.01
```

The global statistic is the reduction in residual sum of squares when
the tested term is restored to the reduced model. Both methods refit at
the same fitted rank (`object$rank`). Freedman–Lane has an extra caveat:
if dropping the tested term lowers the reduced model’s estimable rank
below that value, the reduced fit is clipped while the full fit is not,
so the statistic also absorbs that rank difference. Rotation does not
clip — the reduced design stays in `Q0` and the Haar draw lives in its
complement. P-values use the finite-resample correction
`(1 + exceedances) / (B + 1)`. With `scope = "both"`, the result also
contains featurewise tests; `feature_adjust = "maxT"` controls
multiplicity against the largest feature statistic in each resample.

## What is residual-space rotation doing?

For a tested term, let $`X_0`$ be the design with that term removed.
Write $`Q_0`$ for an orthonormal basis of the nuisance space spanned by
$`X_0`$, and $`Q_e`$ for its orthogonal complement. A null replicate is

``` math
Y^* = Q_0 Q_0^\top Y + Q_e G Q_e^\top Y,
```

where $`G`$ is a random orthogonal matrix. The first term fixes the
fitted nuisance signal. The second rotates only the residual
coordinates. Because the transformation is orthogonal, it preserves
`crossprod(Y)` exactly and retains arbitrary covariance among features.

This produces an exact finite-sample randomization test under a
matrix-normal null with independent, homoscedastic subject rows. It is
only approximate for heteroscedastic, clustered, or non-Gaussian subject
errors. The result records that contract rather than asking you to
remember it:

``` r

unlist(perm_between$resampling_assumptions)
#>                                                                                       exact_under 
#> "matrix-normal null with independent homoscedastic subject rows and arbitrary feature covariance" 
#>                                                                                           weights 
#>                                                                                            "none" 
#>                                                                            exchangeability_blocks 
#>                                                                                          "single" 
#>                                                                                      non_gaussian 
#>                                                                                     "approximate"
```

The implementation evaluates each replicate from compressed
design-by-feature cross-products. It forms a square rotation only in the
residual subject dimension and never constructs a feature-by-feature
covariance matrix.

## Which resampling method should you choose?

Neither method is a universal confirmatory default.

| Situation | Method | Interpretation |
|----|----|----|
| Unweighted fit, one block, Gaussian row-sphericity is plausible | `method = "rotation"` | Sharpest exactness contract; may lose power in small samples. |
| Subject/feature weights or nontrivial exchangeability blocks are required | `method = "freedman_lane"` | Supported legacy approximation; small-sample calibration is imperfect. |
| Heteroscedasticity or dependence is scientifically central | Neither automatically | Treat the result as exploratory or justify a design-specific resampling model. |

`"freedman_lane"` remains the compatibility default, so use `method =`
in analysis scripts to make the inferential choice visible. Missing
block labels (`NA` in `blocks`) are rejected for both methods. Rotation
additionally fails before resampling when weights, multiple blocks, or
fewer than two residual dimensions put the fit outside its supported
scope.

### What did the frozen calibration show?

The tables below are generated from the machine-readable artifacts
shipped with the package. The null rows report rejection frequency at
nominal 0.05 for 20 subjects and 399 resamples.

| Term | Rotation: global null | Rotation: nuisance null | Freedman-Lane: global null | Rotation gate |
|:---|---:|---:|---:|:---|
| group | 0.049 | 0.040 | 0.074 | pass |
| trait | 0.055 | 0.042 | 0.076 | pass |
| group:trait | 0.058 | 0.066 | 0.077 | fail |

Rotation removed the material global-null inflation seen for
Freedman-Lane, but the nuisance-interaction arm rejected 33 of 500 nulls
(0.066), one result above the frozen 0.065 promotion ceiling. Its Wilson
interval still contained 0.05 and its other uniformity diagnostics
passed. The predeclared rule did not permit replacing that seed family
or relaxing the threshold.

| Method        | Detection rate | Wilson lower | Wilson upper |
|:--------------|---------------:|-------------:|-------------:|
| freedman_lane |           0.58 |        0.523 |        0.634 |
| rotation      |           0.43 |        0.375 |        0.487 |

In the frozen strong-interaction simulation (shared DGP, power arm
`B = 199`), rotation detected 43% and Freedman-Lane 58%. That raw gap is
mostly Freedman-Lane size inflation: at matched conditions rotation’s
size is 0.050 and Freedman-Lane’s is 0.086, and size-matched
Freedman-Lane power is about 0.49. The predeclared absolute and relative
gates still fail, so rotation remains opt-in. Increasing `B` improves
Monte Carlo resolution; it does not repair either method’s assumptions
or recover the unadjusted power gap.

## How do you build the target from a DKGE fit?

The more DKGE-native path is to let DKGE construct the subject
representation first, then model subjects in that common target space.

This example starts from heterogeneous subject maps, fits DKGE,
transports one contrast to a medoid subject, and then applies the
between-subject model.

``` r

toy <- dkge_sim_toy(
  factors = list(condition = list(L = 2), phase = list(L = 2)),
  active_terms = c("condition", "phase"),
  S = 6,
  P = c(4, 5, 3, 6, 4, 5),
  snr = 6,
  seed = 9
)

effects <- toy$info$term_names
toy$B_list <- lapply(toy$B_list, function(B) {
  rownames(B) <- effects
  B
})
toy$X_list <- lapply(toy$X_list, function(X) {
  colnames(X) <- effects
  X
})
centroids <- lapply(toy$B_list, function(B) matrix(rnorm(ncol(B) * 3), ncol(B), 3))

fit_dkge <- dkge_fit(
  dkge_data(toy$B_list, designs = toy$X_list, subject_ids = toy$subject_ids),
  K = toy$K,
  rank = 2
)
```

``` r

transport <- dkge_transport_spec(
  centroids = centroids,
  medoid = 1L,
  method = "sinkhorn",
  epsilon = 0.1
)

contrast <- matrix(c(1, 0, 0), ncol = 1, dimnames = list(NULL, "condition"))
rownames(contrast) <- effects

transported_target <- dkge_make_target(
  fit_dkge,
  type = "transported_maps",
  contrast = contrast,
  transport = transport,
  crossfit = "analytic"
)

dim(transported_target$Y)
#> [1] 6 4
```

Now the second-level model looks the same as before.

``` r

subject_meta <- data.frame(
  subject_id = toy$subject_ids,
  group = factor(c("control", "control", "control", "patient", "patient", "patient")),
  trait = c(-0.8, -0.3, 0.2, -0.1, 0.5, 1.0)
)

between_design <- dkge_subject_model(~ group * trait, subject_meta)
fit_transport <- dkge_between_rrr(transported_target, between_design, rank = 1)
perm_transport <- dkge_between_permute(
  fit_transport,
  terms = c("group", "trait", "group:trait"),
  method = "rotation",
  B = 49,
  seed = 7
)

perm_transport$summary[, c("term", "p")]
#>                    term    p
#> group             group 0.76
#> trait             trait 0.80
#> group:trait group:trait 1.00
```

## What should you inspect after fitting?

There are three outputs worth checking first:

- `coef(fit_transport)` gives the full coefficient matrix in target
  space
- `dkge_term_map(fit_transport, "group:trait")` extracts the
  term-specific map
- `perm_transport$summary` gives a global multivariate resampling test
  for each term

If you asked for featurewise inference,
`perm_between$feature_tests[["group:trait"]]` adds per-feature
statistics and adjusted p-values in the same target space.

Report the target construction, fitted rank, tested term, resampling
method, number of resamples, exchangeability assumptions, and
multiplicity adjustment. Do not describe a term-level p-value as
evidence that every feature in its map is individually significant.

## Where to go next

Use `type = "component_scores"` in
[`dkge_make_target()`](https://bbuchsbaum.github.io/dkge/reference/dkge_make_target.md)
when you want a compact latent target, and use
`type = "transported_maps"` when spatial interpretation in a common map
space matters. See
[`vignette("dkge-contrasts-inference")`](https://bbuchsbaum.github.io/dkge/articles/dkge-contrasts-inference.md)
for within-subject contrast inference. For cohorts with partial or
unbalanced effect grids, first see
[`vignette("dkge-partial-effect-spaces")`](https://bbuchsbaum.github.io/dkge/articles/dkge-partial-effect-spaces.md)
for the coverage, precision-weighting, and estimability contracts that
precede this between-subject layer.
