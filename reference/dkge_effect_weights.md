# Specify effect-level reliability weighting

Effect weights control how unequally estimated subject-by-effect rows
enter the pooled second moment. They are distinct from
\[dkge_weights()\], which weights voxels or parcels, and from
\`w_method\`, which supplies one scalar per subject.

## Usage

``` r
dkge_effect_weights(
  method = c("none", "count", "precision", "random_effects"),
  within = c("auto", "noise", "count"),
  tau_method = "DL",
  max_ratio = 10,
  variance_floor = 1e-08,
  tau_floor = 0
)
```

## Arguments

- method:

  Effect precision source. \`"none"\` gives every observed effect equal
  precision, \`"count"\` uses \`effect_n\` stored on each subject,
  \`"precision"\` uses explicitly supplied \`effect_precision\` values,
  and \`"random_effects"\` combines within-subject uncertainty with a
  DerSimonian-Laird between-subject variance estimate.

- within:

  Within-subject variance source for \`method = "random_effects"\`.
  \`"noise"\` uses \`effect_noise_cov\` and residual variances,
  \`"count"\` uses \`1 / effect_n\`, and \`"auto"\` prefers complete
  noise statistics before falling back to complete counts.

- tau_method:

  Between-subject variance estimator. Currently only the non-negative
  DerSimonian-Laird moment estimator (\`"DL"\`) is supported.

- max_ratio:

  Maximum random-effects precision relative to the median positive
  precision for the same effect. This cap prevents a single subject from
  owning an effect geometry; \`Inf\` disables it.

- variance_floor:

  Positive lower bound for within-subject variances.

- tau_floor:

  Non-negative lower bound for the estimated between-subject variance.

## Value

An object of class \`dkge_effect_weights\`, suitable for the
\`effect_weights\` argument of \[dkge_fit()\]. A precision of zero
excludes that subject-effect row. Pairwise reliability is the geometric
mean of the two effect precisions.

## Details

For random effects, effect \`c\` uses \`p_sc = 1 / (v_sc + tau_c^2)\`.
The DerSimonian-Laird estimate is \`max(tau_floor, (Q - (k - 1)) / C,
0)\`, followed by the documented median-relative cap. Because subject
parcel spaces need not be aligned, heterogeneity is estimated from each
effect row's spatial mean.

The package-wide default remains \`"none"\`. \`"random_effects"\` is an
opt-in policy for materially unbalanced trialwise designs.

## Examples

``` r
dkge_effect_weights("none")
#> $method
#> [1] "none"
#> 
#> $within
#> [1] "auto"
#> 
#> $tau_method
#> [1] "DL"
#> 
#> $max_ratio
#> [1] 10
#> 
#> $variance_floor
#> [1] 1e-08
#> 
#> $tau_floor
#> [1] 0
#> 
#> attr(,"class")
#> [1] "dkge_effect_weights"
dkge_effect_weights("count")
#> $method
#> [1] "count"
#> 
#> $within
#> [1] "auto"
#> 
#> $tau_method
#> [1] "DL"
#> 
#> $max_ratio
#> [1] 10
#> 
#> $variance_floor
#> [1] 1e-08
#> 
#> $tau_floor
#> [1] 0
#> 
#> attr(,"class")
#> [1] "dkge_effect_weights"
dkge_effect_weights("random_effects", within = "auto")
#> $method
#> [1] "random_effects"
#> 
#> $within
#> [1] "auto"
#> 
#> $tau_method
#> [1] "DL"
#> 
#> $max_ratio
#> [1] 10
#> 
#> $variance_floor
#> [1] 1e-08
#> 
#> $tau_floor
#> [1] 0
#> 
#> attr(,"class")
#> [1] "dkge_effect_weights"
```
