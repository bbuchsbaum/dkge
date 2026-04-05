# Create a DKGE voxel-weight specification

Constructs a lightweight object describing how voxel-level weights are
resolved and applied when building weighted second moments inside DKGE.
The specification is attached to fits and reused by fold builders to
ensure cross-fitting remains leak-free.

## Usage

``` r
dkge_weights(
  prior = NULL,
  adapt = c("none", "kenergy", "precision", "kenergy_prec"),
  combine = c("product", "sum", "override_adapt", "prefer_prior"),
  mix = 0.6,
  shrink = list(alpha = 0.5, winsor = 0.99, normalize = "mean", roi_smooth = FALSE),
  scope = c("fold", "subject"),
  k_weight = NULL,
  collapse = NULL,
  roi = NULL
)
```

## Arguments

- prior:

  Optional prior weights: numeric vector of length V, logical mask,
  integer indices, or ROI labels (factor/character/integer of length V)
  per voxel. Helpers \[dkge_weights_prior_mask()\] and
  \[dkge_weights_prior_roi()\] ease construction.

- adapt:

  Adaptive weighting rule, one of \`"none"\`, \`"kenergy"\`,
  \`"precision"\`, or \`"kenergy_prec"\`.

- combine:

  How prior and adaptive sources combine: \`"product"\` (default),
  \`"sum"\`, \`"override_adapt"\`, or \`"prefer_prior"\`.

- mix:

  Numeric in \[0,1\] controlling the relative influence of the adaptive
  component. Interpreted in log-space for \`combine = "product"\`.

- shrink:

  List with fields \`alpha\` (shrink towards uniform), \`winsor\` (upper
  quantile cap), \`normalize\` (\`"mean"\` or \`"sum"\`), and optional
  \`roi_smooth = TRUE\` to median-smooth within ROIs.

- scope:

  Either \`"fold"\` (default: compute adapt weights from training
  subjects within each fold) or \`"subject"\` (per-subject adaptive
  weights averaged for fold pooling).

- k_weight:

  Optional effect-space kernel for k-energy rules. When \`NULL\` we
  reuse the kernel stored in the fit, with optional factor \`collapse\`.

- collapse:

  Optional list describing factor collapses (e.g., \`list(time =
  "mean")\` or \`list(time = list(method = "mean", window = 3:8))\`).

- roi:

  Optional ROI labels used when \`shrink\$roi_smooth = TRUE\`.

## Value

Object of class \`"dkge_weights"\`.

## Combine modes

- \`"product"\`:

  Geometric mean of prior and adaptive weights, blended in log-space by
  \`mix\` (default). Smooth and robust.

- \`"sum"\`:

  Linear mixture: \`(1 - mix) \* prior + mix \* adaptive\`.

- \`"override_adapt"\`:

  Adaptive weight is used wherever it is finite and positive; falls back
  to prior otherwise. Adaptive takes priority.

- \`"prefer_prior"\`:

  Prior weight is used wherever it is finite and positive; falls back to
  adaptive otherwise. Prior takes priority.

## Examples

``` r
# Default specification with adaptive k-energy weighting
w <- dkge_weights(adapt = "kenergy")
print(w)
#> dkge weight specification
#>   prior   : uniform 
#>   adapt   : kenergy (scope = fold )
#>   combine : product (mix = 0.60 )
#>   shrink  : alpha =0.50, winsor =0.990, normalize =mean

# Uniform (no adaptive) weights
w_uniform <- dkge_weights(adapt = "none")
```
