# Fit an iterative group functional template on fixed support

Fits all subjects to one identified display support while learning the
target's functional features from a typed alignment-feature channel. The
default initialization is the existing one-pass spatial kNN pooling.
Each iteration alternates fixed-policy Sinkhorn fits with the mass-aware
barycentric update that minimizes feature distortion conditional on the
plans. Every updated target signature is then returned to unit L2 norm,
preventing entropic averaging from silently turning off the functional
cost while preserving the shared kernel-image gauge. Every subject is
refitted once more against the converged template.

## Usage

``` r
dkge_fit_functional_template(
  support,
  alignment_features,
  centroids,
  reference_selection = NULL,
  sizes = NULL,
  support_masses = NULL,
  mapper = dkge_mapper_spec("sinkhorn", epsilon = 0.05, lambda_emb = 1, lambda_spa = 0.5,
    sigma_mm = 15, warm_start = FALSE),
  initialization = c("spatial_knn", "supplied_independent"),
  initial_template = NULL,
  subject_weights = NULL,
  knn_k = 8L,
  knn_sigma = 5,
  max_iter = 25L,
  tolerance = 1e-04,
  objective_tolerance = 1e-04,
  update_rate = 0.5,
  rank_tolerance = 1e-09
)
```

## Arguments

- support:

  A
  [`dkge_reference_support()`](https://bbuchsbaum.github.io/dkge/reference/dkge_reference_support.md)
  identifying the display target.

- alignment_features:

  Typed features from
  [`dkge_alignment_features()`](https://bbuchsbaum.github.io/dkge/reference/dkge_alignment_features.md).

- centroids:

  Named list of subject parcel coordinates.

- reference_selection:

  Optional typed selection receipt when `support` was derived from a
  cohort subject. Selection-free supports must declare fixed ancillary
  provenance.

- sizes:

  Optional positive subject parcel masses.

- support_masses:

  Optional positive target-support masses.

- mapper:

  Fixed Sinkhorn transport specification. Functional cost weight must be
  positive and epsilon calibration is not permitted inside template
  fitting.

- initialization:

  Either one-pass spatial kNN pooling or a supplied template on the same
  support with separately identified independent-data provenance.

- initial_template:

  Required for `"supplied_independent"`.

- subject_weights:

  Optional fixed positive subject weights. Named weights are reordered
  by exact subject ID; unnamed weights are positional in fitted subject
  order. Equal subject weighting is the default; fit-level MFA weights
  are never imported.

- knn_k, knn_sigma:

  Spatial initialization controls.

- max_iter:

  Maximum number of alternating updates.

- tolerance:

  Relative target-feature change required for convergence.

- objective_tolerance:

  Relative objective change required jointly with `tolerance` after at
  least two iterations.

- update_rate:

  Under-relaxation in `(0, 1]` toward each exact plan-conditional
  barycentric update. The default damps small deterministic
  Sinkhorn/template oscillations without changing the fixed point.

- rank_tolerance:

  Relative eigenvalue threshold for rank diagnostics. The
  alignment-feature object's declared tolerance is a non-relaxable
  floor; this argument may make the template gate stricter but not
  weaker.

## Value

An immutable `dkge_functional_template` containing final operators,
plans, objective/scale/rank trajectories, numerical diagnostics,
convergence status, and inferential eligibility.

## Details

A fixed anatomical or MNI support supplies coordinates only. This
function makes it a functional target by learning and recording template
features on that support; passing a bare grid directly is not functional
alignment.
