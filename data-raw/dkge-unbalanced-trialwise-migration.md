# Migrating to reliability-aware trialwise DKGE

This note describes the Tier 0-2 changes that supersede the earlier
partial-effect implementation. The historical defaults remain available, but
the fit now records enough provenance to tell a legacy observed second moment
from a reliability-weighted or noise-corrected estimator.

## What changed

- `dkge_subject()` now honours `observed_rows` even when all subjects initially
  use the same row labels. A zero-filled row is not interpreted as observed
  merely because its label exists.
- `dkge_data()` creates a deterministic, named union of local effect spaces.
  Missing rows are structural zeros with an explicit observation mask.
- `dkge_effect_grid()` defines a global cell order. Its `scope` field is
  metadata for contrast routing and does not change the kernel.
- `dkge_effect_weights()` adds count, explicit-precision, and capped
  DerSimonian-Laird random-effects pooling. No reliability policy is activated
  silently: the compatibility default is still `method = "none"`.
- `dkge_trial_subject()` accepts trialwise beta maps and an IID, covariance, or
  precision error model. `dkge_trial_subject_chunks()` performs the same
  reduction feature block by feature block without materializing the full
  trial-by-feature response.
- `dkge_fit(debias = "analytic")` subtracts the supplied or estimated
  effect-noise moment. `debias = "split_half"` uses a symmetric cross-half
  moment. The stored `moment_estimator` states whether analytic correction was
  IID or covariance-aware and whether a split was independently justified or
  exploratory.
- Reliability-normalized, debiased, ridge, CPCA, and joint-diagonalization fits
  use `representation = "qspace_moment"`. APIs requiring a physical stacked
  block factor fail with an actionable representation error; q-space
  projection, contrasts, fold refits, and supported cell classification remain
  available.
- LOSO, k-fold, and literal Poisson resamples refit the ruler, reliability,
  random-effects variance, noise correction, and missingness policy using only
  the training-subject multiset.

## Storage and scale

The dense trialwise constructor does not retain `Y`, but it does retain the
shared `T x q` design and a `q x P` beta matrix. The chunked constructor also
omits the reconstructible `q x P` effect score, since it equals
`effect_information %*% beta`. Neither constructor nor the group fit forms a
`P x P` matrix.

## Compatibility choices

To reproduce historical fits, keep `effect_weights =
dkge_effect_weights("none")`, `debias = "none"`, and `missingness = "none"`.
These choices preserve equal observed-row weighting and the uncorrected second
moment. They are not recommended merely because trial counts happen to be
available: choose an error model based on the data-generating and first-level
estimation process.

For materially unbalanced cell means, analytic covariance-aware correction is
the primary protection against the diagonal `1 / n_sc` noise term. Count or
random-effects weighting changes influence but does not, by itself, remove
that bias. Random-effects precision is capped to prevent a few high-count
subjects from owning the geometry; its between-subject variance estimate uses
effect-row spatial means because subject parcel coordinates need not align.

## Required review after migration

Inspect `summary(fit)` before downstream inference. It reports effect coverage,
pair effective sample size, negative spectral mass, trial error provenance,
the moment-estimator name, and the representation contract. Negative spectral
mass can legitimately appear after analytic subtraction or split cross-products;
the eigensolver retains only the positive numerical spectrum and reports rank
reduction rather than silently turning negative directions into components.
