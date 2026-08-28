# Functional-alignment efficacy: gate calibration

Date: 2026-08-26. Module: `inst/validation/functional-alignment-benefit/benefit.R`.
Gates asserted in `tests/testthat/test-functional-alignment-benefit.R` on
disjoint seeds 301-304.

## Why a second DGP was needed

The frozen Type-I court cannot demonstrate benefit. In its generator every
parcel sits at `seq(0, 1, length.out = P)` and centroids are that same
coordinate plus 0.01 jitter (`court.R:221-226`), so ground-truth correspondence
is recoverable from geometry alone. Its functional patterns additionally carry
a subject-dependent phase (`0.17 * s`, `court.R:203-205`), which points *away*
from the truth. That is the right design for a validity court and useless for
an efficacy claim.

## A monotone warp is not enough

The first version of this generator used a monotone Beta-CDF reparametrisation.
It failed to create any difficulty at all:

| warp | monotone frac | functional latent err | geometry latent err |
|------|---------------|-----------------------|---------------------|
| 0.00 | 1.00          | 0.0124                | 0.0029              |
| 0.90 | 1.00          | 0.0166                | 0.0054              |

On a 1-D support with balanced uniform masses, entropic OT recovers the
monotone rearrangement, which matches rank to rank -- exactly the true latent
correspondence -- however severe the warp. Geometry therefore scores
near-perfectly and no functional signal can add anything. The generator now
applies a deliberately **non-monotone** displacement; `monotone_fraction` is
reported so this cannot silently regress.

## Calibration run (seeds 201-215, S = 12, 24 base parcels)

Regime B, `warp_amount = 0` (monotone fraction 1.00):

| metric                | functional | geometry |
|-----------------------|------------|----------|
| latent error          | 0.0143     | 0.0031   |
| subject reconstruction| 0.956      | 0.969    |
| diffusion             | 2.08       | 2.08     |

Correlation gain, functional minus geometry: mean -0.014, range
[-0.023, -0.007]. Functional wins 0/15. Geometry is correctly preferred.

Regime A, `warp_amount = 0.25` (monotone fraction 0.70):

| metric                | functional | geometry |
|-----------------------|------------|----------|
| latent error          | 0.1539     | 0.2096   |
| subject reconstruction| 0.622      | 0.307    |
| diffusion             | 2.14       | 2.96     |

Latent-error ratio geometry/functional: mean 1.39, min 1.18. Correlation gain
mean 0.315, range [0.118, 0.485]. Functional wins 15/15. Chance baseline error
0.263.

The functional arm's plans are **sharper** than the geometry arm's (2.14 vs
2.96 effective sources per target), so the advantage is not a smoothing
artefact. That comparison is asserted separately.

## Gates

Set below the observed minima with margin:

| gate                        | value | observed        |
|-----------------------------|-------|-----------------|
| `min_error_ratio`           | 1.10  | min 1.18        |
| `min_cor_gain`              | 0.08  | min 0.118       |
| `max_cor_gain_corresponded` | 0.02  | max -0.007      |
| `max_cor_loss_corresponded` | 0.05  | max loss 0.023  |
| `max_diffusion_excess`      | 0.50  | -0.82 (sharper) |
| `min_chance_ratio`          | 1.50  | 1.71            |
