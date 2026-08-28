# Functional-alignment efficacy study

Companion to `../functional-alignment-court`. The court asks **does alignment
leak?**; this module asks **does alignment do anything?**

## Why the court cannot answer this

The court's generator places every parcel at `seq(0, 1, length.out = P)` and
sets centroids to that coordinate plus 0.01 jitter, so ground-truth
correspondence is recoverable from geometry alone. Its functional patterns also
carry a subject-dependent phase (`0.17 * s`), which points away from the truth.
Both choices are correct for a validity court -- they make any Type-I inflation
unambiguously attributable to leakage -- and both make an efficacy claim
impossible to demonstrate there.

## The design

Ground truth is a latent coordinate `u`; anatomy is `x = u + amount * d_s(u)`
for a subject-specific displacement `d_s`. A functional signature indexed by
`u` and **shared** across subjects supplies the signal alignment can exploit.

Two regimes, and both halves matter:

| regime          | `warp_amount` | monotone frac | expectation                    |
|-----------------|---------------|---------------|--------------------------------|
| corresponded    | 0             | 1.00          | functional must NOT help       |
| idiosyncratic   | 0.25          | ~0.70         | functional must clearly help   |

A method that improves in both regimes is measuring smoothing, not alignment.

## The displacement must be non-monotone

This is the trap the first version fell into. On a 1-D support with balanced
uniform masses, entropic OT recovers the monotone rearrangement, which matches
rank to rank -- exactly the true latent correspondence -- no matter how severe
a *monotone* warp is. Geometry then scores near-perfectly and no functional
signal can add anything. `fab_displace()` is therefore deliberately not
monotone, and `fab_generate()` reports `monotone_fraction` so a regression to
the degenerate design is caught by an assertion rather than by a silently
uninformative result.

## Metrics

- `fab_latent_error()` -- mean absolute latent-coordinate error of each
  operator's induced correspondence. Primary, because it is a *position* and so
  is robust to how much the plan smooths.
- `fab_reconstruction()` -- correlation of the transported field with the truth
  on the reference grid, evaluated on a held-out replicate (`B_eval`) that
  never built the alignment features.
- `fab_operator_diffusion()` -- effective sources per target, `exp(H(row))`.
  Guards against buying reconstruction with extra averaging.
- `fab_operator_symmetry()` -- per-subject diffusion and attenuation, and the
  reference subject's z-score among the others.

Alignment features come from an independent replicate (`B_align`), so efficacy
is measured in the leak-free mode and is not entangled with the court's
question.

## Calibration

Gates in `fab_gates()` were calibrated on seeds 201-215 and are asserted on
disjoint seeds 301-304. See
`data-raw/functional-alignment-benefit/calibration-2026-08-26.md`.
