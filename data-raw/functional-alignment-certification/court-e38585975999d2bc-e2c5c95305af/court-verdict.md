# Functional-alignment v9 court verdict

Status: `complete_approximate_only_inferential_promotion_blocked`.

Exact-oracle gate: **TRUE**. Frozen-plan promotion gate: **FALSE**.

The run contains 5280 arm-level results from 840 fixed cohorts. No seed, threshold, epsilon, or adverse cell was changed. V9 inherits v3's sole predeclared Sinkhorn iteration-ceiling amendment after v2 stopped without statistical outputs; v4-v9 change only candidate/documentation binding, package-gate and supersession provenance, and directory labels.

The table below reports the predeclared S, estimation-rank, kernel-rank, contrast-span, and eigengap comparisons with Wilson 95% intervals. It is an OAT/fractional grid, so differences are descriptive rather than causal factor effects.

| Factor | Setting | Arm | FWER | Wilson 95% CI | Convergence |
|---|---:|---|---:|---:|---:|
| S | 16 | geometry_only | 0.025 | [0.004, 0.129] | 1.000 |
| S | 16 | independent_alignment | 0.075 | [0.026, 0.199] | 1.000 |
| S | 16 | kernel_image_residual_prototype | 0.050 | [0.014, 0.165] | 1.000 |
| S | 16 | full_permutation_reestimation | 0.050 | [0.014, 0.165] | 1.000 |
| S | 32 | geometry_only | 0.025 | [0.004, 0.129] | 1.000 |
| S | 32 | independent_alignment | 0.150 | [0.071, 0.291] | 1.000 |
| S | 32 | kernel_image_residual_prototype | 0.075 | [0.026, 0.199] | 1.000 |
| S | 8 | geometry_only | 0.025 | [0.004, 0.129] | 1.000 |
| S | 8 | independent_alignment | 0.100 | [0.040, 0.231] | 1.000 |
| S | 8 | kernel_image_residual_prototype | 0.025 | [0.004, 0.129] | 1.000 |
| estimation_rank | 2 | geometry_only | 0.000 | [0.000, 0.088] | 1.000 |
| estimation_rank | 2 | independent_alignment | 0.150 | [0.071, 0.291] | 1.000 |
| estimation_rank | 2 | kernel_image_residual_prototype | 0.050 | [0.014, 0.165] | 1.000 |
| estimation_rank | 2 | full_permutation_reestimation | 0.000 | [0.000, 0.088] | 1.000 |
| estimation_rank | 4 | geometry_only | 0.000 | [0.000, 0.088] | 1.000 |
| estimation_rank | 4 | independent_alignment | 0.000 | [0.000, 0.088] | 1.000 |
| estimation_rank | 4 | kernel_image_residual_prototype | 0.000 | [0.000, 0.088] | 1.000 |
| kernel_rank | 4 | geometry_only | 0.025 | [0.004, 0.129] | 1.000 |
| kernel_rank | 4 | independent_alignment | 0.075 | [0.026, 0.199] | 1.000 |
| kernel_rank | 4 | kernel_image_residual_prototype | 0.000 | [0.000, 0.088] | 1.000 |
| kernel_rank | 6 | geometry_only | 0.050 | [0.014, 0.165] | 1.000 |
| kernel_rank | 6 | independent_alignment | 0.075 | [0.026, 0.199] | 1.000 |
| kernel_rank | 6 | kernel_image_residual_prototype | 0.075 | [0.026, 0.199] | 1.000 |
| kernel_rank | 6 | full_permutation_reestimation | 0.050 | [0.014, 0.165] | 1.000 |
| contrast_family_dimension | 1 | geometry_only | 0.050 | [0.014, 0.165] | 1.000 |
| contrast_family_dimension | 1 | independent_alignment | 0.075 | [0.026, 0.199] | 1.000 |
| contrast_family_dimension | 1 | kernel_image_residual_prototype | 0.075 | [0.026, 0.199] | 1.000 |
| contrast_family_dimension | 1 | full_permutation_reestimation | 0.050 | [0.014, 0.165] | 1.000 |
| contrast_family_dimension | 2 | geometry_only | 0.025 | [0.004, 0.129] | 1.000 |
| contrast_family_dimension | 2 | independent_alignment | 0.225 | [0.123, 0.375] | 1.000 |
| contrast_family_dimension | 2 | kernel_image_residual_prototype | 0.075 | [0.026, 0.199] | 1.000 |
| eigengap | large | geometry_only | 0.050 | [0.014, 0.165] | 1.000 |
| eigengap | large | independent_alignment | 0.075 | [0.026, 0.199] | 1.000 |
| eigengap | large | kernel_image_residual_prototype | 0.075 | [0.026, 0.199] | 1.000 |
| eigengap | large | full_permutation_reestimation | 0.050 | [0.014, 0.165] | 1.000 |
| eigengap | small | geometry_only | 0.025 | [0.004, 0.129] | 1.000 |
| eigengap | small | independent_alignment | 0.050 | [0.014, 0.165] | 1.000 |
| eigengap | small | kernel_image_residual_prototype | 0.050 | [0.014, 0.165] | 1.000 |
| eigengap | small | full_permutation_reestimation | 0.050 | [0.014, 0.165] | 1.000 |

Same-data residualized alignment remains approximate regardless of favorable cells. Independent features remove direct reuse in correspondence but do not remove the rank-truncated estimator's latent-span dependence. Geometry-only is an ancillary bound, not evidence of functional alignment.

Held-out known-warp power is reported separately and cannot override this verdict.
