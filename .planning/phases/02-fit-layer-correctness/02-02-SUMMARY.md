---
phase: 02-fit-layer-correctness
plan: 02
subsystem: fit-verification
tags: [dkge_fit, recovery, determinism, testing, SNR, multi-rank]

dependency-graph:
  requires: ["02-01"]
  provides:
    - "Recovery tests with known ground truth at multiple SNR levels"
    - "Multi-rank and multi-factor recovery verification"
    - "Determinism tests for dkge_fit() reproducibility"
  affects: ["03-*"]

tech-stack:
  added: []
  patterns:
    - "dkge_cosines_K() for sign-invariant subspace comparison"
    - "dkge_sim_toy() for ground truth generation"
    - "withr::local_seed() for deterministic test fixtures"

key-files:
  created: []
  modified:
    - tests/testthat/test-toy-recovery.R
    - tests/testthat/test-fit.R

decisions:
  - id: "02-02-snr-thresholds"
    choice: "SNR thresholds: 20->0.98, 8->0.90, 2->0.70"
    reason: "Empirically validated recovery quality at each SNR level"
  - id: "02-02-graceful-degradation"
    choice: "Verify monotonic decrease with 0.05 tolerance for sampling variance"
    reason: "Allows minor variation while ensuring overall trend"

metrics:
  duration: "2 min"
  completed: "2026-01-20"
---

# Phase 02 Plan 02: Recovery and Determinism Tests Summary

**One-liner:** Recovery tests at multiple SNR levels (20/8/2) + multi-rank/multi-factor + determinism verification for all w_methods

## What Was Done

### Task 1: Expand Recovery Tests with Multiple SNR Levels
Added 8 new recovery tests to `test-toy-recovery.R`:

1. **High SNR (20):** Achieves cosine > 0.98
2. **Medium SNR (8):** Achieves cosine > 0.90
3. **Low SNR (2):** Achieves cosine > 0.70 (graceful degradation)
4. **Graceful degradation test:** Verifies monotonic decrease in recovery quality
5. **Multi-rank r=2:** Both components recovered with cosine > 0.90
6. **Multi-rank r=3:** All three components recovered with cosine > 0.85
7. **Multi-factor (A+B):** Two main effects recovered with cosine > 0.85
8. **Multi-factor with interaction (A, B, A:B):** All components including interaction recovered

Key implementation details:
- Used `dkge_cosines_K()` for sign-invariant subspace comparison in K-metric
- Increased S (subjects) and P (clusters) for stability at low SNR
- All tests use `withr::local_seed()` for reproducibility

### Task 2: Add Determinism Tests
Added 6 determinism tests to `test-fit.R`:

1. **w_method='none':** Exact match of U, Chat, R, weights, eigenvalues
2. **w_method='mfa_sigma1':** MFA weighting deterministic with fixed seed
3. **w_method='energy':** Energy weighting deterministic
4. **solver='pooled':** Pooled solver deterministic
5. **With ridge regularization:** Ridge parameter doesn't affect determinism
6. **With Omega_list:** Spatial weights don't affect determinism

All tests verify tolerance=0 exact equality.

### Task 3: Full Test Suite Verification
- **Total tests:** 938 passing
- **Failures:** 0
- **Skips:** 2 (expected: long tests, T4transport package)
- **Warnings:** 148 (pre-existing multivarious::prep() deprecation)

## Commits

| Hash | Description |
|------|-------------|
| 3c1390f | test(02-02): expand recovery tests with multiple SNR levels and multi-rank |
| 92e8b34 | test(02-02): add determinism tests for dkge_fit() |

## Verification Results

Recovery and determinism verified:
- [x] High SNR (20) achieves cosine > 0.98
- [x] Medium SNR (8) achieves cosine > 0.90
- [x] Low SNR (2) achieves cosine > 0.70
- [x] Multi-rank recovery (r >= 2) works correctly
- [x] Determinism: same seed produces identical U, Chat, R, weights
- [x] Tests use withr::local_seed() for reproducibility
- [x] Tests use dkge_cosines_K() for sign-invariant subspace comparison

## Deviations from Plan

None - plan executed exactly as written.

## Test Coverage Added

| Test File | New Tests | Description |
|-----------|-----------|-------------|
| test-toy-recovery.R | 8 | SNR levels, multi-rank, multi-factor recovery |
| test-fit.R | 6 | Determinism across w_methods and configurations |
| **Total** | **14** | Exceeds plan requirement of 8 |

## Next Phase Readiness

Phase 2 (Fit Layer Correctness) is now complete with:
- K-orthonormality property tests (02-01)
- Pooled design verification (02-01)
- Recovery tests with ground truth (02-02)
- Determinism tests (02-02)

Ready to proceed to Phase 3 (Contrast + Inference Correctness).
