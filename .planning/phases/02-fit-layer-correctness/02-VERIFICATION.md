---
phase: 02-fit-layer-correctness
verified: 2026-01-20T03:30:00Z
status: passed
score: 4/4 must-haves verified
must_haves:
  truths:
    - "dkge_fit() output U satisfies K-orthonormality: t(U) %*% K %*% U = I_r within 1e-8 tolerance"
    - "Pooled design computation aggregates subject data correctly: t(R) %*% R = sum_s X_s^T X_s"
    - "Recovery test with known ground truth achieves expected reconstruction error"
    - "Results are deterministic across runs with same seed"
  artifacts:
    - path: "tests/testthat/test-fit.R"
      provides: "K-orthonormality property tests + pooled design tests + determinism tests"
    - path: "tests/testthat/test-toy-recovery.R"
      provides: "Recovery tests at multiple SNR levels and multi-rank"
    - path: "R/dkge-fit-core.R"
      provides: "Core algorithm with K-orthonormal basis computation"
    - path: "R/dkge-sim.R"
      provides: "dkge_sim_toy() and dkge_cosines_K() for ground truth testing"
  key_links:
    - from: "R/dkge-fit-core.R:222"
      to: "tests/testthat/test-fit.R"
      via: "U <- Kihalf %*% U_hat produces K-orthonormal basis"
    - from: "R/dkge-sim.R"
      to: "tests/testthat/test-toy-recovery.R"
      via: "dkge_sim_toy() generates ground truth, dkge_cosines_K() measures subspace agreement"
---

# Phase 2: Fit Layer Correctness Verification Report

**Phase Goal:** Core algorithm produces mathematically valid K-orthonormal bases
**Verified:** 2026-01-20T03:30:00Z
**Status:** passed
**Re-verification:** No - initial verification

## Goal Achievement

### Observable Truths

| # | Truth | Status | Evidence |
|---|-------|--------|----------|
| 1 | dkge_fit() output U satisfies K-orthonormality: t(U) %*% K %*% U = I_r within tolerance | VERIFIED | 9 K-orthonormality tests in test-fit.R all pass with 1e-8 tolerance |
| 2 | Pooled design computation aggregates subject data correctly | VERIFIED | 5 pooled Gram tests verify t(R) %*% R = sum_s X_s^T X_s |
| 3 | Recovery test with known ground truth achieves expected reconstruction error | VERIFIED | 8 recovery tests pass at SNR levels 20/8/2 with cosine thresholds 0.98/0.90/0.70 |
| 4 | Results are deterministic across runs with same seed | VERIFIED | 6 determinism tests verify identical U, Chat, R, weights, eigenvalues |

**Score:** 4/4 truths verified

### Required Artifacts

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| `tests/testthat/test-fit.R` | K-orthonormality + pooled design + determinism tests | VERIFIED | 26 test_that blocks, including 9 K-orthonormality, 5 pooled Gram, 6 determinism |
| `tests/testthat/test-toy-recovery.R` | Recovery tests at multiple SNR levels | VERIFIED | 9 test_that blocks covering high/medium/low SNR + multi-rank + multi-factor |
| `R/dkge-fit-core.R` | Core algorithm with eigendecomposition and K-orthonormalization | VERIFIED | 498 lines, includes .dkge_fit_solve() with eigen() and Kihalf transformation |
| `R/dkge-sim.R` | Ground truth generator and cosine metric | VERIFIED | 150 lines, exports dkge_sim_toy() and dkge_cosines_K() |

### Key Link Verification

| From | To | Via | Status | Details |
|------|-----|-----|--------|---------|
| R/dkge-fit-core.R:222 | tests/testthat/test-fit.R | U <- Kihalf %*% U_hat | WIRED | Implementation computes K-orthonormal basis; tests verify UtKU = I |
| R/dkge-fit-core.R:207 | tests/testthat/test-fit.R | eigen(Chat, symmetric=TRUE) | WIRED | Symmetric eigendecomposition used; tests verify eigenvalues positive |
| R/dkge-sim.R | tests/testthat/test-toy-recovery.R | dkge_sim_toy(), dkge_cosines_K() | WIRED | 8 recovery tests import and call these functions |
| R/dkge-fit-core.R | tests/testthat/test-fit.R | .dkge_compute_shared_ruler() | WIRED | 5 pooled design tests verify Cholesky decomposition R^T R = G_pool |

### Requirements Coverage

Based on ROADMAP success criteria:

| Requirement | Status | Evidence |
|-------------|--------|----------|
| K-orthonormality: t(U) %*% K %*% U = I within tolerance | SATISFIED | 9 tests across identity, RBF, multi-factor kernels + rank edge cases |
| Pooled design aggregates correctly (recovery test) | SATISFIED | 5 explicit Gram tests + 8 recovery tests implicitly verify correct aggregation |
| Recovery with ground truth achieves expected error | SATISFIED | High SNR (20) achieves >0.98, Medium SNR (8) achieves >0.90, Low SNR (2) achieves >0.70 |
| Deterministic results with same seed | SATISFIED | 6 tests verify tolerance=0 exact equality of all outputs |

### Anti-Patterns Found

| File | Line | Pattern | Severity | Impact |
|------|------|---------|----------|--------|
| R/dkge-fit-core.R:435 | 435 | prep() deprecation warning | Info | Pre-existing, does not affect correctness |

No blocking anti-patterns found. The deprecation warning is external (multivarious package) and does not affect DKGE functionality.

### Human Verification Required

None - all success criteria are programmatically verifiable through the test suite.

### Test Coverage Summary

**K-Orthonormality Tests (9):**
1. Identity kernel
2. RBF (ordinal) kernel
3. Multi-factor kernel
4. Rank=1 edge case
5. Rank=q (full rank)
6. Rank=q-1 (near-full)
7. With ridge > 0
8. With MFA subject weights
9. With energy subject weights

**Pooled Design Tests (5):**
1. Cholesky matches Gram matrix (basic)
2. Row standardization (basic)
3. Heterogeneous T_s values
4. Many subjects (S=10)
5. Row standardization scaling

**Determinism Tests (6):**
1. w_method='none'
2. w_method='mfa_sigma1'
3. w_method='energy'
4. solver='pooled'
5. With ridge regularization
6. With Omega_list spatial weights

**Recovery Tests (9):**
1. Single main-effect (original test)
2. High SNR (20) - cosine > 0.98
3. Medium SNR (8) - cosine > 0.90
4. Low SNR (2) - cosine > 0.70
5. Graceful degradation (monotonic)
6. Multi-rank r=2
7. Multi-rank r=3
8. Multi-factor (A+B)
9. Multi-factor with interaction (A:B)

### Full Test Suite Status

```
Duration: 13.2 s
[ FAIL 0 | WARN 148 | SKIP 2 | PASS 938 ]
```

- **Failures:** 0
- **Warnings:** 148 (pre-existing multivarious::prep() deprecation, external to DKGE)
- **Skips:** 2 (long tests disabled, T4transport not available - expected)
- **Passes:** 938

## Verification Conclusion

Phase 2 goal **achieved**. The core algorithm (`dkge_fit()`) produces mathematically valid K-orthonormal bases:

1. **K-orthonormality** is verified across 9 configurations with 1e-8 tolerance
2. **Pooled design computation** correctly aggregates subject data via verified Cholesky decomposition
3. **Recovery tests** demonstrate the algorithm recovers planted signal structure at expected quality levels
4. **Determinism** is verified with tolerance=0 exact equality across 6 configurations

All 4 success criteria from ROADMAP.md are satisfied with comprehensive test coverage.

---
*Verified: 2026-01-20T03:30:00Z*
*Verifier: Claude (gsd-verifier)*
