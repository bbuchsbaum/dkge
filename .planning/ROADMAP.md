# Roadmap: dkge Publication Readiness

## Overview

This roadmap transforms the dkge R package from "works for developer" to "publication ready" through systematic test hardening. We build verification from the computational foundation upward: data/kernel layer, fit layer, cross-fitting, edge cases, transport/inference, and finally integration tests. Each phase delivers observable confidence that the package faithfully implements the mathematical specification in `algo.md`.

## Phases

**Phase Numbering:**
- Integer phases (1, 2, 3): Planned milestone work
- Decimal phases (2.1, 2.2): Urgent insertions (marked with INSERTED)

- [x] **Phase 1: Data + Kernel Foundation** - Verify data constructors and kernel operations ✓
- [ ] **Phase 2: Fit Layer Correctness** - Verify K-orthonormality, eigendecomposition, pooled design
- [ ] **Phase 3: Cross-Fitting Validation** - Verify LOSO/K-fold mechanics and analytic fallback
- [ ] **Phase 4: Numerical Edge Cases** - Verify robustness to degenerate inputs
- [ ] **Phase 5: Transport + Inference** - Verify Sinkhorn convergence and null calibration
- [ ] **Phase 6: Integration + S3 Contracts** - Verify end-to-end pipelines and API contracts

## Phase Details

### Phase 1: Data + Kernel Foundation
**Goal**: All downstream tests can trust that data handling and kernel construction are correct
**Depends on**: Nothing (first phase)
**Requirements**: MATH-01, API-01, COV-01, FRAG-01
**Success Criteria** (what must be TRUE):
  1. `dkge_subject()` and `dkge_data()` reject invalid inputs with clear error messages
  2. Effect alignment produces identical results regardless of effect ordering across subjects
  3. `design_kernel()` output is symmetric and positive semidefinite for all valid inputs
  4. Kernel root reconstruction satisfies `Khalf %*% Khalf = K` within numerical tolerance
**Plans**: 2 plans

Plans:
- [x] 01-01-PLAN.md — Data constructor input validation and ordering invariance tests ✓
- [x] 01-02-PLAN.md — Kernel mathematical invariant and edge case tests ✓

### Phase 2: Fit Layer Correctness
**Goal**: Core algorithm produces mathematically valid K-orthonormal bases
**Depends on**: Phase 1
**Requirements**: MATH-01, PROP-01, REG-01
**Success Criteria** (what must be TRUE):
  1. `dkge_fit()` output U satisfies K-orthonormality: `t(U) %*% K %*% U = I` within tolerance
  2. Pooled design computation aggregates subject data correctly (verified by recovery test)
  3. Recovery test with known ground truth achieves expected reconstruction error
  4. Results are deterministic across runs with same seed
**Plans**: 2 plans

Plans:
- [ ] 02-01-PLAN.md — K-orthonormality property tests and pooled design verification
- [ ] 02-02-PLAN.md — Recovery tests with synthetic ground truth and determinism tests

### Phase 3: Cross-Fitting Validation
**Goal**: LOSO and K-fold cross-fitting are unbiased with no data leakage
**Depends on**: Phase 2
**Requirements**: MATH-01, FRAG-01, COV-01
**Success Criteria** (what must be TRUE):
  1. LOSO basis U^{-s} differs from full basis U (held-out subject excluded from computation)
  2. Injecting extreme values in held-out subject does not affect U^{-s}
  3. K-fold with k=S subjects produces identical results to LOSO
  4. All analytic LOSO fallback paths (eigengap, perturbation_magnitude, solver_not_pooled) have test coverage
  5. Analytic approximation matches iterative LOSO within specified tolerance
**Plans**: TBD

Plans:
- [ ] 03-01: LOSO mechanics and leakage tests
- [ ] 03-02: Analytic fallback path coverage

### Phase 4: Numerical Edge Cases
**Goal**: Package handles degenerate inputs gracefully without silent failures
**Depends on**: Phase 2
**Requirements**: EDGE-01, FRAG-01
**Success Criteria** (what must be TRUE):
  1. Rank-deficient input matrices produce informative errors or graceful degradation
  2. Near-singular matrices do not cause NaN propagation
  3. Tests pass with multiple different random seeds (not seed-dependent)
  4. Effect alignment handles partial overlap between subjects correctly
**Plans**: TBD

Plans:
- [ ] 04-01: Rank-deficient and ill-conditioned inputs
- [ ] 04-02: Multi-seed robustness tests

### Phase 5: Transport + Inference
**Goal**: Sinkhorn transport converges correctly and permutation tests are calibrated
**Depends on**: Phase 2
**Requirements**: MATH-01, PROP-01
**Success Criteria** (what must be TRUE):
  1. Sinkhorn output is doubly stochastic (rows and columns sum to 1) within tolerance
  2. Deterministic transport cases (identity mapping) produce expected exact results
  3. Null distribution p-values are uniform under permutation (chi-square test passes)
  4. `parallel=TRUE` produces identical results to sequential execution
**Plans**: TBD

Plans:
- [ ] 05-01: Sinkhorn convergence tests
- [ ] 05-02: Null calibration tests

### Phase 6: Integration + S3 Contracts
**Goal**: End-to-end workflows and user-facing API behave as documented
**Depends on**: Phases 1-5
**Requirements**: API-01, CHECK-01, COV-01
**Success Criteria** (what must be TRUE):
  1. `dkge_pipeline()` completes successfully on valid multi-subject data
  2. All exported S3 methods (`print`, `predict`, `as.data.frame`) dispatch correctly
  3. `R CMD check` passes with 0 errors, 0 warnings, 0 notes
  4. All exported functions have working `@examples`
  5. Test coverage on exported functions exceeds 80%
**Plans**: TBD

Plans:
- [ ] 06-01: Pipeline integration tests
- [ ] 06-02: S3 method and API contract tests
- [ ] 06-03: R CMD check compliance

## Requirement Mapping

| Requirement | Description | Phase |
|-------------|-------------|-------|
| MATH-01 | Mathematical accuracy verified | 1, 2, 3, 5 |
| API-01 | API contracts verified | 1, 6 |
| CHECK-01 | R CMD check passes | 6 |
| COV-01 | Higher test coverage | 1, 3, 6 |
| PROP-01 | Property-based tests | 2, 5 |
| REG-01 | Regression tests | 2 |
| EDGE-01 | Edge case coverage | 4 |
| FRAG-01 | Fragile areas hardened | 1, 3, 4 |

**Coverage:** 8/8 requirements mapped

## Progress

**Execution Order:**
Phases execute in numeric order: 1 -> 2 -> 3 -> 4 -> 5 -> 6
(Phases 3-5 depend on Phase 2; Phase 6 depends on all others)

| Phase | Plans Complete | Status | Completed |
|-------|----------------|--------|-----------|
| 1. Data + Kernel Foundation | 2/2 | Complete ✓ | 2026-01-19 |
| 2. Fit Layer Correctness | 0/2 | Planned | - |
| 3. Cross-Fitting Validation | 0/2 | Not started | - |
| 4. Numerical Edge Cases | 0/2 | Not started | - |
| 5. Transport + Inference | 0/2 | Not started | - |
| 6. Integration + S3 Contracts | 0/3 | Not started | - |

---
*Created: 2026-01-19*
*Phase 1 planned: 2026-01-19*
*Phase 1 complete: 2026-01-19*
*Phase 2 planned: 2026-01-19*
