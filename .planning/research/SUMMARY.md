# Project Research Summary

**Project:** dkge - Design-Kernel Group Embedding R Package Testing
**Domain:** Scientific R package testing for publication readiness
**Researched:** 2026-01-19
**Confidence:** HIGH

## Executive Summary

The dkge package is a sophisticated statistical/numerical R package implementing kernel-based fMRI analysis methods. Based on research into R package testing best practices, the path to publication readiness requires a layered testing approach that prioritizes mathematical correctness over feature coverage. The package already has a substantial test foundation (~60 test files), but needs systematic completion of invariant-based tests, cross-fitting validation, and edge case coverage.

The recommended approach is to build tests from the computational foundation upward: data layer, kernel layer, fit layer, contrast/cross-fitting layer, transport layer, and finally integration tests. This order respects the dependency structure where higher layers depend on lower layers being correct. The existing test patterns in the codebase (recovery tests, K-orthonormality checks, equivariance tests) provide excellent templates to follow.

The primary risks are subtle numerical correctness issues that could undermine scientific validity: eigenvector sign ambiguity causing false test failures, uncalibrated tolerances masking bugs, and data leakage in LOSO cross-validation invalidating the core methodological claim. These are mitigated by testing mathematical invariants rather than specific values, deriving tolerances from numerical analysis theory, and explicitly verifying held-out subjects are excluded from basis computation.

## Key Findings

### Recommended Stack

The testing stack should build on the existing testthat 3 infrastructure with targeted additions for numerical code validation.

**Core technologies:**
- **testthat 3.3.2**: Primary testing framework (already configured) - Edition 3 with waldo-based comparisons and explicit tolerance support
- **hedgehog 0.2**: Property-based testing for mathematical invariants - essential for testing commutativity, idempotency, and bound verification without manual case enumeration
- **covr 3.6.5**: Code coverage measurement - target 80%+ overall, 90%+ on critical paths (fit, contrast, LOSO)
- **lintr 3.3.0**: Static analysis to catch potential bugs before they become test failures

**CI infrastructure:**
- GitHub Actions with r-lib/actions v2 (check-standard, test-coverage workflows)
- Multi-platform testing: Linux, macOS, Windows
- R version matrix: release, devel, oldrel-1

### Expected Features

**Must have (table stakes):**
- Mathematical invariant tests: K-orthonormality (t(U) %*% K %*% U = I), kernel reconstruction, symmetry, PSD
- API contract tests: input validation with clear error messages, return value structure verification
- Cross-platform CI: R CMD check passes with 0 errors/warnings/notes on all platforms
- Code coverage > 80% on exported functions
- All exported functions have working @examples

**Should have (competitive):**
- Recovery tests: synthetic data with known ground truth demonstrates method correctness
- Property-based tests: hedgehog tests for theoretical bounds and monotonicity
- Reference implementation comparison: already done for Sinkhorn vs T4transport
- Numerical snapshot tests: detect unintended algorithm changes with tolerance-based comparison

**Defer (v2+):**
- Full statistical power analysis (simulation study)
- Cross-language validation (Python/NumPy comparison)
- Extensive performance benchmarks beyond O(q^3) verification
- Interactive test fixtures

### Architecture Approach

Tests should follow a layered hierarchy mirroring the computational dependency structure. Each layer tests its own invariants and assumes lower layers are correct. The existing test organization already demonstrates this pattern well.

**Major components and their test responsibilities:**

1. **Data Layer (Layer 0)**: dkge_subject, dkge_data, effect alignment - test input validation, name alignment, provenance tracking
2. **Kernel Layer (Layer 1)**: design_kernel, kernel_roots - test symmetry, PSD, reconstruction (Khalf %*% Khalf = K)
3. **Fit Layer (Layer 2)**: dkge_fit, eigendecomposition - test K-orthonormality of U, pooled design correctness
4. **Contrast Layer (Layer 3)**: LOSO, K-fold, analytic approximation - test unbiasedness, equivariance, fallback conditions
5. **Transport Layer (Layer 4)**: Sinkhorn OT - test doubly stochastic plans, marginal constraints
6. **Inference Layer (Layer 5)**: permutation tests - test null calibration, FWER control
7. **Pipeline Layer (Layer 6)**: end-to-end integration tests

### Critical Pitfalls

1. **Eigenvector sign/ordering ambiguity** - Never compare eigenvectors directly; test K-orthonormality (t(U) %*% K %*% U = I), reconstruction error, or subspace angles. The existing pattern at test-contrast.R is correct.

2. **Tolerance calibration blindness** - Derive tolerances from theory (condition_number * machine_eps), not magic numbers. Document tolerance rationale in test comments. Use relative tolerances scaled to problem size.

3. **Data leakage in cross-validation** - Explicitly verify that LOSO basis U^{-s} differs from full U and that held-out subject's data was not used. Test by injecting extreme values in held-out subject.

4. **Analytic approximation fallback not tested** - The dkge-analytic.R code has multiple fallback paths (eigengap, perturbation_magnitude, solver_not_pooled, etc.). Each needs a test that triggers it and verifies correctness.

5. **Effect alignment order sensitivity** - Results must not depend on effect ordering. Test with shuffled effect names and verify identical results.

## Implications for Roadmap

Based on research, suggested phase structure:

### Phase 1: Data + Kernel Foundation
**Rationale:** All other tests depend on correct data handling and kernel construction. Build the foundation first.
**Delivers:** Complete test coverage of data constructors and kernel operations
**Addresses:** Input validation, effect alignment, kernel symmetry/PSD, kernel root reconstruction
**Avoids:** Testing higher layers on potentially broken foundation
**Estimated effort:** 1-2 days

### Phase 2: Fit Layer Correctness
**Rationale:** Core algorithm - if fit is wrong, all downstream results are wrong
**Delivers:** Verified eigendecomposition, K-orthonormal bases, pooled design correctness
**Addresses:** Mathematical invariant tests (K-orthonormality), recovery tests with known ground truth
**Avoids:** Eigenvector ambiguity pitfall by testing invariants not values; tolerance calibration blindness
**Estimated effort:** 2-3 days

### Phase 3: Cross-Fitting Validation
**Rationale:** LOSO and K-fold are the primary scientific contribution - must be rigorously correct
**Delivers:** Verified unbiasedness, no data leakage, analytic approximation correctness
**Addresses:** LOSO mechanics, K-fold equivalence (k=S case), analytic fallback paths
**Avoids:** Data leakage pitfall, analytic fallback not tested pitfall
**Estimated effort:** 3-4 days

### Phase 4: Numerical Edge Cases
**Rationale:** Real fMRI data has degenerate cases; tests must cover them
**Delivers:** Robustness to rank-deficient inputs, near-singular matrices, partial effect overlap
**Addresses:** Multi-seed testing, zero input handling, extreme values, NaN propagation
**Avoids:** Seed-dependent tests, missing rank-deficient tests
**Estimated effort:** 2-3 days

### Phase 5: Transport + Inference
**Rationale:** Can run in parallel; both needed for complete package but less critical than core
**Delivers:** Verified Sinkhorn convergence, doubly stochastic plans, null calibration
**Addresses:** Deterministic transport cases, p-value uniformity under null
**Avoids:** Parallel execution non-determinism (test parallel vs sequential equivalence)
**Estimated effort:** 2-3 days

### Phase 6: Integration + S3 Contracts
**Rationale:** Final layer; test complete workflows and user-facing API
**Delivers:** End-to-end pipeline tests, S3 method dispatch, print/predict/as.data.frame
**Addresses:** All exported S3 methods, argument handling, subclass behavior
**Avoids:** S3 method dispatch testing gap, undocumented test dependencies
**Estimated effort:** 2-3 days

### Phase Ordering Rationale

- Layers 0-2 must come first: all higher tests depend on data/kernel/fit being correct
- Phase 3 (cross-fitting) is the core scientific claim - prioritize after fit is verified
- Phases 4-5 can partially parallelize: edge cases and transport/inference are semi-independent
- Phase 6 comes last: integration tests assume all components work

### Research Flags

Phases likely needing deeper research during planning:
- **Phase 3 (Cross-Fitting):** Complex fallback logic in analytic approximation needs careful test design
- **Phase 5 (Inference):** Null calibration requires simulation study design

Phases with standard patterns (skip research-phase):
- **Phase 1 (Data/Kernel):** Well-documented testthat patterns; existing tests provide templates
- **Phase 2 (Fit):** Matrix invariant testing is standard numerical practice
- **Phase 6 (Integration):** Standard R package testing patterns

## Confidence Assessment

| Area | Confidence | Notes |
|------|------------|-------|
| Stack | HIGH | CRAN package versions verified Jan 2026; testthat 3.3.2 and hedgehog 0.2 current |
| Features | HIGH | Based on CRAN submission checklist, rOpenSci guidelines, Bioconductor standards |
| Architecture | HIGH | Based on analysis of existing dkge test suite (60+ files) plus R Packages (2e) |
| Pitfalls | HIGH | Based on R package testing literature, numerical analysis theory, eigen() documentation |

**Overall confidence:** HIGH

### Gaps to Address

- **Analytic fallback coverage:** Current test-analytic.R tests some but not all fallback paths. Need inventory of all paths and dedicated tests.
- **Multi-seed robustness:** Tests use consistent seeds but no systematic multi-seed verification. Add during Phase 4.
- **Parallel/sequential equivalence:** No explicit tests that parallel=TRUE gives same results as sequential. Add during Phase 5.
- **Effect alignment order:** Test at test-project.R:97-102 exists but dkge_align_effects.R has complex ordering logic. Expand during Phase 4.

## Sources

### Primary (HIGH confidence)
- [testthat CRAN](https://cran.r-project.org/package=testthat) - Version 3.3.2, published 2026-01-11
- [covr CRAN](https://cran.r-project.org/package=covr) - Version 3.6.5, published 2025-11-09
- [hedgehog CRAN](https://cran.r-project.org/package=hedgehog) - Version 0.2, published 2025-11-03
- [R Packages (2e) - Testing](https://r-pkgs.org/testing-basics.html) - Canonical guide
- [CRAN Submission Checklist](https://cran.r-project.org/web/packages/submission_checklist.html) - Official requirements
- [R eigen() documentation](https://stat.ethz.ch/R-manual/R-patched/library/base/html/eigen.html) - Eigenvector ambiguity

### Secondary (MEDIUM confidence)
- [rOpenSci Packaging Guide](https://devguide.ropensci.org/pkg_building.html) - Scientific R package best practices
- [Bioconductor Unit Tests](https://contributions.bioconductor.org/tests.html) - Domain-specific testing standards
- [Epiverse-TRACE: Statistical Correctness Testing](https://epiverse-trace.github.io/posts/statistical-correctness/) - Property-based testing for stats

### Tertiary (LOW confidence)
- [R Package Quality: Code Quality (R-bloggers)](https://www.r-bloggers.com/2025/07/r-package-quality-code-quality/) - Community practices
- [fpCompare package](https://cran.r-project.org/web/packages/fpCompare/vignettes/fpCompare.html) - Tolerance calibration reference

---
*Research completed: 2026-01-19*
*Ready for roadmap: yes*
