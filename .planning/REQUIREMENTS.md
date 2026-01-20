# Requirements: dkge Publication Readiness

## Overview

Quality assurance requirements for the dkge R package to achieve publication readiness.

## v1 Requirements (Active)

### Mathematical Correctness

| ID | Requirement | Priority |
|----|-------------|----------|
| MATH-01 | Mathematical accuracy verified - algorithms match `algo.md` specification | Must |

### API Quality

| ID | Requirement | Priority |
|----|-------------|----------|
| API-01 | API contracts verified - all exported functions behave as documented | Must |
| CHECK-01 | R CMD check passes - zero errors, zero warnings, zero notes | Must |

### Test Coverage

| ID | Requirement | Priority |
|----|-------------|----------|
| COV-01 | Higher test coverage - more code paths exercised (target 80%+ on exports) | Must |
| PROP-01 | Property-based tests - mathematical invariants verified (K-orthonormality, symmetry, PSD) | Should |
| REG-01 | Regression tests - known-good outputs locked in to catch drift | Should |

### Robustness

| ID | Requirement | Priority |
|----|-------------|----------|
| EDGE-01 | Edge case coverage - degenerate inputs, boundary conditions, ill-conditioned matrices | Must |
| FRAG-01 | Fragile areas hardened - effect alignment, analytic LOSO fallback, kernel roots, target weights | Must |

## v2 Requirements (Deferred)

| ID | Requirement | Rationale for Deferral |
|----|-------------|------------------------|
| PERF-01 | Performance optimization | Correctness first, speed later |
| STREAM-01 | Streaming implementation | Valuable but not required for publication |
| FUTURE-01 | `future/` directory experiments | Not exported, defer to post-publication |

## Traceability

| Requirement | Phase | Status |
|-------------|-------|--------|
| MATH-01 | 1, 2, 3, 5 | Pending |
| API-01 | 1, 6 | Pending |
| CHECK-01 | 6 | Pending |
| COV-01 | 1, 3, 6 | Pending |
| PROP-01 | 2, 5 | Pending |
| REG-01 | 2 | Pending |
| EDGE-01 | 4 | Pending |
| FRAG-01 | 1, 3, 4 | Pending |

---
*Created: 2026-01-19*
