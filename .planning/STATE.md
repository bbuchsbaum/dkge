# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-01-19)

**Core value:** Implementation must faithfully match the mathematical specification in `algo.md` - users need to trust the results before publishing research based on them.
**Current focus:** Phase 4 - Numerical Edge Cases

## Current Position

Phase: 4 of 6 (Numerical Edge Cases)
Plan: 1 of 2 in current phase
Status: In progress
Last activity: 2026-01-20 - Completed 04-01-PLAN.md (numerical robustness)

Progress: [██████░░░░] 58%

## Performance Metrics

**Velocity:**
- Total plans completed: 7
- Average duration: 3.6 min
- Total execution time: 0.42 hours

**By Phase:**

| Phase | Plans | Total | Avg/Plan |
|-------|-------|-------|----------|
| 01-data-kernel-foundation | 2 | 8min | 4min |
| 02-fit-layer-correctness | 2 | 4min | 2min |
| 03-cross-fitting-validation | 2 | 7min | 3.5min |
| 04-numerical-edge-cases | 1 | 6min | 6min |

**Recent Trend:**
- Last 5 plans: 02-02 (2min), 03-01 (3min), 03-02 (4min), 04-01 (6min)
- Trend: Slightly increasing (edge case tests require more thorough verification)

*Updated after each plan completion*

## Accumulated Context

### Decisions

Decisions are logged in PROJECT.md Key Decisions table.
Recent decisions affecting current work:

- [Init]: Focus on core R/ only - `future/` is experimental, not exported
- [Init]: R CMD check as success gate - standard R package quality bar
- [Init]: Mathematical accuracy priority - publication requires trusted results
- [01-01]: NA/Inf values in beta matrices are accepted (not rejected) - documented behavior
- [01-01]: Effect order follows first subject's design column order
- [01-01]: Zero-cluster subjects accepted without error
- [01-01]: Duplicate effect names accepted - may need future hardening
- [01-02]: Circular kernels use l=0.5 in PSD tests (short length-scale ensures PSD)
- [01-02]: Use pmin instead of min for vectorized circular distance computation (bug fix)
- [02-01]: K-orthonormality tolerance: 1e-8 (consistent with Phase 1)
- [02-01]: Gram matrix verification tolerance: 1e-10 (tighter for exact linear algebra)
- [02-02]: SNR thresholds for recovery: 20->0.98, 8->0.90, 2->0.70
- [02-02]: Graceful degradation allows 0.05 tolerance for sampling variance
- [03-01]: Chat_minus tolerance 1e-7 due to floating-point accumulation
- [03-01]: Basis comparison via cosine similarity for sign-flip invariance
- [03-02]: Basis rotation ambiguity - use K-orthonormality check instead of per-column cosine
- [03-02]: Perturbation magnitude test requires small gaps + large off-diagonal coupling
- [04-01]: Minimum 2 subjects required for group analysis (error on single subject)
- [04-01]: Condition number threshold 1e8 for ill-conditioning warnings
- [04-01]: Sparse subject threshold >50% missing effects for warnings
- [04-01]: effective_rank and rank_reduced metadata stored in fit objects

### Pending Todos

None.

### Blockers/Concerns

- Current R CMD check status unknown - needs baseline run
- GitHub-only dependencies may complicate CRAN submission
- Pre-existing deprecation warning from multivarious::prep() - affects test output but not functionality

## Session Continuity

Last session: 2026-01-20
Stopped at: Completed 04-01-PLAN.md (numerical robustness infrastructure)
Resume file: None
