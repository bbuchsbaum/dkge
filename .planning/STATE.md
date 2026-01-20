# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-01-19)

**Core value:** Implementation must faithfully match the mathematical specification in `algo.md` - users need to trust the results before publishing research based on them.
**Current focus:** Phase 6 - Integration + S3 Contracts

## Current Position

Phase: 6 of 6 (Integration + S3 Contracts)
Plan: 0 of 3 in current phase
Status: Ready to plan
Last activity: 2026-01-20 - Phase 5 verified and complete

Progress: [████████░░] 83%

## Performance Metrics

**Velocity:**
- Total plans completed: 10
- Average duration: 4.3 min
- Total execution time: 0.72 hours

**By Phase:**

| Phase | Plans | Total | Avg/Plan |
|-------|-------|-------|----------|
| 01-data-kernel-foundation | 2 | 8min | 4min |
| 02-fit-layer-correctness | 2 | 4min | 2min |
| 03-cross-fitting-validation | 2 | 7min | 3.5min |
| 04-numerical-edge-cases | 2 | 11min | 5.5min |
| 05-transport-inference | 2 | 11min | 5.5min |

**Recent Trend:**
- Last 5 plans: 04-01 (6min), 04-02 (5min), 05-01 (7min), 05-02 (6min)
- Trend: Stable at ~5-7min per plan for testing phases

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
- [04-02]: Use 5 diverse seeds for single-rank tests: {1, 42, 123, 999, 2024}
- [04-02]: Use nominal kernels for multi-rank recovery (better conditioned than ordinal)
- [04-02]: CV threshold 10% for recovery stability, 15-25% for multi-rank
- [05-01]: Sinkhorn tolerance 1e-6 for marginal convergence tests
- [05-01]: R/CPP equivalence at 1e-7 tolerance (accounts for floating-point accumulation)
- [05-01]: Non-convergence returns best-effort result without warning (documented actual behavior)
- [05-01]: Near-diagonal threshold: >0.8-0.9 diagonal mass for similar embeddings
- [05-02]: Test FWER control instead of per-cluster uniformity for max-T procedure

### Pending Todos

None.

### Blockers/Concerns

- Current R CMD check status unknown - needs baseline run
- GitHub-only dependencies may complicate CRAN submission
- Pre-existing deprecation warning from multivarious::prep() - affects test output but not functionality

## Session Continuity

Last session: 2026-01-20
Stopped at: Phase 5 complete, ready for Phase 6
Resume file: None
