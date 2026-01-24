# dkge Publication Readiness

## What This Is

A quality assurance initiative for the dkge R package — ensuring mathematical correctness, API contract compliance, and comprehensive test coverage before publication. The package implements Design-Kernel Group Embedding for cluster-level fMRI analysis. Now publication-ready with verified mathematical correctness, 1382 passing tests, and R CMD check compliance.

## Core Value

The implementation must faithfully match the mathematical specification in `algo.md` — users need to trust the results before publishing research based on them.

## Requirements

### Validated

- Data layer: `dkge_subject()`, `dkge_data()` bundle multi-subject betas — v1.0
- Kernel layer: `design_kernel()` constructs factorial similarity kernels — v1.0
- Fit layer: `dkge_fit()` estimates shared latent basis U via eigendecomposition — v1.0
- Contrast layer: `dkge_contrast()` computes LOSO/K-fold/analytic cross-fitted contrasts — v1.0
- Transport layer: `dkge_transport_contrasts_to_medoid()` aligns subjects via Sinkhorn OT — v1.0
- Inference layer: `dkge_signflip_maxT()` provides permutation-based FWER control — v1.0
- Prediction layer: `dkge_freeze()`, `predict.dkge()` apply frozen basis to new subjects — v1.0
- Classification layer: `dkge_classify()` performs cross-validated decoding — v1.0
- Pipeline layer: `dkge_pipeline()` orchestrates full workflows — v1.0
- MATH-01: Mathematical accuracy verified — algorithms match `algo.md` specification — v1.0
- API-01: API contracts verified — all exported functions behave as documented — v1.0
- CHECK-01: R CMD check passes — zero errors, zero warnings (except system), zero notes — v1.0
- COV-01: Higher test coverage — 1382 tests, 80 test files — v1.0
- PROP-01: Property-based tests — K-orthonormality, symmetry, PSD verified — v1.0
- REG-01: Regression tests — known-good outputs locked in — v1.0
- EDGE-01: Edge case coverage — degenerate inputs, boundary conditions handled — v1.0
- FRAG-01: Fragile areas hardened — effect alignment, analytic LOSO, kernel roots — v1.0

### Active

- [ ] CRAN submission workflow — coordinate with dependency maintainers
- [ ] Performance optimization — post-publication priority
- [ ] Streaming implementation — for very large datasets

### Out of Scope

- `future/` directory experiments — not exported, defer to post-publication
- Mobile/web interfaces — R package only

## Context

**Shipped v1.0** with comprehensive test hardening:
- 1382 tests passing across 80 test files
- 53 R source files, ~55K lines of R code
- R CMD check: 0 errors, 1 WARNING (system only), 0 notes
- All 14 vignettes building successfully
- 43/43 exported functions with @examples

**Tech stack:** R, testthat3, roxygen2, knitr, Rcpp

**Dependencies:** neuroim2, fmridesign, fmrireg (GitHub-only — may require CRAN coordination)

**Known issues:**
- GitHub-only dependencies complicate CRAN submission
- Deprecation warning from multivarious::prep() — external dependency
- covr instrumentation fails on Rcpp packages

## Constraints

- **Scope**: Core exported package only (`R/` directory) — `future/` is out of scope
- **Language**: R package conventions (roxygen2, testthat, R CMD check)
- **Dependencies**: GitHub-only deps may require CRAN coordination

## Key Decisions

| Decision | Rationale | Outcome |
|----------|-----------|---------|
| Focus on core R/ only | `future/` is experimental, not exported | Good |
| R CMD check as success gate | Standard R package quality bar | Good |
| Mathematical accuracy priority | Publication requires trusted results | Good |
| K-orthonormality tolerance 1e-8 | Consistent numerical precision | Good |
| Minimum 2 subjects for group analysis | Single subject has no group structure | Good |
| Condition number threshold 1e8 | Warns on ill-conditioned matrices | Good |
| Accept NA/Inf in beta matrices | Documented behavior, downstream handling | Pending review |
| Accept duplicate effect names | May need future hardening | Pending review |
| Replace \dontrun with \donttest | Slow examples still run in --as-cran | Good |

---
*Last updated: 2026-01-22 after v1.0 milestone*
