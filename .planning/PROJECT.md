# dkge Publication Readiness

## What This Is

A quality assurance initiative for the dkge R package — ensuring mathematical correctness, API contract compliance, and comprehensive test coverage before publication. The package implements Design-Kernel Group Embedding for cluster-level fMRI analysis.

## Core Value

The implementation must faithfully match the mathematical specification in `algo.md` — users need to trust the results before publishing research based on them.

## Requirements

### Validated

Existing functionality that works (from codebase analysis):

- ✓ Data layer: `dkge_subject()`, `dkge_data()` bundle multi-subject betas — existing
- ✓ Kernel layer: `design_kernel()` constructs factorial similarity kernels — existing
- ✓ Fit layer: `dkge_fit()` estimates shared latent basis U via eigendecomposition — existing
- ✓ Contrast layer: `dkge_contrast()` computes LOSO/K-fold/analytic cross-fitted contrasts — existing
- ✓ Transport layer: `dkge_transport_contrasts_to_medoid()` aligns subjects via Sinkhorn OT — existing
- ✓ Inference layer: `dkge_signflip_maxT()` provides permutation-based FWER control — existing
- ✓ Prediction layer: `dkge_freeze()`, `predict.dkge()` apply frozen basis to new subjects — existing
- ✓ Classification layer: `dkge_classify()` performs cross-validated decoding — existing
- ✓ Pipeline layer: `dkge_pipeline()` orchestrates full workflows — existing

### Active

Quality assurance goals for publication:

- [ ] Mathematical accuracy verified — algorithms match `algo.md` specification
- [ ] API contracts verified — all exported functions behave as documented
- [ ] R CMD check passes — zero errors, zero warnings, zero notes
- [ ] Higher test coverage — more code paths exercised
- [ ] Property-based tests — mathematical invariants verified (orthonormality, symmetry, K-orthogonality)
- [ ] Regression tests — known-good outputs locked in to catch drift
- [ ] Edge case coverage — degenerate inputs, boundary conditions, ill-conditioned matrices
- [ ] Fragile areas hardened — effect alignment, analytic LOSO fallback, kernel roots, target weights

### Out of Scope

- `future/` directory experiments — not exported, defer to post-publication
- New features — focus is verification of existing functionality
- Performance optimization — correctness first, speed later
- Streaming implementation — valuable but not required for publication

## Context

**Publication target:** Method paper describing DKGE for fMRI community

**Existing test infrastructure:**
- testthat edition 3
- 67 test files in `tests/testthat/`
- Helper files for fixtures and toy data
- Some coverage gaps identified in CONCERNS.md

**Known fragile areas (from codebase analysis):**
- Effect alignment across subjects (`R/dkge-data.R`)
- Analytic LOSO fallback logic (`R/dkge-analytic.R`)
- Kernel root computation (`R/design-kernel.R`)
- Target weight matrix construction (`R/dkge-targets.R`)

**Mathematical specification:** `data-raw/algo.md` defines the algorithm

**Current R CMD check status:** Unknown — needs to be run

## Constraints

- **Scope**: Core exported package only (`R/` directory) — `future/` is out of scope
- **Language**: R package conventions (roxygen2, testthat, R CMD check)
- **Dependencies**: GitHub-only deps (neuroim2, fmridesign, fmrireg) may complicate CRAN submission

## Key Decisions

| Decision | Rationale | Outcome |
|----------|-----------|---------|
| Focus on core R/ only | `future/` is experimental, not exported | — Pending |
| R CMD check as success gate | Standard R package quality bar | — Pending |
| Mathematical accuracy priority | Publication requires trusted results | — Pending |

---
*Last updated: 2026-01-19 after initialization*
