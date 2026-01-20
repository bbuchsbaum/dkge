# Phase 6: Integration + S3 Contracts - Context

**Gathered:** 2026-01-20
**Status:** Ready for planning

<domain>
## Phase Boundary

Verify end-to-end pipelines and user-facing API behave as documented. This phase ensures the package works as a complete product — from data in to results out — with proper R conventions (S3 methods, R CMD check, examples). All exported functions have working examples and test coverage exceeds 80%.

</domain>

<decisions>
## Implementation Decisions

### Pipeline scenarios
- Minimum test case: 2 subjects (smallest valid multi-subject case)
- Data dimensions: Both small synthetic (q=5, P=20) for fast CI tests AND realistic fMRI-like (q=10-20, P=1000+) for comprehensive suite
- Kernel types: Nominal only — most common case, sufficient for integration verification
- Pipeline scope: Full pipeline — fit + LOSO + transport + inference, all chained together

### S3 method contracts
- print(): Minimal summary — just dimensions, rank, n_subjects (one-liner style)
- predict(): Returns projected coordinates — matrix of q-space projections for new subjects
- as.data.frame(): Wide format — one row per subject, components as columns
- Test approach: Verify output correctness only — trust R's dispatch mechanics

### R CMD check strictness
- Target: 0 errors, 0 warnings, 0 notes (CRAN submission ready)
- Deprecation warning: Fix upstream — update multivarious dependency or work around the deprecated prep() call
- GitHub dependencies: Document in DESCRIPTION using Remotes: field
- Check mode: Standard R CMD check (not --as-cran)

### Example coverage
- Style: Minimal runnable — smallest code that demonstrates the function works
- Slow operations: Wrap in \donttest{} to keep R CMD check fast
- Coverage requirement: All exported functions must have working examples (no \dontrun{})
- Data source: Package datasets in data/ for reuse across examples

### Claude's Discretion
- Exact test data generation approach
- How to structure the package dataset(s)
- Which specific assertions for S3 method output verification
- Order of addressing R CMD check issues

</decisions>

<specifics>
## Specific Ideas

- Full pipeline test chains: dkge_data() → dkge_fit() → LOSO → transport → inference
- Package dataset should be small enough for fast examples but realistic enough to demonstrate the workflow
- Examples should be self-documenting — a user reading the example should understand what the function does

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope

</deferred>

---

*Phase: 06-integration-s3-contracts*
*Context gathered: 2026-01-20*
