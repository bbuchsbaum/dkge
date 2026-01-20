# Phase 5: Transport + Inference - Context

**Gathered:** 2026-01-20
**Status:** Ready for planning

<domain>
## Phase Boundary

Verify Sinkhorn optimal transport convergence and permutation test null calibration. Tests confirm the transport algorithm produces valid doubly-stochastic matrices and that permutation-based inference yields calibrated p-values under the null hypothesis.

</domain>

<decisions>
## Implementation Decisions

### Convergence tolerances
- Doubly-stochastic tolerance: 1e-6 (standard numerical tolerance for iterative algorithms)
- Verify final result only, not convergence progress during iterations
- Non-convergence behavior: warning + best result (not error)
- Tests should verify warning is issued when max iterations reached without convergence

### Transport verification
- Verify doubly-stochastic property: all row sums = 1, all column sums = 1, all entries >= 0
- Also verify sparsity pattern for near-diagonal/identity transport cases
- Check that transport plan has expected structure when source and target are similar

### Permutation test design
- Use 500 permutations for null calibration tests (balanced precision vs runtime)
- Null scenario constructed via shuffled labels (not pure noise) - break true association in real-structured data
- Single seed for reproducibility, trust the algorithm

### Claude's Discretion
- Chi-square significance threshold for uniformity testing
- Specific sparsity thresholds for near-diagonal transport verification
- Exact construction of shuffled-label test scenarios
- Parallel vs sequential equivalence tolerance (not discussed, but in success criteria)
- Deterministic edge case selection (identity mappings)

</decisions>

<specifics>
## Specific Ideas

- Doubly-stochastic verification should match the 1e-6 tolerance used by the Sinkhorn algorithm itself
- Shuffled labels approach mirrors standard permutation testing practice in neuroimaging
- Warning-based non-convergence handling allows graceful degradation rather than hard failure

</specifics>

<deferred>
## Deferred Ideas

None - discussion stayed within phase scope

</deferred>

---

*Phase: 05-transport-inference*
*Context gathered: 2026-01-20*
