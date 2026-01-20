# Phase 4: Numerical Edge Cases - Context

**Gathered:** 2026-01-19
**Status:** Ready for planning

<domain>
## Phase Boundary

Verify the package handles degenerate inputs gracefully without silent failures. This includes rank-deficient matrices, near-singular conditions, NaN propagation prevention, multi-seed robustness, and partial effect overlap between subjects.

</domain>

<decisions>
## Implementation Decisions

### Error vs Degradation Policy
- Rank-deficient matrices: warn and reduce rank (not error, not silent)
- Warnings must identify culprits: include subject IDs and/or effect names
- Minimum viable rank threshold required: error if rank drops below threshold (Claude determines reasonable threshold based on requested rank)
- Detection timing: early (at data construction in dkge_subject()/dkge_data())
- Single problematic subject: drop and continue (with warning), don't fail entire fit
- Minimum subjects for group analysis: at least 2 subjects required
- Result metadata: store degradation info (dropped subjects, reduced rank) in returned object
- LOSO with problematic held-out subject: skip that fold (return NA), others continue
- Distinguish ridge regularization from natural rank deficiency in messaging
- No strict mode flag: single consistent behavior path
- Default stance for unanticipated edge cases: fail safe (error)
- Test scope: verify behavior only, not message wording

### NaN/Inf Handling Strategy
- Input NA/NaN in beta matrices: exclude affected voxels from computation
- Runtime NaN/Inf (from computation): warn and continue with degraded results
- No proactive checks: handle after operation, detect NaN in results
- Exclusion reporting: store indices/count of excluded voxels in result object

### Tolerance Thresholds
- Near-singular condition number threshold: 1e8 (moderate)
- Existing tolerances (1e-8 K-orthonormality, 1e-10 Gram): fixed in code, not configurable

### Partial Overlap Behavior
- Effect union handling: Claude decides most sensible approach for reasonable results
- Sparse subjects (>50% missing effects): warn above threshold but include
- Design kernel K with missing effects: Claude decides approach that preserves K's mathematical properties
- Effect presence tracking: store subject × effect presence matrix in result object

### Claude's Discretion
- Exact minimum viable rank threshold heuristic
- Eigenvalue "effectively zero" threshold determination
- Test assertion tolerance decisions per scenario
- Best approach for effect union that produces reasonable results
- Kernel handling for partially-missing effects (preserving PSD property)

</decisions>

<specifics>
## Specific Ideas

- Mathematical accuracy for publication is the core value - users need to trust results before publishing research
- Behavior verification (not message wording) is sufficient for test coverage
- Earlier phases established tolerances: 1e-8 for K-orthonormality, 1e-10 for Gram matrices, 1e-7 for Chat_minus floating-point accumulation

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope

</deferred>

---

*Phase: 04-numerical-edge-cases*
*Context gathered: 2026-01-19*
