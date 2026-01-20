# Codebase Concerns

**Analysis Date:** 2025-01-19

## Tech Debt

**Streaming Implementation Remains in `future/` Directory:**
- Issue: The streaming two-pass DKGE fit (`dkge_fit_streamed()`) is implemented but not integrated into the main package
- Files: `future/dkge-stream.R`
- Impact: Users cannot benefit from memory-light streaming estimators for large datasets; main `dkge_fit()` loads all subject data into memory
- Fix approach: Move streaming functions to `R/`, integrate with the `dkge_data` constructor, add tests and documentation

**Reliability Adaptive Weighting Not Implemented:**
- Issue: `dkge_weights(adapt = "reliability")` is declared but throws `stop("adapt = 'reliability' not implemented")`
- Files: `R/dkge-weights.R` (line 342)
- Impact: Documented feature that cannot be used; users expecting reliability-based weighting will get runtime errors
- Fix approach: Implement reliability estimation from subject betas or remove from the enum

**Freedman-Lane Inference Requires External Adapters:**
- Issue: `dkge_freedman_lane()` is a scaffold requiring user-supplied `fit_glm`, `resample_resid`, and `transport_and_stat` adapter functions
- Files: `R/dkge-inference.R` (lines 442-500)
- Impact: Time-series level permutation tests cannot be run without significant user effort to implement adapters
- Fix approach: Provide default adapters using `fmrireg` for GLM fitting and residual resampling

**Future Directory Contains Unintegrated Code:**
- Issue: Several R files in `future/` are not exported or tested: CPCA variants, CV helpers, mixed-effects, Sinkhorn alternatives, neuroim2 streaming
- Files: `future/dkge-cpca.R`, `future/dkge-cv.R`, `future/dkge-mixed.R`, `future/dkge-sinkhorn.R`, `future/dkge-sinkhorn-cpp.R`, `future/dkge-neuroim2-stream.R`, `future/dkge-write.R`
- Impact: Promising features sit unused; divergence between main code and future prototypes
- Fix approach: Evaluate each future module for production-readiness; integrate or archive

**Deprecated Arguments Still Present:**
- Issue: `transported` argument in `dkge_infer()` and `transport` argument in `dkge_contrast()` are deprecated but remain in function signatures
- Files: `R/dkge-inference.R` (line 75, 152), `R/dkge-contrast.R` (line 116)
- Impact: API clutter; users may rely on deprecated pathways
- Fix approach: Remove deprecated arguments after a suitable deprecation period; document migration in NEWS

## Known Bugs

**None currently documented as TODO/FIXME/BUG in source code.**

The codebase is notably clean of inline bug markers.

## Security Considerations

**Environment Manipulation for Cache:**
- Risk: `dkge_sinkhorn_cache` uses a module-level environment to cache dual variables across calls
- Files: `R/dkge-transport.R` (lines 47-72)
- Current mitigation: Cache is bounded to 64 entries and can be manually cleared via `dkge_clear_sinkhorn_cache()`
- Recommendations: Document cache behavior in user-facing functions; consider using `memoise` for a more standard approach

**Random Seed Side Effects:**
- Risk: Several functions call `set.seed()` directly when a `seed` argument is provided, affecting global RNG state
- Files: `R/dkge-classify.R` (line 77), `R/dkge-bootstrap.R` (line 41, 138), `R/dkge-inference.R` (line 468)
- Current mitigation: Seed argument is optional; users control when to set
- Recommendations: Use `withr::local_seed()` to avoid side effects on global RNG state

**No Input Sanitization for File Paths:**
- Risk: Functions like `dkge_write_group_map()` accept file paths without validation
- Files: `R/dkge-write.R`
- Current mitigation: Functions rely on underlying I/O libraries for path handling
- Recommendations: Add explicit path validation for user-supplied output paths

## Performance Bottlenecks

**Classification Module Builds Feature Matrices Per-Fold:**
- Problem: `dkge_classify()` recomputes `Z <- target$weight_matrix %*% loader$Y` for every fold and subject combination
- Files: `R/dkge-classify.R` (lines 370-382, 440-447)
- Cause: Feature cache is local to each `run_cv()` call; no cross-fold sharing
- Improvement path: Pre-compute features once per subject before cross-validation loop; store in shared cache

**Large File: dkge-classify.R (1307 lines):**
- Problem: The classification module is monolithic with multiple internal classifiers, target handling, and metric computation
- Files: `R/dkge-classify.R`
- Cause: Organic growth as features were added
- Improvement path: Factor out LDA/logit classifiers into separate module; extract target preparation; split cell vs delta mode logic

**Sinkhorn Transport Recomputed Per-Component:**
- Problem: `dkge_transport_loadings_to_medoid()` builds transport plans incrementally but still iterates per component
- Files: `R/dkge-transport.R` (lines 469-488)
- Cause: Historical API design; transport cache helps but is built lazily
- Improvement path: Use `dkge_prepare_transport()` proactively in pipelines; consider batched transport

**Eigendecomposition Repeated in Bootstrap:**
- Problem: `dkge_bootstrap_qspace()` performs full `eigen()` for each of B replicates
- Files: `R/dkge-bootstrap.R` (line 186)
- Cause: Multiplier weights change Chat each iteration; no incremental update
- Improvement path: Implement rank-one updates for small perturbations; batch eigensolves

## Fragile Areas

**Effect Alignment Across Subjects:**
- Files: `R/dkge-data.R` (lines 313-350), `R/dkge-align-effects.R`
- Why fragile: Subject betas must have matching row names; misalignment triggers complex union-alignment logic that can mask missing data
- Safe modification: Always validate effect names before calling `dkge_data()`; use explicit `subject_ids` parameter
- Test coverage: `tests/testthat/test-data.R` covers basic cases but not all edge cases around partial overlap

**Analytic LOSO Fallback Logic:**
- Files: `R/dkge-analytic.R` (lines 59-150)
- Why fragile: Multiple conditions trigger fallback to full eigensolve (solver type, voxel weights, eigengaps); any change to fit structure can break detection
- Safe modification: Check `result$method` to verify which path was taken; add new conditions at top of function
- Test coverage: `tests/testthat/test-analytic.R` covers main paths

**Kernel Root Computation:**
- Files: `R/design-kernel.R` (lines 236-258)
- Why fragile: Eigenvalue clamping with warning when >1% are clamped; symmetrization warning for non-symmetric input
- Safe modification: Always pass symmetric kernels; monitor warnings for clamping events
- Test coverage: `tests/testthat/test-design-kernel.R` covers basic kernels

**Target Weight Matrix Construction:**
- Files: `R/dkge-targets.R`
- Why fragile: Complex logic for residualization, factor collapse, and interaction terms; depends on `kernel_info$map` and `kernel_info$levels`
- Safe modification: Use helper functions like `dkge_anchor_targets_from_prototypes()` for anchor fits; validate targets before classification
- Test coverage: Partial; `tests/testthat/test-classify.R` covers some paths

## Scaling Limits

**In-Memory Beta Storage:**
- Current capacity: All subject betas loaded into `fit$Btil` list
- Limit: Memory proportional to S * q * max(P_s) where S = subjects, q = effects, P_s = clusters
- Scaling path: Implement streaming from `future/dkge-stream.R`; use disk-backed storage via neuroim2 on-disk arrays

**Transport Plan Memory:**
- Current capacity: Full Q x Q transport matrices stored per subject
- Limit: Quadratic in reference parcellation size
- Scaling path: Sparse transport plans; approximate transport with fewer iterations; discard plans after use

## Dependencies at Risk

**GitHub-Only Dependencies:**
- Risk: `neuroim2`, `fmridesign`, `fmriAR`, `fmrireg` are installed from GitHub (bbuchsbaum repos)
- Impact: Package installation fails if GitHub rate-limited or repos reorganized
- Migration plan: Publish dependencies to CRAN or R-universe; add fallback installation instructions

**Rcpp/RcppArmadillo for C++ Code:**
- Risk: C++ compilation required at install time
- Impact: Users without compiler toolchain cannot install; OpenMP optional but affects performance
- Migration plan: Provide pre-compiled binaries via r-universe; document toolchain requirements

## Missing Critical Features

**No Parallel Support in Classification:**
- Problem: `dkge_classify()` has `parallel = FALSE` parameter but parallelism is not implemented
- Files: `R/dkge-classify.R` (line 69)
- Blocks: Large-scale permutation testing is slow; users cannot leverage multiple cores

**No Progress Reporting:**
- Problem: Long-running operations (bootstrap, permutation tests) provide no progress feedback
- Files: `R/dkge-bootstrap.R`, `R/dkge-inference.R`, `R/dkge-classify.R`
- Blocks: Users cannot estimate completion time for large analyses

**No Checkpoint/Resume for Long Analyses:**
- Problem: Bootstrap and permutation loops cannot be resumed after interruption
- Files: `R/dkge-bootstrap.R`, `R/dkge-inference.R`
- Blocks: Analyses must restart from scratch if interrupted

## Test Coverage Gaps

**No Tests for `future/` Directory:**
- What's not tested: Streaming fit, mixed effects, advanced CV, alternative Sinkhorn implementations
- Files: All files in `future/`
- Risk: Future code may have regressions when integrated
- Priority: Low (not exported)

**Limited Transport Edge Cases:**
- What's not tested: Transport with very different cluster counts across subjects; degenerate parcellations
- Files: `R/dkge-transport.R`
- Risk: Runtime errors or numerical instability with edge-case inputs
- Priority: Medium

**Sparse Inference Testing:**
- What's not tested: Freedman-Lane permutations (only scaffold exists); interaction between transport and inference
- Files: `R/dkge-inference.R`
- Risk: Users may attempt unsupported workflows
- Priority: Medium

**Missing Vignette Coverage for Advanced Features:**
- What's not tested: Anchor-based fits, validated contrasts, custom mappers
- Files: `vignettes/`
- Risk: Features exist but users cannot discover or learn them
- Priority: Medium

---

*Concerns audit: 2025-01-19*
