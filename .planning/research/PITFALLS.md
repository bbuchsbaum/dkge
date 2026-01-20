# Domain Pitfalls: Testing Scientific/Numerical R Packages

**Domain:** Scientific R package testing (matrix decomposition, statistical inference, numerical algorithms)
**Package:** dkge - Design-Kernel Group Embedding for cluster-level fMRI analysis
**Researched:** 2026-01-19

---

## Critical Pitfalls

Mistakes that cause false confidence in correctness or missed bugs that escape to publication.

### Pitfall 1: Eigenvector Sign and Ordering Ambiguity

**What goes wrong:** Tests fail intermittently or pass when they should fail because eigenvector signs and orderings are not mathematically unique. Two correct implementations can return eigenvectors that differ by sign flips or permutations when eigenvalues are close.

**Why it happens:** R's `eigen()` function explicitly documents that "eigenvectors may differ in sign and (in the asymmetric case) in normalization. They may also differ between methods and between platforms." Tests that directly compare eigenvector values will fail across platforms or R versions despite both being mathematically correct.

**Consequences:**
- CRAN check failures on different architectures (M1 Mac vs x86 Linux)
- False test failures that erode trust in the test suite
- Alternatively, tests that pass despite incorrect implementations

**Warning signs:**
- Tests pass locally but fail on GitHub Actions or CRAN
- Adding `set.seed()` fixes a test temporarily
- Tests comparing eigenvectors with `expect_equal()` directly

**Prevention:**
1. Never compare eigenvectors directly; compare derived invariants:
   - Subspace angles (principal angles between column spaces)
   - Reconstruction error: `norm(A - U %*% diag(d) %*% t(V))`
   - Gram matrices: `t(U) %*% K %*% U` should equal identity
   - Projection matrices: `U %*% t(U)` is sign-invariant
2. Use absolute values when comparing loadings: `abs(cor(v1, v2)) > 0.99`
3. For ordered comparisons, canonicalize sign by requiring positive first element
4. Test K-orthonormality constraint rather than specific vectors

**Phase mapping:** Phase 1 (Core fit correctness) - establish invariant-based testing patterns early

**DKGE-specific context:** The `dkge_procrustes_K()` and `dkge_align_bases_K()` functions exist precisely to handle this. Tests should verify Procrustes alignment quality (cosines near 1) rather than raw basis equality. The existing test at `test-contrast.R:119` checking `gram <- t(U_s) %*% K %*% U_s; expect_equal(gram, diag(r), tolerance = 1e-10)` is the correct pattern.

**Sources:**
- [R eigen() documentation](https://stat.ethz.ch/R-manual/R-patched/library/base/html/eigen.html)
- [RSpectra vignette on eigenvalue decomposition](https://cran.r-project.org/web/packages/RSpectra/vignettes/introduction.html)

---

### Pitfall 2: Tolerance Calibration Blindness

**What goes wrong:** Tests use arbitrary tolerances (1e-6, 1e-10) without understanding what tolerance is actually achievable for the algorithm. Tests either pass when they should fail (tolerance too loose) or fail on correct code (tolerance too tight).

**Why it happens:** Developers copy tolerance values from examples without considering:
- Machine epsilon propagation through matrix operations
- Condition number of the input matrices
- Accumulation of rounding errors in iterative algorithms
- Whether the algorithm is forward-stable, backward-stable, or neither

**Consequences:**
- Bugs in numerical code go undetected for years
- Correct implementations fail tests on edge cases
- False sense of precision verification

**Warning signs:**
- Tests with magic tolerance numbers (1e-6) without justification
- Different tolerances for the same operation in different tests
- Tests that "just pass" with current tolerance
- Using `expect_identical()` for floating-point results

**Prevention:**
1. **Derive tolerances from theory:** For eigendecomposition of a well-conditioned matrix, expect ~`condition_number * .Machine$double.eps` relative error
2. **Use relative tolerances scaled to the problem:**
   ```r
   expect_equal(result, expected, tolerance = 1e-10 * max(abs(expected), 1))
   ```
3. **Document tolerance rationale in test comments:**
   ```r
   # Tolerance: 1e-8 allows for O(n^2) float ops on well-conditioned input
   expect_equal(reconstructed, K, tolerance = 1e-8)
   ```
4. **Test tolerance boundaries:** Verify that tightening tolerance by 10x causes expected failures
5. **Use `testthat::testthat_tolerance()` (default sqrt(.Machine$double.eps)) as baseline

**Phase mapping:** Phase 1 (Core fit correctness) - calibrate tolerances during initial test development

**DKGE-specific context:** The `.dkge_kernel_roots()` function adds jitter of 1e-10 to eigenvalues; tests should expect at least this much numerical noise. The analytic approximation in `dkge_analytic_loso()` uses `gap_tol = max(tol, 1e-8)` and `perturb_tol = 0.1` - these thresholds should be verified in tests.

**Sources:**
- [fpCompare package](https://cran.r-project.org/web/packages/fpCompare/vignettes/fpCompare.html) - default tolerance is `.Machine$double.eps^0.5`
- [Numerical Stability in Eigenvalue Decomposition](https://www.numberanalytics.com/blog/numerical-stability-eigenvalue-decomposition-techniques-best-practices)

---

### Pitfall 3: Testing the Wrong Invariant

**What goes wrong:** Tests verify properties that hold for correct code but also hold for many incorrect implementations. The test provides false confidence.

**Why it happens:** Developers test what's easy to test (dimensions, types, non-null values) rather than what matters (mathematical correctness). Tests verify necessary conditions but not sufficient conditions.

**Consequences:**
- Major algorithmic bugs escape testing
- Refactoring introduces silent regressions
- Publication of incorrect results

**Warning signs:**
- Tests that only check `expect_equal(dim(result), c(n, k))`
- Tests that verify output is finite without checking values
- No tests that would fail if core algorithm returned random noise
- All tests pass when a key computation is commented out

**Prevention:**
1. **Test mathematical identities that uniquely constrain the solution:**
   - For DKGE: `t(U) %*% K %*% U == I` (K-orthonormality)
   - For eigendecomposition: `A %*% V == V %*% diag(d)` (eigenvalue equation)
   - For Procrustes: `t(R) %*% R == I` and optimality condition
2. **Use recovery tests with known ground truth:**
   ```r
   # Generate data FROM the model, verify we recover parameters
   U_true <- known_orthonormal_basis()
   data <- generate_from_model(U_true)
   fit <- dkge_fit(data)
   subspace_distance(fit$U, U_true) < tolerance
   ```
3. **Test failure modes:** Verify code fails appropriately on invalid input
4. **Mutation testing:** Temporarily break core algorithm, verify tests fail

**Phase mapping:** Phase 1 (Core fit correctness) - establish recovery tests; Phase 2 (LOSO/cross-fitting) - test bias properties

**DKGE-specific context:** The test at `test-toy-recovery.R` demonstrates the correct pattern: generate data with known `U_true`, fit model, verify recovery via `dkge_cosines_K()`. This should be the template for all core algorithm tests.

---

### Pitfall 4: Data Leakage in Cross-Validation Tests

**What goes wrong:** Tests for leave-one-out or k-fold methods don't actually verify that held-out data is excluded. The test passes but the implementation leaks information.

**Why it happens:** Cross-validation correctness is subtle. Common leakage sources:
- Preprocessing (centering, scaling) computed on full data before split
- Feature selection using all data
- Hyperparameter tuning that peeks at test fold
- Pooled statistics (like the design Cholesky R) incorrectly computed

**Consequences:**
- LOSO contrasts are biased (the primary scientific claim of the method)
- Overly optimistic performance estimates
- Non-reproducible findings in publications

**Warning signs:**
- LOSO results identical to in-sample results
- Cross-validation improves results (should generally be worse or equal)
- No explicit test that held-out subject data is unused in basis computation

**Prevention:**
1. **Verify basis differs:** For each held-out subject s, check that the LOSO basis U^{-s} differs from the full basis U:
   ```r
   # Bases should be different (not identical to full U)
   expect_false(isTRUE(all.equal(U_loso, fit$U)))
   ```
2. **Test prediction orthogonality:** Held-out predictions should be orthogonal to training-only derived statistics
3. **Inject detectable leakage:** Add extreme values to held-out subject, verify they don't affect held-out basis
4. **Compare to oracle:** Compare LOSO basis to basis fitted without that subject's data entirely

**Phase mapping:** Phase 2 (LOSO/cross-fitting correctness) - critical for scientific validity

**DKGE-specific context:** The test at `test-contrast.R:109-120` verifies bases differ and are K-orthonormal. But it doesn't verify the *specific* property that subject s's data wasn't used. Need explicit test that: (1) recomputes Chat without subject s, (2) recomputes eigen, (3) verifies LOSO basis matches.

**Sources:**
- [Data Leakage in Cross-Validation](https://scikit-learn.org/stable/common_pitfalls.html)
- [Being Aware of Data Leakage and Cross-Validation Scaling](https://analyticalsciencejournals.onlinelibrary.wiley.com/doi/full/10.1002/cem.70026)

---

### Pitfall 5: Analytic Approximation Fallback Not Tested

**What goes wrong:** The analytic LOSO approximation has fallback conditions (eigengap too small, perturbation too large) but tests don't verify these conditions trigger correctly or that fallback produces correct results.

**Why it happens:** Approximation algorithms have complex validity conditions. Developers test the happy path (approximation works) but not the edge cases where fallback is needed.

**Consequences:**
- Approximation used when it's inaccurate
- Fallback never triggered despite invalid conditions
- Subtle numerical errors in edge cases

**Warning signs:**
- No tests where `method == "fallback"` is expected
- Fallback rate always 0% in test suite
- No tests with degenerate eigenvalue gaps

**Prevention:**
1. **Test fallback triggers:** Create inputs that should trigger each fallback condition:
   ```r
   # Create data with near-degenerate eigenvalues
   fit_degen <- create_degenerate_eigenvalue_case()
   result <- dkge_analytic_loso(fit_degen, s=1, c)
   expect_equal(result$method, "fallback")
   expect_equal(result$diagnostic$reason, "eigengap")
   ```
2. **Test fallback correctness:** Verify fallback gives same result as explicit recomputation
3. **Test boundary conditions:** Tolerance thresholds, epsilon values
4. **Verify diagnostic metadata:** Check that diagnostic info is populated correctly

**Phase mapping:** Phase 2 (LOSO/cross-fitting correctness) - test analytic approximation boundaries

**DKGE-specific context:** The code in `dkge-analytic.R` has multiple fallback paths: `solver_not_pooled`, `nonuniform_voxel_weights`, `missing_full_decomposition`, `dimension_mismatch`, `eigengap`, `perturbation_magnitude`. Each needs a test. Current `test-analytic.R` tests some but not all.

---

## Moderate Pitfalls

Mistakes that cause technical debt or delayed bug discovery.

### Pitfall 6: Seed-Dependent Tests Without Robustness Verification

**What goes wrong:** Tests pass with `set.seed(123)` but fail with different seeds. The test covers one lucky case, not general correctness.

**Why it happens:** Random test data generators can produce degenerate cases (rank-deficient matrices, extreme values) that break algorithms. Developers tune seed until tests pass rather than fixing underlying issues.

**Consequences:**
- Bugs discovered only when users encounter unlucky random states
- Fragile test suite that breaks on R updates
- False confidence in algorithm robustness

**Warning signs:**
- Changing `set.seed()` value breaks tests
- Tests fail "randomly" in CI
- Helper functions like `make_fit_fixture()` always use same seed

**Prevention:**
1. **Multi-seed testing:** Run critical tests across multiple seeds:
   ```r
   for (seed in c(1, 42, 123, 999, 2024)) {
     set.seed(seed)
     data <- generate_test_data()
     expect_true(verify_property(data))
   }
   ```
2. **Property-based testing:** Use packages like `hedgehog` or `quickcheck` for automatic input generation
3. **Explicit edge case generation:** Separately test degenerate cases rather than hoping random generation hits them
4. **Document seed sensitivity:** Note which tests are seed-sensitive and why

**Phase mapping:** Phase 4 (Numerical edge cases) - systematic robustness testing

**DKGE-specific context:** Tests use consistent seeds (123, 456, 789) but no multi-seed verification. The `dkge_sim_toy()` function should be tested across seeds.

---

### Pitfall 7: Missing Rank-Deficient/Degenerate Input Tests

**What goes wrong:** Algorithm works on well-conditioned inputs but fails silently or produces garbage on rank-deficient matrices, zero inputs, or extreme values.

**Why it happens:** Test data generators produce "nice" matrices. Real fMRI data can have:
- Collinear design columns
- Zero-variance voxels
- Extreme outliers
- Missing data patterns that create rank deficiency

**Consequences:**
- Runtime errors in production
- Silent numerical instability
- NaN/Inf propagation

**Warning signs:**
- No tests with `matrix(0, ...)` or near-zero inputs
- No tests with rank < expected dimensions
- No tests with Inf or NaN handling
- Tests always use orthogonalized design matrices

**Prevention:**
1. **Zero input tests:**
   ```r
   expect_no_error(dkge_fit(zero_betas, designs, K))
   # or expect specific error
   expect_error(dkge_fit(zero_betas, designs, K), "rank deficient")
   ```
2. **Near-singular tests:** Matrices with condition number > 1e10
3. **Extreme value tests:** Include outliers 1e6 times larger than typical values
4. **NaN propagation tests:** Verify NaN in one voxel doesn't corrupt all results
5. **Partial overlap tests:** Subjects with different effect sets (per `dkge_align_effects`)

**Phase mapping:** Phase 4 (Numerical edge cases) - dedicated edge case phase

**DKGE-specific context:** The `dkge_align_effects.R` code handles partial effect overlap, which creates rank-deficient completed kernels. Test `test-align-effects.R` should cover: (1) all subjects observe all effects, (2) minimal overlap, (3) one subject with no overlap, (4) prior-based completion.

---

### Pitfall 8: S3 Method Dispatch Testing Gap

**What goes wrong:** S3 generic methods (`print`, `predict`, `as.matrix`) work for one class but fail for subclasses or when called with unexpected argument combinations.

**Why it happens:** S3 dispatch is implicit. Developers test the primary path but not:
- Subclass variations
- Missing optional arguments
- Non-standard argument orders
- Conflicting method signatures

**Consequences:**
- User-facing functions fail in edge cases
- `as.data.frame()` works but `as.data.frame(stringsAsFactors=TRUE)` fails
- `predict()` with new data format causes cryptic errors

**Warning signs:**
- Methods only tested with minimal arguments
- No tests of inheritance (`inherits(x, "parent_class")`)
- Print methods tested only for non-error, not output content

**Prevention:**
1. **Test all exported S3 methods explicitly:**
   ```r
   test_that("as.data.frame.dkge_contrasts handles all arguments", {
     result <- as.data.frame(contrast_obj,
                            row.names = c("a", "b"),
                            optional = TRUE,
                            stringsAsFactors = FALSE)
     expect_s3_class(result, "data.frame")
   })
   ```
2. **Test subclass behavior:** If `dkge` inherits from `multiblock_biprojector`, test both dispatch paths
3. **Verify print output:** Use `capture.output()` and check for expected strings
4. **Test error messages:** Verify user-friendly errors for bad input

**Phase mapping:** Phase 5 (Pipeline integration and S3 contracts) - API surface testing

**DKGE-specific context:** The package has `print.dkge_contrasts`, `as.matrix.dkge_contrasts`, `as.data.frame.dkge_contrasts`, and various other S3 methods. The `as.matrix` method throws a custom condition class `dkge_transport_needed` - this error path should be tested (and is, at `test-contrast.R:281-289`).

---

### Pitfall 9: Effect Alignment Order Sensitivity

**What goes wrong:** Results depend on the order of effects in the design matrix or kernel, even though the mathematical model should be order-invariant.

**Why it happens:** Code assumes effects are in a specific order, uses numeric indices instead of names, or has implicit ordering from `match()` or `intersect()` operations.

**Consequences:**
- Results differ between users who order effects differently
- Subtle bugs when effects are added/removed
- Incorrect contrast application

**Warning signs:**
- Tests always use effects in alphabetical order
- Code uses `effects[[1]]` instead of `effects[["effect_name"]]`
- No tests with shuffled effect orders

**Prevention:**
1. **Shuffle tests:**
   ```r
   test_that("predict respects effect reordering", {
     shuffled <- lapply(betas, function(mat)
       mat[rev(seq_len(nrow(mat))), , drop = FALSE])
     preds_shuffled <- dkge_predict_loadings(fit, shuffled)
     preds <- dkge_predict_loadings(fit, betas)
     expect_equal(preds_shuffled, preds, tolerance = 1e-10)
   })
   ```
2. **Test non-alphabetical naming:** Use effect names like `c("zebra", "alpha", "middle")`
3. **Test partial overlap with different orderings**
4. **Use named vectors for contrasts, verify name matching**

**Phase mapping:** Phase 3 (Effect alignment and kernel handling) - alignment robustness

**DKGE-specific context:** The test at `test-project.R:97-102` verifies effect reordering works. But `dkge_align_effects.R` has complex order-dependent logic (`order(obs)`, `sort(obs_idx)`) that needs thorough testing.

---

### Pitfall 10: Parallel Execution Non-Determinism

**What goes wrong:** Tests pass sequentially but fail when `parallel = TRUE`. Results are numerically different or order-dependent.

**Why it happens:** Parallel execution changes:
- Random number generation (each worker has different stream)
- Floating-point accumulation order (sum is not associative)
- Result collection order

**Consequences:**
- Unreproducible results for users with different core counts
- CI flakiness
- Subtle numerical differences between runs

**Warning signs:**
- Tests skip `parallel = TRUE` code paths
- `set.seed()` doesn't ensure reproducibility
- Results differ between `future::plan("sequential")` and `future::plan("multicore")`

**Prevention:**
1. **Test parallel explicitly:**
   ```r
   test_that("parallel and sequential give same results", {
     future::plan("sequential")
     result_seq <- dkge_contrast(fit, c, parallel = TRUE)  # uses sequential plan
     future::plan("multicore")
     result_par <- dkge_contrast(fit, c, parallel = TRUE)
     expect_equal(result_seq$values, result_par$values, tolerance = 1e-10)
   })
   ```
2. **Use `future.apply` which handles RNG properly**
3. **Test with small and large inputs** (parallelization thresholds)
4. **Verify accumulation order doesn't matter:** If order matters, document or fix

**Phase mapping:** Phase 5 (Pipeline integration) - test parallel code paths

**DKGE-specific context:** The `.dkge_apply()` helper wraps parallel execution. Contrast computation loops at `dkge-contrast.R:206-219` and `dkge-analytic.R:306-324` use this. Need tests that verify `parallel = TRUE` gives same results.

---

## Minor Pitfalls

Issues that cause annoyance but are quickly fixable.

### Pitfall 11: Test Helper Pollution

**What goes wrong:** Helper functions in `helper-*.R` files define objects that leak into test namespace, causing tests to depend on execution order.

**Prevention:** Use local test fixtures created within each test, or use `testthat::local_*` for cleanup.

**Phase mapping:** All phases - maintain test hygiene

---

### Pitfall 12: Missing `skip_on_cran()` for Slow Tests

**What goes wrong:** CRAN rejects package because tests take too long, or slow tests are skipped everywhere including local development.

**Prevention:** Use `skip_on_cran()` for slow tests, `skip_on_ci()` only for truly environment-specific tests. Aim for CRAN tests under 2 minutes.

**Phase mapping:** All phases - tag slow tests

---

### Pitfall 13: Undocumented Test Dependencies

**What goes wrong:** Tests require packages not in `Suggests`, causing failures on minimal installations.

**Prevention:** Use `skip_if_not_installed("package")` or ensure all test dependencies are in `Suggests`.

**Phase mapping:** Phase 5 (Pipeline integration) - dependency management

---

## Phase-Specific Warning Summary

| Testing Phase | Key Pitfalls to Watch | Priority Checks |
|--------------|----------------------|-----------------|
| Phase 1: Core fit | Eigenvector ambiguity, tolerance calibration, wrong invariant | Recovery tests, K-orthonormality, reconstruction error |
| Phase 2: LOSO/Cross-fitting | Data leakage, analytic fallback | Basis differs, held-out exclusion verified, fallback triggers |
| Phase 3: Effect alignment | Order sensitivity, partial overlap | Shuffle tests, minimal overlap tests, prior completion |
| Phase 4: Edge cases | Seed dependence, rank deficiency | Multi-seed, zero inputs, near-singular matrices |
| Phase 5: Integration | S3 dispatch, parallel execution | All arguments, parallel/sequential equivalence |

---

## Sources

### R Package Testing
- [testthat documentation (CRAN, 2026)](https://cran.r-project.org/web/packages/testthat/testthat.pdf)
- [R Packages (2e) - Testing Basics](https://r-pkgs.org/testing-basics.html)
- [testthat 3.3.0 release notes](https://tidyverse.org/blog/2025/11/testthat-3-3-0/)

### Floating-Point Comparison
- [fpCompare package](https://cran.r-project.org/web/packages/fpCompare/index.html)
- [Comparing Floating Point Numbers in R](https://gcdi.commons.gc.cuny.edu/2023/03/15/comparing-floating-point-numbers-in-r/)
- [Numerical pitfalls in computing variance](https://www.r-bloggers.com/2016/05/numerical-pitfalls-in-computing-variance/)

### Numerical Algorithm Testing
- [Drake: Accuracy, Tolerance, and Precision](https://drake.mit.edu/doxygen_cxx/group__accuracy__and__tolerance.html)
- [Testing Math Functions in Microsoft Cloud Numerics](https://learn.microsoft.com/en-us/archive/msdn-magazine/2012/october/numerics-testing-math-functions-in-microsoft-cloud-numerics)
- [Numerical Stability in Eigenvalue Decomposition](https://www.numberanalytics.com/blog/numerical-stability-eigenvalue-decomposition-techniques-best-practices)

### Cross-Validation and Data Leakage
- [scikit-learn: Common Pitfalls](https://scikit-learn.org/stable/common_pitfalls.html)
- [Being Aware of Data Leakage and Cross-Validation Scaling](https://analyticalsciencejournals.onlinelibrary.wiley.com/doi/full/10.1002/cem.70026)
- [LightGBM data leakage issue](https://github.com/microsoft/LightGBM/issues/4319)

### Verification and Validation
- [NAS: Verification in Mathematical Models](https://nap.nationalacademies.org/read/13395/chapter/5)
- [Archive of Numerical Software: Testing Numerical Code](https://journals.ub.uni-heidelberg.de/index.php/ans/article/download/27447/29545)
