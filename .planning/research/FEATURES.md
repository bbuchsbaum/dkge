# Feature Landscape: Testing for Publication-Ready Scientific R Packages

**Domain:** Scientific R package testing (statistical/mathematical algorithms)
**Researched:** 2026-01-19
**Package:** dkge - Design-Kernel Group Embedding for fMRI analysis

## Table Stakes

Features users and reviewers expect. Missing = package lacks credibility for publication.

### Unit Testing Framework

| Feature | Why Expected | Complexity | Notes |
|---------|--------------|------------|-------|
| testthat 3rd edition | Industry standard for R packages; required by CRAN implicitly | Low | Already configured in DESCRIPTION |
| Test organization by function/module | Reviewers expect logical test structure | Low | Group tests by layer: data, kernel, fit, contrast, etc. |
| Deterministic test seeds | Reproducibility of stochastic tests | Low | Use `set.seed()` at test file or block level |
| Tolerance-based numerical comparisons | Floating-point arithmetic varies by platform | Low | Use `expect_equal(x, y, tolerance = 1e-8)` not `expect_identical()` |

### Mathematical Invariant Tests

| Feature | Why Expected | Complexity | Notes |
|---------|--------------|------------|-------|
| Orthonormality verification | U^T K U = I is core contract per algo.md | Low | Already in test-fit.R |
| Matrix reconstruction tests | K^{1/2} K^{1/2} = K, K^{1/2} K^{-1/2} = I | Low | Already in test-fit.R |
| Symmetry verification | Chat, K must be symmetric | Low | `expect_equal(M, t(M))` |
| Positive semi-definiteness | K, Chat eigenvalues >= 0 | Low | Check `all(eigen(M)$values >= -tol)` |
| Marginal constraints (OT) | Transport plan row/col sums match marginals | Medium | Already in test-transport-sinkhorn.R |

### API Contract Tests

| Feature | Why Expected | Complexity | Notes |
|---------|--------------|------------|-------|
| Input validation error messages | Functions should fail gracefully with clear errors | Low | Use `expect_error(..., regexp = "meaningful message")` |
| Dimension mismatch detection | Users make mistakes; early failure saves debugging | Low | Already in test-fit.R |
| Return value structure | S3 classes return expected components | Low | Check class, names, dimensions |
| NULL/missing handling | Edge cases must not crash | Medium | Test with NULL Omega, missing weights, etc. |

### Cross-Platform Compatibility

| Feature | Why Expected | Complexity | Notes |
|---------|--------------|------------|-------|
| R CMD check passes | CRAN requirement; no NOTEs, WARNINGs, ERRORs | Medium | Run with `--as-cran` flag |
| Multi-platform CI | CRAN tests on Linux, macOS, Windows | Medium | Use r-lib/actions GitHub workflows |
| R version matrix | CRAN tests on release, devel, oldrel | Medium | Test R 4.2+, current, devel |

### Documentation Tests

| Feature | Why Expected | Complexity | Notes |
|---------|--------------|------------|-------|
| All exported functions have @examples | CRAN requirement; examples are run during check | Medium | Examples must complete in < 5 sec |
| Examples actually execute | `devtools::run_examples()` must pass | Low | Mark slow examples with `\donttest{}` |
| Vignette code runs | Vignettes are executed during build | High | Vignettes often most time-consuming |

### Code Coverage

| Feature | Why Expected | Complexity | Notes |
|---------|--------------|------------|-------|
| Coverage measurement with covr | Standard practice for quality packages | Low | `covr::package_coverage()` |
| Coverage badge in README | Signals quality to users | Low | Integrate with Codecov or Coveralls |
| Exported functions covered | All public API should have tests | Medium | Target 80%+ line coverage |

## Differentiators

Features that set package apart. Not expected, but valued for publication credibility.

### Reference Implementation Comparison

| Feature | Value Proposition | Complexity | Notes |
|---------|-------------------|------------|-------|
| Cross-check against external packages | Validates correctness against known implementations | Medium | Already done: T4transport for Sinkhorn |
| Cross-language validation | Python/MATLAB comparison builds trust | High | Consider NumPy/SciPy for eigen, OT |
| Published algorithm reproduction | Match figures/tables from papers | High | Reproduce algo.md equations exactly |

### Property-Based Tests

| Feature | Value Proposition | Complexity | Notes |
|---------|-------------------|------------|-------|
| Theoretical bound verification | Results within mathematical limits | Medium | E.g., contrast values bounded by input range |
| Monotonicity tests | Rank k+1 explains >= rank k variance | Medium | Eigenvalue ordering, cumulative variance |
| Sensitivity analysis | Single parameter changes produce expected direction | Medium | More subjects -> more stable estimates |

### Recovery Tests (Ground Truth)

| Feature | Value Proposition | Complexity | Notes |
|---------|-------------------|------------|-------|
| Simulated data recovery | Recover known signal from synthetic data | Medium | Already done: dkge_sim_toy tests |
| SNR relationship | Higher SNR -> better recovery | Medium | Already in test-toy-scenarios.R |
| Component recovery | Cosine similarity to true U > threshold | Medium | Already in test-toy-recovery.R |

### Snapshot/Golden Tests

| Feature | Value Proposition | Complexity | Notes |
|---------|-------------------|------------|-------|
| Numeric result snapshots | Detect unintended changes in algorithms | Medium | Use `expect_snapshot_value()` with tolerance |
| Output format snapshots | Detect changes in print/summary methods | Low | Use `expect_snapshot()` for console output |
| Reference data files | Permanent test fixtures with known answers | Medium | Store in tests/testthat/fixtures/ |

### Statistical Inference Tests

| Feature | Value Proposition | Complexity | Notes |
|---------|-------------------|------------|-------|
| Type I error control | Sign-flip p-values uniform under null | High | Simulate null, check p-value distribution |
| FWER control | Max-T correction controls family-wise error | High | Multiple comparison simulation |
| Power analysis | Method detects effects at expected rate | High | Requires simulation study |

### Performance Benchmarks

| Feature | Value Proposition | Complexity | Notes |
|---------|-------------------|------------|-------|
| Scaling tests | O(q^3) not O(P^3) as claimed | Medium | Time with varying q, P |
| Memory profiling | No large P x P matrices formed | Medium | Use bench::mark() or profvis |
| Streaming equivalence | Streaming matches batch within FP tolerance | Medium | Already implied in algo.md |

### Numerical Stability Tests

| Feature | Value Proposition | Complexity | Notes |
|---------|-------------------|------------|-------|
| Near-singular K handling | Ridge prevents divide-by-zero | Medium | Test with ill-conditioned kernels |
| Small eigenvalue gaps | Basis stable despite near-degenerate eigenvalues | High | Perturb Chat, check U stability |
| Double precision fidelity | Results don't change with redundant operations | Low | Symmetrize, check consistency |

### LOSO Cross-Validation Tests

| Feature | Value Proposition | Complexity | Notes |
|---------|-------------------|------------|-------|
| Unbiased estimation | LOSO values differ from full-sample | Medium | Statistical difference test |
| Leave-out mechanics | Chat^(-s) = Chat - w_s * S_s exactly | Low | Already testable from fit object |
| Analytic approximation accuracy | Approximation close to exact when perturbation small | Medium | Already in test-analytic.R |

## Anti-Features

Testing approaches to explicitly AVOID for scientific code.

| Anti-Feature | Why Avoid | What to Do Instead |
|--------------|-----------|-------------------|
| Testing only happy path | Misses edge cases that cause real failures | Test boundaries: n=1, n=0, singular, NA |
| Exact floating-point equality | Fails randomly across platforms | Always use tolerance for numerics |
| Overly loose tolerances | Masks real bugs | Use tightest tolerance that passes reliably |
| Testing implementation details | Tests break on refactor | Test API contracts and mathematical properties |
| Stochastic tests without seeds | Non-reproducible failures | Set seed at beginning of each test block |
| Mocking mathematical operations | Defeats purpose of correctness tests | Mock I/O and external services, not math |
| Testing only with ideal data | Real data has issues | Include tests with missing data, outliers |
| Snapshot tests for numerics | Fragile across platforms, R versions | Use tolerance-based comparison instead |
| Long-running tests in main suite | Slows development cycle | Use `skip_on_cran()`, separate slow test suite |
| Testing without reference | No way to know if answer is correct | Compare to published results or alternative implementations |

## Feature Dependencies

```
Code Coverage
    |
    v
Unit Tests (testthat) --> API Contract Tests --> Input Validation Tests
    |                          |
    v                          v
Mathematical Invariant Tests   Error Handling Tests
    |
    v
Property-Based Tests --> Recovery Tests (need invariants first)
    |
    v
Reference Implementation Comparison (needs working property tests)
    |
    v
Statistical Inference Tests (needs all above for credibility)
    |
    v
Performance Benchmarks (only meaningful after correctness established)
```

## MVP Recommendation

For publication-ready minimum:

1. **Mathematical Invariants** (HIGH PRIORITY)
   - All algo.md equations verified: U^T K U = I, Chat reconstruction, LOSO mechanics
   - These are the core scientific claims

2. **API Contracts** (HIGH PRIORITY)
   - Input validation with clear errors
   - Return value structure verification
   - All exported functions have working examples

3. **Recovery Tests** (HIGH PRIORITY)
   - Synthetic data with known ground truth
   - Demonstrate method works as intended

4. **Cross-Platform CI** (MEDIUM PRIORITY)
   - GitHub Actions for Linux, macOS, Windows
   - R release, devel, oldrel-1

5. **Code Coverage > 80%** (MEDIUM PRIORITY)
   - All exported functions covered
   - Critical paths (fit, contrast, LOSO) > 90%

**Defer to post-publication:**
- Extensive performance benchmarks
- Full statistical power analysis
- Cross-language validation
- Interactive test fixtures

## Complexity Summary

| Category | Table Stakes | Differentiators | Total Effort |
|----------|-------------|-----------------|--------------|
| Unit Testing Framework | Low | - | Low |
| Mathematical Invariants | Low-Medium | Medium-High | Medium |
| API Contracts | Low-Medium | - | Low-Medium |
| Cross-Platform CI | Medium | - | Medium |
| Documentation Tests | Medium | - | Medium |
| Code Coverage | Low-Medium | - | Low-Medium |
| Reference Comparison | - | Medium-High | Medium-High |
| Property-Based Tests | - | Medium | Medium |
| Recovery Tests | - | Medium (mostly done) | Low (extend) |
| Statistical Inference | - | High | High |
| Performance Benchmarks | - | Medium | Medium |

**Estimated total effort for publication-ready:** Medium-High (most table stakes exist; need systematic completion)

## Sources

### Authoritative References
- [testthat R Package Documentation](https://cran.r-project.org/web/packages/testthat/testthat.pdf) - Unit testing for R, version 3.3.2 (Jan 2026)
- [R Packages (2e) - Testing Basics](https://r-pkgs.org/testing-basics.html) - Hadley Wickham's canonical guide
- [covr - Test Coverage for Packages](https://cran.r-project.org/web/packages/covr/index.html) - Code coverage measurement

### Standards and Guidelines
- [CRAN Submission Checklist](https://cran.r-project.org/web/packages/submission_checklist.html) - Official CRAN requirements
- [rOpenSci Packaging Guide](https://devguide.ropensci.org/pkg_building.html) - Best practices for scientific R packages
- [Bioconductor Unit Tests](https://contributions.bioconductor.org/tests.html) - Bioconductor testing standards
- [rOpenSci CI Best Practices](https://devguide.ropensci.org/pkg_ci.html) - Continuous integration guidelines

### Testing Statistical Correctness
- [Epiverse-TRACE: Statistical Correctness Testing](https://epiverse-trace.github.io/posts/statistical-correctness/) - Reference implementation comparison, property testing
- [testthat Snapshot Testing](https://testthat.r-lib.org/articles/snapshotting.html) - Golden tests documentation
- [testthat Equality Expectations](https://testthat.r-lib.org/reference/equality-expectations.html) - Numerical tolerance handling

### Package Quality Assessment
- [R Package Quality: Code Quality](https://www.r-bloggers.com/2025/07/r-package-quality-code-quality/) - Quality indicators (July 2025)
- [goodpractice Package](https://cran.r-project.org/web/packages/goodpractice/vignettes/goodpractice.html) - Automated quality checks
- [R Package Validation in Pharma](https://www.r-bloggers.com/2024/10/a-guide-to-r-package-validation-in-pharma/) - Validation requirements

### Reproducibility
- [Ten Simple Rules for R Packages](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009884) - PLOS Computational Biology guidance
- [Primer on Reproducible Research in R](https://pmc.ncbi.nlm.nih.gov/articles/PMC10969410/) - PMC reproducibility guide
