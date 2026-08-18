# DKGE Tier 0-2 validation report

Date: 2026-08-15  
Branch: `codex/dkge-trialwise-trust`  
Base: `ca48d70`  
Mote epic: `bd-01M03KRRJ77XHPH0DDQ9V0YK6F`

## Release disposition

The Tier 0-2 implementation is numerically and structurally much stronger than
the historical path, but the package is **not yet release-ready**. The new
algebraic, leakage, simulation, performance, and API gates pass as described
below. Package-wide release gates remain open because of pre-existing test and
vignette failures, pre-existing long-null calibration misses, one performance
dispersion warning, one unavailable optional dependency, and the absence of an
independent reviewer for this final diff.

No commit, push, hosted CI run, or release was performed. Local evidence must
not be interpreted as publication evidence.

## Implemented estimator contract

- A named global effect grid and deterministic union alignment preserve local
  empty-cell provenance without treating zero padding as observation.
- Subject-by-effect reliability supports counts, explicit precision, and
  capped DerSimonian-Laird random-effects precision.
- Dense and feature-chunked trialwise constructors support IID OLS and
  covariance/precision-aware GLS, effect covariance, residual scale, analytic
  noise moments, and independently audited or exploratory split halves.
- The canonical pooled estimator remains a symmetric eigendecomposition of a
  `q x q` moment. Negative directions introduced by correction are diagnosed;
  only the positive numerical spectrum is retained. No `P x P` eigensolver is
  needed or formed.
- Reliability-normalized, debiased, ridge, CPCA, and JD fits advertise a
  q-space-only representation. Physical block projection APIs fail closed;
  supported q-space projection, cell classification, contrasts, folds, and
  resampling use the explicit representation boundary.
- LOSO, k-fold, and literal Poisson resamples refit the ruler, reliability,
  random-effects variance, debiasing, and missingness using the training
  subject multiset.
- Fit diagnostics name the active moment estimator and report coverage, pair
  effective sample size, negative spectral mass, error-model provenance,
  solver, and representation.

## Numerical and computational evidence

Independent dense oracles cover raw and corrected subject moments, pairwise
precision normalization, covariance-aware analytic subtraction, the pooled
ruler, the kernel transform, the positive-spectrum eigenspace, projections,
and permutation/metamorphic invariants. Literal subject duplication agrees with
Poisson multiplicity refits. Leakage-injection tests show held-out reliability
and count perturbations cannot affect training folds.

The frozen 3 by 5 by 4 simulation contains 10 scenarios, 20 fixed seeds, 18
subjects, 60 effects, and seven estimators (1,400 replicate rows). Artifacts
were reproduced byte-for-byte in two clean executions:

- replicates SHA-256:
  `1b43362d26822007ed288c767318a4f364434eabb934a0bbad7f46ca6fc36030`
- summary SHA-256:
  `b5113617038b5d778be4f7d5410fe8676e71633b99299b9798758e54111c9cda`

Ten of 14 predeclared scientific gates pass. Important results:

- covariance-aware analytic correction versus legacy under hard noise:
  subspace-angle ratio 0.725 and contrast-RMSE ratio 0.184;
- covariance-aware analytic and split-half inverse-count-correlation
  improvements: -0.188 and -0.162;
- null false-positive rates: 0.05 for both corrected estimators;
- IID equivalence maximum angle difference: `1.07e-14` degrees;
- failures retained: explicit precision and random-effects weighting alone do
  not reduce the inverse-count artifact; random-effects high-count influence
  ratio is 1.0 under this zero-spatial-mean data-generating process;
  covariance-analytic interval coverage is 0.694 against the 0.80 floor.

The frozen q=60 benchmark ran five workloads in three fresh processes each.
It records elapsed time, peak RSS, total and largest allocations, returned
object size, and retained design size separately. No workload detected a
`P x P` allocation. At `T=240`, `q=60`, `P=100,000`, the chunked trial path
returns a 75.4 MB object, retains a 119 KB design, and has a 48.0 MB largest
allocation. Compact literal repools retain about 0.31 MB. Four of five gates
pass. The preserved miss is scheduling dispersion in the Poisson workload:
median 5.10 s, IQR 6.62 s; its median peak RSS (1.23 GB), elapsed ceiling, and
allocation gates pass. Artifact SHA-256:
`a8832319fe65f2c9d9ad15e7fd81afa04579a6cf976112ce7bcf8a26d2af6cd7`.

## Package and test gates

| Gate | Outcome |
|---|---|
| Focused Tier 0-2 tests | Pass |
| Full `devtools::test()` | Completes with two failures independently reproduced on clean `ca48d70`; no Tier 0-2 failure |
| Installed-tarball tests | 1,839 pass, 2 fail, 51 warnings, 50 skip under R CMD check's CRAN mode |
| R code/static checks | Syntax, load/unload, dependencies, S3 consistency, possible problems, Rd, examples, and compiled-code checks pass |
| New 3x5x4 vignette | Renders successfully standalone and during the full build |
| Exact full vignette build | Blocked by pre-existing `dkge-concepts.Rmd` call to absent `dkge_aggregate_fit()` |
| Optional T4transport path | Not exercised; dependency unavailable |

The final audit added the missing `tests/testthat.R` entry point. Before this
repair, `R CMD check` installed the tarball but silently ran none of the 80+
test files. The tarball check now exposes the same two clean-base failures:

1. `test-between-rrr.R:336`: a stochastic Sinkhorn non-convergence warning is
   asserted away;
2. `test-coverage-gaps.R:241`: stale expected names omit the already returned
   `p_unadj` field.

R CMD check skips 48 tests marked `skip_on_cran()`, the unavailable
T4transport test, and the inverse optional-dependency test. Those 48 tests were
exercised by the local `devtools::test()` run except the two explicit long-null
tests, which were run separately.

## Long-null calibration

The opt-in assertion harness used unsupported `info=` arguments with
`expect_gt()` and `expect_lt()`; that test bug was repaired without changing a
statistical threshold. The long tests then revealed deterministic,
independently reproduced clean-base deficits:

| Test | KS p | Maximum quantile deviation | Gate |
|---|---:|---:|---|
| between group | 0.0694 | 0.0875 | pass |
| between trait | 0.0447 | 0.0844 | fail KS > 0.05 |
| between group:trait | 0.9419 | 0.0378 | pass |
| adaptive none | 0.7456 | 0.0490 | pass |
| adaptive kenergy | 0.0933 | 0.1025 | fail deviation < 0.08 |
| adaptive kenergy_prec | 0.0735 | 0.0810 | fail deviation < 0.08 |

Thresholds were not changed and failed outcomes were not rerun away.

## Documentation and claim boundary

The executable vignette and migration note explicitly state that the shared
`T x q` design is retained, that `scope` is metadata rather than a kernel
modifier, and that three-factor kernel terms must be enumerated when two-way
interactions are intended. They distinguish weighting from noise-bias removal,
state the separable trial/spatial covariance assumption, explain when
trialwise first-level betas violate IID errors, and label deterministic
within-cell splits exploratory unless independence is justified.

The evidence supports covariance-aware moment correction and independent split
cross-moments as the principal Tier 2 protections against unequal-count noise
bias. It does not support claims that count or random-effects weighting alone
removes that bias, that covariance-analytic intervals are calibrated in the
frozen scenarios, or that all package-level release gates are green.

## Remaining release blockers

1. Repair or remove the stale aggregate-fit vignette call, then complete a full
   vignette-bearing tarball build and ordinary R CMD check.
2. Resolve or explicitly re-baseline the two clean-base unit-test failures.
3. Investigate the trait, `kenergy`, and `kenergy_prec` long-null calibration
   misses; do not relax thresholds without a new predeclared plan.
4. Exercise the T4transport path in an environment that supplies it, or accept
   the optional-path gap explicitly.
5. Obtain a genuinely independent review of model assumptions, q-space
   algebra, leakage boundaries, fail-closed APIs, and claim wording.
6. Repeat the Poisson performance workload on a quiet reference host if a
   green dispersion gate is required.
