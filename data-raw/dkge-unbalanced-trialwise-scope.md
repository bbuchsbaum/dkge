# DKGE unbalanced trialwise hardening scope

This document is the path and integration manifest for Mote epic
`bd-01M03KRRJ77XHPH0DDQ9V0YK6F`.

## Integration baseline

- Isolated worktree: `/tmp/dkge-trust-hardening`
- Branch: `codex/dkge-trialwise-trust`
- Base: local `origin/main` at `ca48d70`
- Source checkout: `/Users/bbuchsbaum/code/dkge`, preserved unchanged
- Historical implementation epic: `bd-01KR6VV5ZG9ZRHABD832ZFDZ8K`

The isolated branch includes the upstream correctness fixes already present on
the local `origin/main` reference: deterministic subject weighting,
scale-relative eigenvalue tolerance, non-finite beta validation, safe empty
block indexing, and fold-local reliability subsetting.

## Included paths and owners

| Paths | Owning ticket | Purpose |
|---|---|---|
| `R/dkge-align-data.R`, `R/dkge-data.R` | `bd-01M03KTTSNK674XG3FWBV9B5H6` | Observed-row and union-alignment contract |
| `R/dkge-fit.R`, `R/dkge-fit-core.R`, `R/dkge-moments.R`, `R/dkge-folds.R`, `tests/testthat/test-partial-coverage-algebra.R` | `bd-01M03KTTSNK674XG3FWBV9B5H6` | Canonical raw-mask, right-weight, ruler, kernel, and fold convention |
| `R/dkge-effect-grid.R`, `R/design-kernel.R` | `bd-01M03KTTSNK674XG3FWBV9B5H6` | Named global effect grid and kernel metadata |
| `R/dkge-effect-weights.R`, `R/dkge-moments.R`, `R/dkge-fit-core.R`, `tests/testthat/test-random-effects-reliability.R` | `bd-01M03KTW6MEB3SBKDH8DJA7PFZ` | Count, explicit, and DL random-effects reliability; structural-zero equivalence; ESS and dominance diagnostics |
| `R/dkge-trial-subject.R`, `R/dkge-data.R`, `R/dkge-align-data.R`, `R/dkge-moments.R`, `R/dkge-fit.R`, `R/dkge-fit-core.R`, `tests/testthat/test-covariance-aware-trialwise.R` | `bd-01M03KTVCXRHQ19XER1YGT7WVB` | IID/GLS trialwise sufficient statistics, covariance-aware ruler, analytic noise provenance, and AR(1) bias oracle |
| `R/dkge-trial-subject.R`, `R/dkge-align-data.R`, `R/dkge-data.R`, `tests/testthat/test-split-half-reliability.R` | `bd-01M03KTVWBFDW0TAVR7EY87680` | Run-aware and explicit split provenance, independence audit, covariance-consistent half fitting, and bounded split reliability export |
| `R/dkge-fit.R`, `R/dkge-fit-core.R`, `R/dkge-project.R`, `data-raw/dkge-qspace-algebra-contract.md`, `tests/testthat/test-fit-algebra-contract.R` | `bd-01M03KTTKQHAE2BYW6C7RBQY2E` | Fit algebra, eigensolve, representation boundary, and fail-closed projection |
| `R/dkge-predict.R`, `R/dkge-project.R`, `R/dkge-latent-utils.R`, `R/dkge-loso.R`, `R/dkge-analytic.R`, `R/dkge-bootstrap.R`, `R/dkge-plot.R`, `tests/testthat/test-advanced-fit-representation.R` | `bd-01M03KTV4QZQHHJBGT6C8P63KS` | Coherent q-space projection semantics and downstream support/fail-closed boundary |
| `R/dkge-folds.R`, `R/dkge-kfold.R`, `R/dkge-loso.R`, `R/dkge-contrast.R`, `R/dkge-analytic.R`, `R/dkge-bootstrap.R`, `R/dkge-moments.R`, `tests/testthat/test-resampling-advanced-moments.R` | `bd-01M03KTX246WMJTV4WHSCGVZDD` | Leakage-free fold refits, literal Poisson multiplicities, analytic fallback, estimator provenance, and null calibration |
| `tests/testthat/test-data.R`, `tests/testthat/test-fit.R`, `tests/testthat/test-effect-reliability.R`, `tests/testthat/test-trial-subject.R`, `tests/testthat/test-design-kernel.R`, `tests/testthat/test-independent-moment-oracles.R` | `bd-01M03KTWQAS4D2EA9WK8CKJ3TW` | Regression tests plus independent dense moment, debiasing, eigenspace, projection, and metamorphic oracles |
| `data-raw/dkge-unbalanced-trialwise-simulation-plan.md`, `dev/simulate-dkge-unbalanced-trialwise.R`, `inst/extdata/dkge-unbalanced-trialwise-*.csv`, `tests/testthat/test-unbalanced-3x5x4-simulation.R` | `bd-01M03KTXAZ6PR3GS62PKGPCG2V` | Frozen 20-seed 3 x 5 x 4 scientific validation, machine-readable evidence, gates, and golden regression |
| `R/dkge-trial-subject.R`, `R/dkge-moments.R`, `R/dkge-bootstrap.R`, `data-raw/dkge-q60-performance-plan.md`, `dev/benchmark-dkge-q60.R`, `inst/extdata/dkge-q60-performance.csv`, `tests/testthat/test-trial-subject-chunks.R`, `tests/testthat/test-q60-performance-contract.R` | `bd-01M03KTXQK53KJD0HF1K8Y39CM` | Feature-chunked trial reduction, compact literal repools, frozen q = 60 performance evidence, and no-P-by-P regression gates |
| `R/dkge-diagnostics.R`, `vignettes/dkge-partial-effect-spaces.Rmd`, `data-raw/dkge-unbalanced-trialwise-migration.md`, and generated documentation | `bd-01M03KTY3XBM2MD0MVA7VTC3T4` | Public diagnostics, API contract, assumptions, examples, and migration guidance |
| `tests/testthat.R`, `tests/testthat/test-between-permute-null.R`, `data-raw/dkge-tier0-2-validation-report.md` | `bd-01M03KTYFG9YCN003FBSVPK7VF` | Tarball test entry point, repaired long-null assertion harness, and classified final validation record |

`NAMESPACE` is shared generated metadata. Only exports belonging to the paths
above are in scope.

## Explicitly deferred or excluded

The following dirty-checkout work is not part of this branch unless a later
ticket separately justifies it:

- aggregate-target APIs and tests;
- new plotting suites and SDAM report artifacts;
- unrelated between-subject RRR and inference changes;
- contrast-scope routing changes;
- test data and rendered PDFs;
- local `.mote`, `.omc`, and `.claude` state.

The partial-effect vignette is deferred until the documentation ticket so that
the integration baseline can be tested without importing unrelated contrast or
reporting features.

## Required compatibility evidence

Before this integration ticket closes:

1. focused data, fit, reliability, trialwise, and design-kernel tests pass;
2. full-suite results are recorded, with any baseline failures independently
   reproduced on the clean integration base;
3. a full vignette build succeeds from the isolated worktree;
4. `git diff --check` is clean;
5. the dirty source checkout has no new code or generated-file changes.
