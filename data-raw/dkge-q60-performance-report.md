# DKGE q=60 performance and dispersion report

Executed 2026-08-16 with the unchanged runner, seeds, workloads, and thresholds
predeclared in `data-raw/dkge-q60-performance-plan.md`. All 15 repetitions are
retained in `inst/extdata/dkge-q60-performance.csv`.

## Decision

All five workload gates pass. The earlier Poisson-bootstrap outlier does not
reproduce: the new three-run range is 1.313–1.376 seconds, median 1.364 seconds,
and IQR 0.0315 seconds, compared with the earlier range 4.406–17.651 seconds.
Because the exact workload and seeds are unchanged and the new run was made on
a host with substantial unrelated CPU load, the result strongly supports a
host/session-variance explanation rather than a reproducible DKGE runtime tail.
It does not satisfy the tracker's separate requirement for a quiet reference-
host repetition, so that ticket remains open for the literal final acceptance
run.

The first sandboxed run retained valid timing and allocation evidence but could
not collect RSS because macOS `/usr/bin/time -l` was denied its `sysctl`
query. It was not used for the final gate. The full unchanged runner was then
repeated with OS counter access; that second, complete 15-run artifact is the
versioned result below.

## Results

| Workload | Median elapsed (s) | IQR (s) | Median peak RSS (GB) | Largest allocation (MB) | Returned object (MB) | Pass |
|---|---:|---:|---:|---:|---:|---:|
| trial_chunks | 7.439 | 0.3815 | 1.013 | 48.0 | 75.4 | yes |
| group_legacy | 3.307 | 0.8980 | 1.210 | 86.4 | 232.4 | yes |
| group_advanced | 3.145 | 0.1440 | 1.256 | 86.4 | 233.1 | yes |
| fold_refit | 3.306 | 0.1300 | 1.259 | 86.4 | 0.392 | yes |
| poisson_bootstrap | 1.364 | 0.0315 | 0.747 | 17.3 | 0.185 | yes |

No run allocated a P-by-P object. At P=100,000 the chunked trialwise path kept
the largest single allocation to 48.0 MB, versus 80 GB for a dense P-by-P
double matrix. Its retained design is 119,480 bytes and its returned object is
75.4 MB, both below the frozen limits.

## Reconciliation found by the gate

The first benchmark attempt failed before measurement because this checkout
did not expose `dkge_trial_subject_chunks()`, although the predeclared Tier-2
imaging workload and an isolated development branch relied on it. The public
chunked OLS constructor was restored with the same estimator and split-half
semantics as the dense constructor present in this checkout. New tests compare
coefficients, effect covariance, residual variance, noise trace, and split
halves against the dense oracle, and test function-source termination and
malformed input. The benchmark itself was not weakened or rewritten.

This restoration is deliberately limited to the IID OLS model implemented by
the current dense constructor. Covariance-aware GLS chunking from the isolated
development branch has not been silently imported and is not claimed here.

## Provenance

- Plan SHA-256:
  `b9d5892b1d4efdd8096afcbf2bcce3a7f04628bf0718f56722a9571e9d73f5b6`
- Runner SHA-256:
  `6294565d593da8d9857e66a4fb647e5805575118f6da10be0c4df1605170c8e6`
- Result SHA-256:
  `16393f583966edccce10bea77dde4970be1970e0e920961663fdbeed27718b7c`
