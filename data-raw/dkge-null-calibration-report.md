# DKGE null-calibration audit: version 1

Executed 2026-08-16 under the predeclared plan in
`data-raw/dkge-null-calibration-plan.md`. No thresholds, seeds, estimands, or
simulation settings were changed after observing results.

## Decision

- Adaptive LOSO sign-flip inference is calibrated in this experiment. The
  earlier `kenergy` and `kenergy_prec` failures are classified as Monte Carlo
  variation: all three 1,000-replicate arms pass the original KS and quantile
  gates, with empirical size 0.040 (`none`), 0.051 (`kenergy`), and 0.052
  (`kenergy_prec`).
- Between-subject Freedman-Lane inference is materially anti-conservative in
  the primary small-sample arm. With 20 subjects, empirical size is 0.074 for
  `group`, 0.076 for `trait`, and 0.077 for `group:trait`; every 95% Wilson
  lower bound exceeds 0.05 and every original uniformity gate fails.
- The between-subject distortion attenuates with sample size but is not
  monotone term by term. At 40 subjects the sizes are 0.050, 0.066, and 0.070;
  at 80 they are 0.030, 0.048, and 0.052. Together with the independent
  fast-versus-explicit-refit oracle already in the test suite, this supports a
  finite-sample method limitation rather than a discrepancy in the optimized
  implementation.

The package therefore retains the estimator but now states the observed
small-sample limitation in `dkge_between_permute()` documentation. Increasing
the permutation count improves p-value resolution; it does not correct the
residual-exchangeability distortion. Small-sample p-values close to alpha must
not be represented as fully calibrated confirmatory evidence.

## Frozen evidence

| Family and setting | Estimand | Replicates | Empirical size | Wilson 95% interval | KS p | Max quantile deviation | Original gate | Classification |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| Adaptive, n=16 | none | 1,000 | 0.040 | 0.0295–0.0540 | 0.935 | 0.0120 | pass | calibrated |
| Adaptive, n=16 | kenergy | 1,000 | 0.051 | 0.0390–0.0664 | 0.902 | 0.0095 | pass | calibrated |
| Adaptive, n=16 | kenergy_prec | 1,000 | 0.052 | 0.0399–0.0676 | 0.111 | 0.0260 | pass | calibrated |
| Between, n=20 | group | 1,000 | 0.074 | 0.0594–0.0919 | 3.76e-7 | 0.0775 | fail | anti-conservative |
| Between, n=20 | trait | 1,000 | 0.076 | 0.0611–0.0941 | 2.60e-5 | 0.0713 | fail | anti-conservative |
| Between, n=20 | group:trait | 1,000 | 0.077 | 0.0620–0.0952 | 0.00189 | 0.0500 | fail | anti-conservative |

The complete 12-row diagnostic table is
`inst/extdata/dkge-null-calibration-summary.csv`; all 9,000 underlying p-values
are in `inst/extdata/dkge-null-calibration-replicates.csv`.

## Provenance

- Calibration identifier: `null-calibration-v1`
- Frozen plan SHA-256:
  `d671401031967ce8260ee4399c2cfe592af1ca4b0dd17289fbdd013fc31c7a74`
- Runner SHA-256:
  `394f52fc1e49bfa5627c4eef8596ac217756defb89d757a505313ea0c60ad886`
- Replicate artifact SHA-256:
  `276abe7a057c1702c9930d7d4e3e05b2fd22a989db099b75163fe2f3e1b7c0b3`
- Summary artifact SHA-256:
  `1706d05375f026cce15884acc1d761e6aee1380a8911e9fbcc52582df043cf2e`
- Metadata artifact SHA-256:
  `4d4b9a5404c835603f3a9486e232a9e324222dd230f27cda7f8036d44527bad6`
- Recorded source HEAD: `953d2a64ff356d55f4a6fd3278a7f11db018fea6`
- Runtime: R 4.5.1 on `aarch64-apple-darwin20`; 1,445.0 seconds.
