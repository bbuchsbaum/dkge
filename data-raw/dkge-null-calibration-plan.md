# Predeclared DKGE null-calibration audit

Status: frozen before execution for Mote ticket
`bd-01M0404VC90A4Q9T1ZQV4G7E5B`.

## Purpose

The existing single-batch long tests gave borderline or failed results for the
between-subject `trait` term and the adaptive `kenergy` and `kenergy_prec`
procedures. This audit distinguishes Monte Carlo variation from a persistent
method limitation without changing the existing test thresholds after seeing
new results.

The package's fast closed-form Freedman-Lane implementation already has a
differential test against the explicit refit algorithm. This audit therefore
tests calibration of the stated algorithm, not equivalence of its two
implementations.

## Frozen command and artifacts

- Command: `Rscript dev/calibrate-dkge-null.R`.
- Replicate p-values:
  `inst/extdata/dkge-null-calibration-replicates.csv`.
- Diagnostics:
  `inst/extdata/dkge-null-calibration-summary.csv`.
- Runtime and provenance:
  `inst/extdata/dkge-null-calibration-metadata.csv`.
- All p-values use the package's existing finite-permutation `+1` correction.

## Primary arms

Each arm uses five independent deterministic batches of 200 null data sets,
for 1,000 p-values per estimand.

1. Between-subject Freedman-Lane: `n = 20`, `p = 24`, requested RRR rank 2,
   `B = 399`, formula `~ group * trait + age`, and terms `group`, `trait`, and
   `group:trait`. Batch base seeds are 800, 100800, 200800, 300800, and 400800.
2. LOSO adaptive sign-flip statistic: `nsub = 16`, `V = 128`, `B = 499`, for
   `none`, `kenergy`, and `kenergy_prec`. The three methods use disjoint seed
   families beginning at 101, 20101, and 40101 and separated by 100000 between
   batches.

The simulation helpers and estimators are unchanged from the existing opt-in
long tests.

## Diagnostic arms

To test whether any between-subject distortion is a small-sample residual-
exchangeability limitation, repeat the same between-subject experiment at
`n = 40` and `n = 80`. Each diagnostic arm uses five batches of 100 data sets,
`B = 399`, and base seeds beginning at 600800 and 1100800 respectively.

Evidence of attenuation toward nominal size as `n` grows, together with the
existing fast-versus-explicit-refit oracle, is classified as a method-level
small-sample limitation rather than an implementation discrepancy. Persistence
or worsening with `n` is treated as possible implementation bias.

## Frozen diagnostics

For every term or weighting method, report:

- empirical size at alpha 0.05;
- Monte Carlo standard error and a 95% Wilson interval;
- batch-specific empirical sizes and their range;
- mean and median p-value;
- Kolmogorov-Smirnov D and p-value against Uniform(0,1);
- empirical p-value quantiles at 0.10, 0.25, 0.50, 0.75, and 0.90;
- maximum absolute deviation of those quantiles from their targets.

The Wilson interval is descriptive uncertainty, not a replacement for the
existing gates.

## Frozen gates and classification

The existing gates remain unchanged:

- between-subject: KS p-value greater than 0.05 and maximum quantile deviation
  less than 0.09;
- adaptive weighting: KS p-value greater than 0.05 and maximum quantile
  deviation less than 0.08.

The audit additionally labels a pooled arm:

- `anti_conservative` when the Wilson lower bound exceeds 0.05;
- `conservative` when the Wilson upper bound is below 0.05;
- `calibrated` when the Wilson interval contains 0.05 and both existing gates
  pass;
- `inconclusive` otherwise.

No failed gate will be rerun away or weakened. If a primary arm remains
anti-conservative, the package must either change the estimator and rerun this
exact plan under a new implementation identifier, or retain the result as a
versioned documented method limitation.

