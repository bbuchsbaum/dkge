# Predeclared 3 x 5 x 4 unbalanced-trial validation

Status: frozen before the first simulation run for Mote ticket
`bd-01M03KTXAZ6PR3GS62PKGPCG2V`.

## Reproducibility contract

- Canonical command: `Rscript dev/simulate-dkge-unbalanced-trialwise.R`.
- Replicate seeds: `73101:73120` (20 replicates per scenario).
- The tracked summary is written to
  `inst/extdata/dkge-unbalanced-trialwise-summary.csv`; replicate-level rows
  are written to `inst/extdata/dkge-unbalanced-trialwise-replicates.csv`.
- `DKGE_SIM_REPS` may reduce the number of leading seeds for local smoke tests,
  but only the full 20-seed run may be used for the gates below.
- This file, including thresholds, must be committed before results are
  interpreted. Thresholds are not changed after inspecting results.

## Data-generating design

Each dataset has 18 subjects, 12 spatial features, and a 3 x 5 within-subject
factorial crossed with ordinal response levels 1, 2, 3, and 4 (`q = 60`). Trial
counts are subject- and cell-specific Poisson draws with subject exposure and
response-level modifiers. Counts can be zero. Only observed cells enter each
subject's local one-hot second-stage design; `dkge_data()` aligns them into the
named 60-effect union.

Four independent acquisition runs are generated. Noise is Gaussian and is one
of:

1. `iid`: common trial variance;
2. `heteroscedastic`: response- and subject-dependent trial variance;
3. `ar`: AR(1), rho = 0.60, independently within run;
4. `unequal_precision`: feature and subject variances spanning a four-fold
   range, combined with cell-count imbalance.

Missingness is either MCAR cell dropout or response-dependent dropout, with
larger empty-cell probability for response 4. Every retained subject has at
least two observations per represented local effect and a full-rank local
one-hot design.

Signal regimes are:

- `null`: no systematic effect signal;
- `ordinal`: a centered linear response trend;
- `nonlinear`: a U-shaped response profile orthogonal to the linear trend;
- `interaction`: response trend whose sign and magnitude depend on both
  factorial variables;
- `low_rank`: two independent effect profiles and spatial loading vectors.

The fixed ten-scenario grid crosses every signal, noise, and missingness regime
without claiming a complete Cartesian design:

| scenario | signal | noise | missingness |
|---|---|---|---|
| 1 | null | iid | MCAR |
| 2 | null | heteroscedastic | response-dependent |
| 3 | ordinal | iid | MCAR |
| 4 | ordinal | heteroscedastic | response-dependent |
| 5 | nonlinear | ar | MCAR |
| 6 | interaction | ar | response-dependent |
| 7 | low_rank | unequal_precision | MCAR |
| 8 | ordinal | unequal_precision | response-dependent |
| 9 | nonlinear | heteroscedastic | response-dependent |
| 10 | low_rank | ar | MCAR |

## Methods

Every dataset is analyzed with the same named ordinal design kernel and rank
equal to the true signal rank (rank one under the null, solely to define the
noise diagnostic):

1. `legacy`: cell means, no effect weighting or debiasing;
2. `count`: precision proportional to cell count;
3. `explicit_precision`: count divided by the known cell-average trial
   variance;
4. `random_effects`: DL random-effects precision with the documented cap;
5. `iid_analytic`: OLS cell means and IID analytic noise subtraction;
6. `covariance_analytic`: GLS cell means using the known trial covariance and
   analytic subtraction;
7. `split_half`: run-disjoint cross-half moments with split-derived precision.

All methods use pair-normalized pooling (`missingness = "rescale"`) where
coverage is incomplete. No method receives the true signal coefficients.

## Metrics

Metrics are computed per dataset and method, then summarized by the median
across the 20 fixed seeds unless stated otherwise.

- `subspace_angle_deg`: largest principal angle between the estimated and
  noiseless target subspaces in the K-whitened coordinate.
- `inverse_count_correlation`: absolute correlation between squared leading
  K-whitened loading and mean inverse observed cell count.
- `contrast_bias`: signed error in the normalized quadratic response-trend
  estimand, `c' M c / tr(M)`.
- `contrast_rmse`: root mean squared contrast error across seeds.
- `interval_coverage`: percentile subject-bootstrap coverage of the noiseless
  response-trend estimand (399 fixed resamples per dataset).
- `false_positive_rate`: proportion of null datasets with a sign-flip
  response-trend p-value <= .05 (999 flips).
- `high_count_influence`: maximum subject share of the response-trend pair
  weight, averaged over cells with nonzero support.

## Predeclared gates

The full run passes only if all confirmatory gates hold:

1. In AR and unequal-precision non-null scenarios, `covariance_analytic` has
   median subspace angle and contrast RMSE no greater than 90% of `legacy`.
2. Across non-null scenarios, the median inverse-count correlation for each of
   `explicit_precision`, `random_effects`, `covariance_analytic`, and
   `split_half` is at least 0.05 below `legacy`, or is already <= 0.10.
3. Across non-null scenarios, absolute median contrast bias for
   `covariance_analytic` and `split_half` is no larger than for `legacy`.
4. Null false-positive rate for supported sign-flip inference lies in
   `[0.01, 0.12]` for `covariance_analytic` and `split_half`.
5. Interval coverage lies in `[0.80, 1.00]` for `covariance_analytic` and
   `split_half` in non-null scenarios.
6. In high-imbalance scenarios, `random_effects` maximum subject influence is
   no greater than 90% of `count` influence.
7. In IID non-null scenarios, `iid_analytic` and `covariance_analytic` median
   subspace angles differ by no more than 5 degrees.

All scenario-level rows and gate outcomes are retained. A failed gate is a
scientific limitation to report, not permission to revise this plan.
