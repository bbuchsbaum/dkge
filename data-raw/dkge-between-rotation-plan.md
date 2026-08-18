# Predeclared residual-space rotation inference validation

Status: frozen before implementation results for Mote ticket
`bd-01M06PBXPPMQR37M52PFMCDDXQ`.

## Method contract

For a term-specific reduced design `X0`, let `Q0` span `col(X0)` and let `Qe`
span its orthogonal complement. A rotation replicate is

`Y* = Q0 Q0' Y + Qe G Qe' Y`,

where `G` is Haar-uniform on the orthogonal group in the residual dimension.
The observed and rotated samples use the existing reduced-minus-full RRR
residual-SSE statistic at the fitted rank. The implementation may compute the
same quantities from compressed cross-products, but it must not change the
statistic.

Under a matrix-normal null with independent homoscedastic subject rows and an
arbitrary feature covariance, the transformation fixes the null mean and
preserves the error law. It also preserves `crossprod(Y)` exactly. Normality and
row sphericity are method assumptions, not conclusions drawn from the data.

Version 1 supports unweighted fits and a single exchangeability block. It must
fail closed for subject/feature weights, multiple blocks, or fewer than two
residual dimensions. Freedman-Lane remains available explicitly.

## Implementation gates

- A dense reconstruction oracle and the compressed implementation agree for
  full/reduced scores, total sum of squares, statistics, feature statistics,
  and p-values to `1e-10` relative/absolute tolerance.
- Every generated `G` satisfies `crossprod(G) = I` to `1e-12` at ordinary
  dimensions; rotated data preserve `crossprod(Y)` to `1e-10`.
- Seeded serial and parallel paths are exactly equal and restore the caller's
  RNG state.
- Row/feature naming, rank tolerance, `+1` p-value correction, and maxT/FDR
  behavior are unchanged.
- No subject-by-subject or feature-by-feature square matrix is formed; the
  largest new square object is residual-dimension by residual-dimension.

## Frozen simulations

All primary calibration arms use `n = 20`, `p = 24`, RRR rank 2, 399 rotations,
five deterministic batches, and the same `~ group * trait + age` design used by
`null-calibration-v1`. Seeds are disjoint from all prior calibration families.

1. **Gaussian global null:** five batches of 200 data sets, base seeds
   1,700,000 through 2,100,000 by 100,000; 1,000 p-values per term.
2. **Gaussian partial null with nonzero nuisance signal:** separately for
   `group`, `trait`, and `group:trait`, five batches of 100 data sets per term.
   The true coefficient matrix has rank at most 2, the tested term is exactly
   zero, and other eligible design rows carry signal. Seed families begin at
   2,300,000, 2,900,000, and 3,500,000.
3. **Non-Gaussian robustness:** multivariate t with 5 degrees of freedom and a
   centered-skewed innovation arm, each using five batches of 100 global-null
   data sets. Seed families begin at 4,100,000 and 4,700,000. These are
   robustness tests, not exactness claims.
4. **Strong-signal power:** 300 data sets with a rank-one interaction
   coefficient of 1.0, 199 resamples, seed family beginning at 5,300,000.
   Rotation and Freedman-Lane use the same generated data.

## Frozen diagnostics and gates

For every calibration arm report empirical size at 0.05, Monte Carlo SE,
Wilson 95% interval, mean/median p, KS D and p, batch range, five p-value
quantiles, and maximum quantile deviation.

Primary Gaussian promotion gates, applied separately to every term:

- empirical size in `[0.035, 0.065]`;
- Wilson interval contains 0.05;
- KS p-value greater than 0.01;
- maximum quantile deviation below 0.06.

Non-Gaussian robustness gates:

- empirical size in `[0.025, 0.075]`;
- Wilson lower bound does not exceed 0.05;
- KS p-value greater than 0.005;
- maximum quantile deviation below 0.08.

Power gate: rotation rejection at 0.05 is at least 0.75 and no more than 0.12
below Freedman-Lane. This comparison is descriptive because Freedman-Lane's
null size is already inflated.

## Promotion rule

If every Gaussian global-null and nonzero-nuisance gate passes, rotation becomes
the default `dkge_between_permute()` method for its supported unweighted,
single-block scope. Otherwise it remains opt-in. Non-Gaussian failure narrows
the documented scope but does not overturn Gaussian exactness; no failed gate
will be rerun away or relaxed after results are observed.

## Artifacts

- Runner: `dev/calibrate-dkge-between-rotation.R`
- Replicates: `inst/extdata/dkge-between-rotation-replicates.csv`
- Summary: `inst/extdata/dkge-between-rotation-summary.csv`
- Metadata: `inst/extdata/dkge-between-rotation-metadata.csv`
- Interpretation: `data-raw/dkge-between-rotation-report.md`
