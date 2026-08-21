# Residual-space rotation calibration report

Originally executed 2026-08-16/17 and independently rerun from clean source on
2026-08-20 under the frozen plan in
`data-raw/dkge-between-rotation-plan.md`. The method, seeds, simulation arms,
gates, and promotion rule were fixed and versioned before the formal run. No
failed gate was relaxed or rerun with a replacement seed family.

## Certification status

The original run recorded source `HEAD 953d2a6` but did not capture its dirty
working-tree digest, so that run alone is historical, non-certifying evidence.
The complete frozen runner was rerun without modification in a disposable,
clean Git snapshot of the intended 2026-08-20 source. Snapshot commit
`5b185e9fd92e1e68a2984eaf33f9421a38fa51eb` identifies the exact rerun source;
it is deliberately not a commit in the primary checkout. All 8,100 replicate
rows and all 14 summary rows are byte-identical to the original uncompressed
artifacts, so the numerical decisions below are independently reproduced while
the missing historical digest is not fabricated.

## Decision

Residual-space rotation is now implemented as an opt-in method for
`dkge_between_permute()`, with a deliberately narrow validated scope:

- unweighted between-subject fits;
- one exchangeability block;
- at least two reduced-model residual dimensions; and
- a matrix-normal null with independent, homoscedastic subject rows and an
  arbitrary feature covariance for the finite-sample exactness claim.

It does **not** become the default. The frozen promotion rule required every
Gaussian null gate and the power gate to pass. The nonzero-nuisance interaction
arm rejected 33 of 500 nulls (0.066), one rejection above the predeclared 0.065
ceiling. Its Wilson interval contains 0.05 and its other uniformity diagnostics
pass, so this is a failed promotion gate rather than evidence of the material
inflation previously seen for Freedman-Lane. More decisively, rotation had
0.430 power in the frozen strong-interaction arm versus 0.580 for
Freedman-Lane: below the absolute 0.75 gate and 0.15 lower than the comparator,
exceeding the allowed 0.12 gap.

Freedman-Lane therefore remains the compatibility default and the only current
method for weighted or blocked analyses, but its documented small-sample
anti-conservatism remains. Users should choose rotation explicitly only when
its Gaussian row-sphericity assumptions are scientifically defensible, and
should expect a possible power cost at the tested rank and sample size.

## Size-adjusted comparison

The raw 0.430 versus 0.580 power comparison pits an exact test against one
whose matched-condition size is 0.086 (rotation size 0.050). At a
size-matched threshold, Freedman-Lane power falls to about 0.49, and the
predeclared relative gate would pass (`0.430 >= 0.367`). The frozen
absolute and relative gates are **not** relaxed: they still fail, and
rotation stays opt-in. Freedman-Lane's raw power advantage is mostly size
inflation. The matched-null diagnostic lives in
`dev/calibrate-dkge-between-rotation.R` rather than as a second shipped
artifact.

## Frozen evidence

| Arm | Term | Replicates | Rejection at 0.05 | Wilson 95% interval | KS p | Max quantile deviation | Gate |
|---|---|---:|---:|---:|---:|---:|---|
| Gaussian global null | group | 1,000 | 0.049 | 0.0373-0.0642 | 0.436 | 0.0100 | pass |
| Gaussian global null | trait | 1,000 | 0.055 | 0.0425-0.0709 | 0.560 | 0.0100 | pass |
| Gaussian global null | group:trait | 1,000 | 0.058 | 0.0451-0.0742 | 0.292 | 0.0225 | pass |
| Gaussian partial null | group | 500 | 0.040 | 0.0260-0.0610 | 0.811 | 0.0125 | pass |
| Gaussian partial null | trait | 500 | 0.042 | 0.0276-0.0634 | 0.134 | 0.0375 | pass |
| Gaussian partial null | group:trait | 500 | 0.066 | 0.0474-0.0912 | 0.069 | 0.0406 | **fail** |
| Multivariate t5 global null | group | 500 | 0.054 | 0.0374-0.0774 | 0.844 | 0.0150 | pass |
| Multivariate t5 global null | trait | 500 | 0.052 | 0.0357-0.0751 | 0.536 | 0.0213 | pass |
| Multivariate t5 global null | group:trait | 500 | 0.060 | 0.0423-0.0844 | 0.536 | 0.0150 | pass |
| Centered-skew global null | group | 500 | 0.038 | 0.0245-0.0586 | 0.008 | 0.0694 | pass |
| Centered-skew global null | trait | 500 | 0.064 | 0.0457-0.0890 | 0.008 | 0.0663 | pass |
| Centered-skew global null | group:trait | 500 | 0.066 | 0.0474-0.0912 | 0.416 | 0.0238 | pass |

The non-Gaussian gates were intentionally wider and do not establish exactness
outside the Gaussian model. The complete summary also contains the paired
power arm. All 8,100 term-level results are shipped in
`inst/extdata/dkge-between-rotation-replicates.csv`.

## Implementation evidence

- The compressed evaluator matches explicit dense reconstruction and refitting
  to `1e-10`.
- Haar draws are orthogonal to `1e-12`; rotated data preserve
  `crossprod(Y)` to `1e-10`.
- Seeded serial and future-based parallel paths are bit-identical and restore
  the caller's random-number state.
- Unsupported weights, multiple blocks, missing block labels, and residual
  dimension below two fail before resampling.
- The implementation stores only effect/feature cross-products and a
  residual-dimension-square rotation; it does not form a feature-square
  covariance.

## Provenance

- Mote ticket: `bd-01M06PBXPPMQR37M52PFMCDDXQ`
- Frozen plan SHA-256:
  `95c3258340128be36e99675512a051f6fd3abd8eccf3baf187657cba15a6e735`
- Frozen runner SHA-256:
  `60a99d76e813e0a164c8323a0e664fd3e8a2ce5c8c0a49f61de37f9159ef7736`
- Replicate artifact SHA-256:
  `c1be2a843e7bd06c4d98265fa249a6a660734a1099288f05b8413074259ead3d`
- Summary artifact SHA-256:
  `efef6f51cd4f9c0e9429aac191cce570c50e07824d106edf465618f3aa634b95`
- Metadata artifact SHA-256:
  `e267ab454cb59dd4b977f0a75ef64ca8cd2b4d8119ed4340e192a9b5db5de90c`
- Exact clean-source snapshot commit:
  `5b185e9fd92e1e68a2984eaf33f9421a38fa51eb`
- Historical source record: `HEAD 953d2a64ff356d55f4a6fd3278a7f11db018fea6`;
  working-tree digest not captured, explicitly non-certifying on its own.
- The replicate CSV is packaged as `.csv.gz`; its uncompressed content is
  byte-identical to the clean rerun and has the SHA-256 above.
- Clean-rerun runtime: 663 seconds with R 4.5.1 on
  `aarch64-apple-darwin20`.
