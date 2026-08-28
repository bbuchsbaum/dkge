# Frozen functional-alignment court verdict

The frozen v1 court is complete and the candidate is **invalid for inferential
promotion**. The full permutation-re-estimation oracle passed its predeclared
gate, but the negative-control gate failed. No threshold, seed, arm, or solver
setting was changed after formal results.

## Execution and integrity

- 840 deterministic cohorts produced 5,280 raw arm rows: six frozen-plan arms
  in all 21 cells and the exact arm in six cells.
- All 840 seeds are unique and every required numeric result is finite. The
  3,600 `NA` values in `exact_elapsed_seconds` are structural: those rows are
  from cells where the exact arm was not scheduled.
- The run requested four R workers, but the managed process reported one
  physical core through `parallel::detectCores(logical = FALSE)` and therefore
  ran one R task at a time. Its numerical libraries used multiple threads. The
  actual elapsed time was 4,996.343 seconds versus the four-worker pilot
  projection of 1,929.44 seconds. This is a pilot extrapolation defect, not a
  change to the statistical workload.
- Formal source-tree, script, and protocol SHA-256 values recompute exactly to
  the values in `formal-manifest.txt`.

## Calibration result

- The exact full-re-estimation arm had FWER 0 or 0.05 in every scheduled cell;
  each Wilson 95% interval contained 0.05. Its minimum cell-level convergence
  was 0.9966.
- Geometry-only FWER ranged from 0 to 0.075 and no geometry cell had a Wilson
  lower bound above 0.05. Nevertheless, the frozen negative-control gate also
  required every cell's Wilson upper bound to be at most 0.12. At `n = 40`, a
  single rejection already has upper bound 0.1288, so the gate requires zero
  rejections in every control cell. It failed as written and is not weakened
  post hoc.
- Independently fitted *alignment features* did not rescue the same-data
  rank-truncated contrast estimator. FWER reached 0.425 under strong nuisance;
  six cells had Wilson lower bounds above 0.05. Independent correspondence data
  therefore cannot by itself be labelled inferentially eligible.
- The current fold-loading arm was materially inflated in numerous cells,
  reaching 0.35 under strong nuisance. The legacy and L1-only arms reached
  0.325. L1/L2 repairs are necessary provenance fixes but are not sufficient
  inferential repairs.
- The kernel-image residual prototype ranged from 0 to 0.10 and triggered no
  predeclared material-inflation flag. It remains approximate because its
  generating direction depends on the estimated rank-truncated fold span.
- The proposed `1 / (S * relative_eigengap)` predictor did not explain excess
  rejection in the two candidate same-data arms (linear coefficient 0.0058,
  p = 0.63; logistic coefficient 0.0674, p = 0.63). Its operator-sensitivity
  coefficient was -0.0202 (p = 0.0061), opposite the simple positive
  first-order prediction. The residual dependence is therefore not summarized
  as uniformly `O(1/S)` or as validated by this grid.
- The minimum cell/arm convergence was 0.953125 (geometry-only at low epsilon),
  so every cell met the frozen 95% solver gate. Individual cohorts did include
  sparse nonconvergence, which accounts for the accumulated mapper warnings.

## Consequences for FA-5 through FA-8

1. Exact same-data inference requires full re-estimation under every valid raw
   sign action for the frozen estimand.
2. `independent` describes the source of correspondence features, not a blanket
   exactness guarantee. Overall eligibility must also account for dependence in
   the contrast estimator and latent span.
3. `same_data_residualized` may be exposed only as approximate and opt-in. Its
   favorable court behavior cannot promote it to exact.
4. Geometry-only remains an ancillary correspondence bound, not evidence for
   functional alignment.
5. Group inference must refuse ineligible fitted alignments by default and make
   any approximate override explicit in the returned estimand/provenance.

Machine-readable evidence is in `formal-raw.csv`, `formal-summary.csv`,
`formal-inflation-flags.csv`, `formal-audit.csv`, the latent-span model files,
and `formal-manifest.txt`.
