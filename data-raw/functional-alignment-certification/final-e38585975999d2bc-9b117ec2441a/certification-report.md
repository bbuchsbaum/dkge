# Functional-alignment certification

- Status: `complete_approximate_only`
- Source-tree SHA-256: `e38585975999d2bc6bdc1b76747e90a3ab89328940c92d676592494255831400`
- Tarball SHA-256: `9b117ec2441a4d6d94d23c39b977828aab8a5fea6c7ce09b96ef741d4c5cbc59`
- Protocol SHA-256: `4e6b828f463d9ae57de5a100f440dc7e96a376da3fdcae34509a42362223ac0a`
- Inferential promotion: **blocked by the frozen court**
- Scope: exact candidate only; no merge, tag, hosted publication, or release claim

## Gate result

All provenance, deterministic-execution, algebraic/regression, source-test, built-tarball, documentation, and independent-review gates passed. The exact full-pipeline re-estimation exact null-action comparator passed. The negative-control promotion gate did not pass, so aligned inference remains explicitly approximate and fail-closed.

## Formal court: baseline

The superseded v2 protocol, failed-run budget, source-test manifest, results, and session receipt are retained under evidence bundle `ba8afb6e29f53fd705646086a9eaabcaff161c4793ea0b71eb591332281e8cfa`. V3 changed only the Sinkhorn iteration ceiling and its directory label. V3 is preserved under history bundle `d5f255d6234f28acc543b7ff5cd19e78cf8d5958c4078258fc4bb144a5b215f3`. V4 completed computation after its candidate drifted, but its analyzer failed before reading statistical outputs; the complete inadmissible history is bound under `12b64ae352a266f897b30e2ceedf17dea3ba428ddaed08a1f85f9d084c3e08bf`. V5 and V6 were superseded before execution under `819f0f5e6e78f3abe7c8d8f508cffac9fb07ea348aa1cd34ae2f1f069123afd4` and `88c0e583d3012b75839e766d746d20fd2c58934026301b9c873ed9e3de30dae5`. V7 completed source, court, and power execution but was superseded before collection after the exact package check exposed nonstatistical documentation and test-harness defects; its protocol and supersession are bound under `cf53ca7ab03080f5e952ea8c2dc0980f05af5e806b0d5b674f6ebc6f0874d04a`. V8 completed source, court, power, and exact package-check execution but failed the separately required pkgdown gate before collection; every V8 output and the failure supersession are bound under `8e89ab1b2a55a4bfac61b1fdfcb773a3ec0ec3d1840128bed8ce5cbc7d56816f`. V9 changes only the 19 vignette comparator namespaces, validation and supersession provenance, and the output-directory label.

- `full_permutation_reestimation`: FWER 0.050 (Wilson 95% CI 0.014-0.165)
- `geometry_only`: FWER 0.050 (Wilson 95% CI 0.014-0.165)
- `independent_alignment`: FWER 0.075 (Wilson 95% CI 0.026-0.199)
- `kernel_image_residual_prototype`: FWER 0.075 (Wilson 95% CI 0.026-0.199)

The complete factor grid, uncertainty intervals, adverse cells, and latent-span
models are retained in the court run directory; no observed result was removed
and no threshold was changed after inspection.

## Held-out known-warp efficacy

- `geometry_only`: power 0.025 (Wilson 95% CI 0.004-0.129), correlation -0.016, RMSE 0.371, amplitude ratio 1.150, point spread 1.000
- `iterative_template`: power 0.025 (Wilson 95% CI 0.004-0.129), correlation 0.225, RMSE 0.373, amplitude ratio 1.452, point spread 1.033
- `ordinary_mni_average`: power 0.025 (Wilson 95% CI 0.004-0.129), correlation -0.016, RMSE 0.371, amplitude ratio 1.150, point spread 1.000
- `raw_functional_medoid`: power 0.075 (Wilson 95% CI 0.026-0.199), correlation 0.279, RMSE 0.336, amplitude ratio 1.279, point spread 1.170

These power and recovery metrics are descriptive. They cannot override the
failed Type-I promotion gate.
Template outer-loop convergence was 39/40; every scheduled cohort, including non-convergence, is retained.

## Package evidence

- Full source suite: 4338 expectations, 0 failures; 14 warning events recorded with generic classifications in the source-test manifest; 4 manifest-counted skips
- Known-warp serial/parallel comparison: identical within 1e-12
- Built-tarball R CMD check: `Status: OK`
- R CMD check policy: vignettes excluded only through `--ignore-vignettes`; installation, tests, examples, and codoc retained
- Checked package source is byte-bound to the tarball contents
- pkgdown required pages: 7/7 present and byte-bound to the exact tarball
- Independent review verdict: `pass_with_limitations` with 0 blockers

See `certification-manifest.json`, `artifact-inventory.csv`, and
`evidence-map.json` for the machine-readable record.
