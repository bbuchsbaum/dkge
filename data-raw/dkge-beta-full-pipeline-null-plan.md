# Predeclared DKGE public-beta full-pipeline null calibration

Status: frozen before the formal run for Mote
`bd-01M0GADXACGVAS314W6P7C7HDY`.

## Question and scope

This court asks whether the public-beta inferential paths remain calibrated
when the fitted representation, cross-fitting, transport, classification, and
multiplicity machinery are exercised together under a generated global null.
It is a narrow software-and-method calibration, not a universal validity claim
for neuroimaging data.

The formal run must start from a clean Git commit. The runner records the exact
40-character source commit, SHA-256 digests of this plan and the runner, the R
and package versions, platform, seeds, and all frozen settings. Output is
written only after the clean-source check, so generated evidence may be
committed separately without pretending it was part of the evaluated source.

## Frozen data-generating process

- 200 independent null data sets.
- Eight subjects, three effect rows, and subject parcel counts alternating
  between three and four.
- Every beta entry is independent standard Gaussian noise. Effect rows are
  exchangeable and have no class or contrast signal.
- The design kernel is `diag(c(1, 1, 0))`. A requested rank of three must be
  capped to the two-dimensional kernel range; every leave-one-subject-out
  training fit must retain rank two.
- Subject centroids are deterministic functions of parcel index and subject.
  They do not depend on beta values or randomized signs.

Replicate seeds are `8600000 + replicate_id * 1000`. Inference and
classification consume disjoint offsets within each replicate.

## Arm 1: transported contrast inference

- Contrast: `c(1, -1, 0)`.
- LOSO cross-fitting with fold alignment enabled.
- Geometry-only Sinkhorn transport to subject 1 with `lambda_emb = 0`,
  `lambda_spa = 1`, `epsilon = 0.2`, and no warm start.
- Provenance class: `geometry_only`. The fixed operator is sign-invariant by
  construction; loading distance has exactly zero cost weight.
- 199 sign flips and max-T correction across all medoid parcels.
- Replicate rejection: any max-T adjusted parcel p-value at or below 0.05.

## Arm 2: prospective classification

- Two direct cell targets select effect rows 1 and 2.
- `mode = "cell_cross"`; each representation excludes the held-out subject.
- LDA with one externally preselected `lambda = 0.001` and 199 within-subject
  label permutations. No observed-data lambda grid or lambda-selection callback
  is permitted.
- Every randomized label configuration permutes the corresponding subject
  effect rows, refits DKGE (including numerical rank selection), rebuilds all
  outer-fold representations, refits the classifier, and reevaluates log loss.
  Duplicate configurations may reuse an exact cached result within a replicate.
- Metric: log loss, a continuous proper score that avoids the severe tie
  discreteness of accuracy with only 16 held-out cell rows.
- Replicate rejection: the empirical log-loss p-value at or below 0.05.

The lambda-selection policy is itself frozen: the penalty is selected outside
the inferential sample and held fixed. The public API must reject a missing
lambda, `lambda_grid`, or `lambda_fun` whenever permutations are requested.

## Cross-arm multiplicity

For each replicate, the transported max-T p-value and classification p-value
are combined with a two-test Bonferroni correction. Family rejection occurs
when `min(1, 2 * min(p_transport, p_classification)) <= 0.05`. This is in
addition to max-T correction within the transported spatial family.

## Frozen gates

For each arm separately:

- empirical rejection at 0.05 must lie in `[0.015, 0.085]`;
- the two-sided 95% Wilson interval must contain 0.05; and
- no replicate may violate the rank, provenance, finite-result, or prospective
  claim-scope contract.

For the Bonferroni family:

- empirical rejection must not exceed 0.075; and
- the Wilson lower bound must not exceed 0.05.

All 200 replicates must complete. A failed gate is reported as a failure; no
seed family, threshold, arm, or failed replicate may be replaced after results
are observed. Reduced pilot runs may be used only for runner debugging and may
not be shipped as certification evidence.

## Pre-formal revision record

The first 20-replicate runner diagnostic used accuracy and produced no
rejections with a mean p-value of 0.716. Inspection showed that accuracy has
only 17 attainable values for 16 held-out rows, so inclusive tie handling makes
the randomization p-value strongly conservative. Before any formal run, the
primary classification metric was changed to continuous log loss. The data
generator, exact full-pipeline recomputation, seeds, replicate count,
randomization count, thresholds, and all other gates were unchanged. The
accuracy pilot is non-certifying and remains recorded in the Mote history.

## Artifacts

- Runner: `tools/calibration/run-beta-full-pipeline-null.R`
- Replicates: `inst/extdata/dkge-beta-full-pipeline-null-replicates.csv.gz`
- Summary: `inst/extdata/dkge-beta-full-pipeline-null-summary.csv`
- Metadata: `inst/extdata/dkge-beta-full-pipeline-null-metadata.csv`
- Interpretation: `data-raw/dkge-beta-full-pipeline-null-report.md`
