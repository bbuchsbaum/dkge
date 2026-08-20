# DKGE public-beta full-pipeline null calibration report

## Decision

The preregistered formal calibration passed. All 200 generated null data sets
completed under the frozen 199-randomization design, every structural contract
passed, and all three statistical gates passed. These results support the two
tested public-beta paths; they do not establish universal calibration for other
data-generating processes, transport provenance classes, classifiers, tuning
policies, or inferential procedures.

## Exact evaluated source

- Source commit: `3fdbac18ef9b5db682f6d27658421ee0fb972e70`
- Plan SHA-256: `8f1cbb6549c17801a159f1db16e5bcf7b9ddd437b743a15d7f789ea917c986e3`
- Runner SHA-256: `297ee165d414de8d371bb3165783233a035dfd9c42740682eb05a6e05ef049da`
- R: 4.5.1 on `aarch64-apple-darwin20`
- Started: 2026-08-20 21:52:51 UTC
- Completed: 2026-08-20 22:02:40 UTC
- Runtime: 588.526 seconds

The runner required a clean worktree before beginning and wrote the evidence
only after resolving the exact 40-character source commit.

## Statistical results

| Arm | Rejections | Rate | 95% Wilson interval | Frozen rate gate | Result |
|---|---:|---:|---:|---:|---|
| Geometry-only transported max-T | 5/200 | 0.025 | [0.0107, 0.0572] | [0.015, 0.085] and interval contains 0.05 | Pass |
| Prospective `cell_cross` log loss | 9/200 | 0.045 | [0.0239, 0.0833] | [0.015, 0.085] and interval contains 0.05 | Pass |
| Two-arm Bonferroni family | 5/200 | 0.025 | [0.0107, 0.0572] | [0, 0.075] and lower bound <= 0.05 | Pass |

Mean p-values were 0.5295, 0.5111, and 0.6050, respectively. The formal run
did not alter seeds, thresholds, metrics, arms, or failed replicates after the
results were observed.

## Structural contracts

All 200 replicates satisfied every frozen contract:

- singular-PSD kernel rank and effective rank were both two;
- every held-out training fit retained rank two;
- fold alignment metadata was present;
- transported inference used `geometry_only` provenance and reported
  `randomization_exact_sign_invariant_operator`;
- classification reported a prospective held-out-subject claim with a basis
  cross-fitted without the held-out subject;
- the externally preselected penalty remained exactly 0.001; and
- all transported and classification p-values were finite.

## Archived machine-readable evidence

- `inst/extdata/dkge-beta-full-pipeline-null-replicates.csv.gz`
  - SHA-256: `4147c5d68e96ea9a115882fa55c7640ac6096a8ce8f2928b292fb2169d2751b6`
- `inst/extdata/dkge-beta-full-pipeline-null-summary.csv`
  - SHA-256: `185863f63ca88880091aaf6616f142482d43edd5fe74d11202f515e104fd95ac`
- `inst/extdata/dkge-beta-full-pipeline-null-metadata.csv`
  - contains the exact source, plan, runner, environment, settings, contract
    result, and artifact digests.

The committed plan remains the authoritative statement of the null generator,
full-pipeline recomputation boundary, multiplicity policy, and interpretation
limits.
