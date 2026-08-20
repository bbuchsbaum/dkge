# DKGE vignette gauntlet progress

Benchmark: *Model to Meaning*, “Counterfactual comparisons”

## Quality bar

- A newcomer can state the input, task, output, and first valid workflow from the opening screen.
- Examples use the current exported API and execute against the checkout.
- Statistical claims distinguish estimands, weighting choices, estimability, and scientific limits.
- Visible output teaches interpretation; validation and build mechanics remain quiet.
- The vignette sequence has one clear entry point and purposeful next steps.

## Current state

- [x] Preserve and inventory the pre-existing dirty worktree.
- [x] Identify the 18 source vignettes and current pkgdown reading order.
- [x] Select the external comparison artifact.
- [x] Complete independent writer, critic, and API-verification passes.
- [x] Apply the first-round correctness, evidence, and navigation revisions.
- [x] Render and execute every changed vignette.
- [x] Build the source package and rebuild every vignette from the tarball.
- [x] Receive the final blind comparison and close all remaining blockers.

## Evidence log

- Source files, exports, package metadata, tests, and existing vignette diffs are under review.
- No generated `inst/doc` artifacts will be edited by hand.
- Baseline render: all 18 source vignettes rendered from an isolated package installation before the gauntlet edits; the dense-rendering article was slow but completed.
- Critic round 1 verdict: benchmark better. Main gaps were a false default-classification guarantee, an inverted CV-score explanation, unmatched-cluster averaging, unsupported adaptive-weighting and PLS claims, random attribution panels, and stale architecture documentation.
- Builder round 1 rewrote the onboarding path (`dkge`, `dkge-workflow`, and `dkge-concepts`). Targeted fixes now cover every high-severity finding from both independent audits.
- Critic round 2 still found overclaims around pooled anchors and CPCA, plus
  examples that printed output without tying it to one scientific estimand.
- Builder round 2 corrected the anchor preprocessing boundary, restricted CPCA
  claims to K-orthogonality, and revised classification and contrast examples
  into estimand -> number -> interpretation -> limitation chains.
- Focused renders passed for anchors, CPCA, classification, components,
  contrasts/inference, workflow, and the main introduction. The seeded contrast
  example recovered the planted positive direction in 8/8 subjects, with an
  active-region estimate of 0.808 and a conditional bootstrap interval of
  0.714 to 0.940. The seeded classification example produced accuracy 1.00 and
  log loss 1e-12, explicitly framed as strong toy-data recovery.
- `R CMD build --no-manual` succeeded and rebuilt all 18 source vignettes.
- `R CMD check --no-manual` passed tests, examples, package-vignette checks, and
  the tarball vignette rebuild. Final status: one environment/compiler warning
  (`-Wfixed-enum-extension` is unknown to the active Homebrew clang); the
  unavailable suggested package `T4transport` was reported as INFO.
- A final prose scan removed residual shorthand such as “fully
  cross-validated,” “fold-safe,” and solver-stability claims. The five affected
  articles rendered successfully after those edits; `git diff --check` remains
  clean.
- Final blind verdict: the DKGE suite now meets or beats the benchmark for its
  package purpose, with no blocker findings. The three remaining soft-polish
  suggestions were also implemented: seeded workflow values are interpreted in
  prose, the CPCA paragraph is limited to the output it prints, and the
  conditional bootstrap endpoints are restated and bounded. All three articles
  rendered successfully after the final edits.
