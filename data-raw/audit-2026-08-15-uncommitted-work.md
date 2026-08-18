# Audit of uncommitted work — 2026-08-15

Scope: everything in the working tree since `953d2a6` (35 modified files,
+2432/−264; new `R/dkge-aggregate.R`, `dkge-moments.R`, `dkge-effect-grid.R`,
`dkge-effect-weights.R`, `dkge-trial-subject.R`, three new test files, one new
vignette). Four independent reviewers each took a disjoint slice; every
CONFIRMED item below was reproduced by running code against the working tree,
and the blockers were re-checked by hand in source.

Baseline: `devtools::test()` → 1483 pass / 0 fail / 48 skip.
`R CMD check` (no tests/vignettes) → 0 errors, 2 NOTEs (see §6), plus a
local-toolchain install WARNING that is not package code.

## 0. Verdict

The **core linear algebra is right** and was verified to machine precision:

- New moment-pooling engine reproduces the legacy `Chat` exactly on the
  default path (max |Chat − Σ w_s contrib_s| = 2e-13; `UᵀKU = I`).
- `dkge_aggregate_fit` matches `svd(K^{1/2} Y)` to 1e-15.
- Closed-form Freedman–Lane between-subject permutation matches the old refit
  loop to 1e-13 (weighted and unweighted, with exchangeability blocks) at
  ~5× speed.
- Folds rewrite reproduces legacy fold `Chat` (0 diff full-train, 1e-13 subset).
- Analytic noise trace, split-half cross-moment, Kish ESS, and
  `dkge_trial_subject` sufficient statistics all match hand-computed oracles.

The **problems are at the seams** — resampling bookkeeping, partial-coverage
conventions, and the new plotting layer. Nine findings produce silently wrong
numbers on reachable inputs. The suite is green because the tests that touch
those paths either re-derive the implementation or exercise only the case where
the bug is invisible. **Do not commit the resampling / plotting / partial-
coverage layers as-is.**

## 1. Blockers (silent wrong answers, reachable via defaults)

| # | Where | What | Fix |
|---|-------|------|-----|
| B1 | `R/dkge-aggregate.R:839` | `sample(ii, length(ii), replace=TRUE)` with a length-1 stratum samples `1:ii`, not `ii`. Default `strata = interaction(group_vars)`, so any singleton cell mis-stratifies the bootstrap (verified: singleton stratum drew subject 3 from stratum "a"; 6 singletons → hard error). | `ii[sample.int(length(ii), length(ii), TRUE)]`; error/warn on strata < 2. |
| B2 | `R/dkge-aggregate.R:641` | After K-Procrustes alignment, `U`, `saliences`, `scores_feature` are rotated but `singular_values` are only sliced. Rank-3 bootstrap: 27/30 draws wrong by >5%; swap case reports 5 where truth is 1. Contaminates `component_score_contrast`, `between_group_contrast`, `component_scale="score"`, both permutation null and bootstrap. | `out$singular_values <- sqrt(colSums(out$scores_feature^2))` (equals `d_k` exactly when unaligned). |
| B3 | `R/dkge-aggregate.R:563` | `abs()` on every contrast statistic → bootstrap percentile interval of a folded quantity can never contain 0. On pure-null data: interval [0.02, 1.67], "excludes zero" = TRUE. Undocumented; inconsistent with the signed `excludes_zero` at :979. | Return signed; take `abs` only inside the permutation p-value. |
| B4 | `R/dkge-fit.R:186-207`, `R/dkge-data.R:513` | Duplicate subject IDs are not rejected; `.dkge_obs_masks_from_provenance` reorders by `match(subject_ids, names(masks))` so both duplicates take subject 1's mask. Verified: subject 2's observed row silently zeroed (`effect_moment[3,3]` = 0, truth 7.15). | Reject duplicates in `.normalize_subject_ids`; reorder only when names unique. |
| B5 | `R/dkge-fit-core.R:596-601` | Masked X blocks: `block[obs,] <- Khalf[obs,obs] %*% …` breaks `Chat == X Xᵀ` (diff 47.5, 6.4% rel) and desyncs training blocks from `dkge_transform_block()` (diff 1.15) under partial coverage. Multiblock loadings and out-of-sample projection are inconsistent with the fit. Also masks rows of `Btil = RᵀB`, which are already R-mixed. | Revert to `Khalf %*% block` (betas already zero-filled upstream). |
| B6 | `R/dkge-fit.R:369` | `missingness="shrink"` blends `rescaled` with `diag(Chat)` (un-rescaled). Cell 2,2: reports 7.2, should be 4. Twin in `dkge-moments.R:224` does it right → the two shrink paths disagree. | `diag_part <- diag(diag(rescaled), …)`. |
| B7 | `R/dkge-plot.R:82-100` | `dkge_principal_angles_K` (exported) assumes K-orthonormal inputs, doesn't check or document it; `pmin(sv,1)` clamp turns garbage into "perfectly aligned". `U2 = 0.5·U1` → 60°. | K-orthonormalize inputs or `stop()` when `UᵀKU ≠ I`. |
| B8 | `R/dkge-plot.R:567-572` | `as.character(groups)` strips names, so a named group vector is matched by position; verified labels exactly inverted for reverse-ordered input. data.frame branch yields NA for missing subjects with no warning. | Preserve names; error on incomplete coverage. |
| B9 | `R/dkge-plot.R:610-690` | `dkge_subject_component_projections` returns `‖BwᵀKu‖²` (non-negative, sign-blind: negated betas give identical values) while docs and `geom_hline(0)` promise a signed score. Pooled mode is bit-identical to existing `.dkge_energy_by_subject_component()` (2.8e-14) computed more expensively. | Either compute a signed brain score or rename to energy and delegate. |

## 2. Design decision needed: partial-coverage masking convention

Three mutually inconsistent conventions for "masked" now coexist:

1. **Fit path** (`dkge-moments.R:144-149`): zero-fill unobserved rows/cols in
   raw effect space, then apply the **full** `Khalf` congruence.
2. **Subject weights + dead accumulator** (`dkge-fit.R:136, 255`): `Khalf[obs,obs]`
   — which is *not* a root of `K[obs,obs]` and operates on R-mixed rows.
3. **`dkge-project.R:46,67`**: unmasked.

Reviewers disagreed on which is "right". My recommendation: **convention 1**.
It is what the docstring states ("policies act in raw effect space, before R
and K mix rows") and it is the only one that is mathematically coherent — the
K-metric embedding of a masked effect-space covariance. The consequence to
accept and document is that a *coupling* K (no `block_factors`) will spread
observed energy into unobserved coordinates; that is the metric doing its job,
not a leak. But two things then follow:

- `.dkge_subject_weights` (`dkge-fit.R:131-137`) must switch to the same
  convention (its result already differs from pre-diff under partial coverage:
  1.5318 vs 1.5276), and the empty-mask branch `weights[s] <- 1/1e-12` must go
  (an all-unobserved subject currently gets the *largest* weight, verified
  w = (0.3, 0.3, 2.4)).
- `vignettes/dkge-partial-effect-spaces.Rmd:140-142` claims placeholders "do
  not … leak through `K^{1/2}`". That is true only for block-diagonal K (the
  vignette's own example). Reword.

If instead you want strict non-propagation, re-apply the mask to `Chat`
after the transform — but then say so, and drop `Khalf[obs,obs]` anyway.

## 3. Significant (wrong or misleading, narrower reach)

- `R/dkge-moments.R:207-233`: `missingness="rescale"` is a silent no-op under
  `effect_weights` (verified `identical(Chat)`); other policies dispatch.
- `R/dkge-moments.R:20, 57, 135`: `x[-integer(0),] <- 0` trap — all-FALSE mask
  returns the *unmasked* moment. Unreachable via `dkge_data()` today, reachable
  from `.dkge_fold_weight_context`.
- `R/dkge-fit-core.R:103`, `R/dkge-align-data.R:22-28`: kernel dimnames never
  validated against `dat$effects`; `dkge_data()` unions effects in first-
  appearance order and has no `effects=` to pin to `grid$cell_labels`.
  Reordering the subject list silently reindexes the kernel. Now much likelier
  to bite because `dkge_effect_grid` exists precisely to define an order the
  union doesn't reproduce.
- `R/dkge-contrast.R:179-213`: estimability classification is name-matching
  only; the `blocks` structural fallback is dead (real cell-basis kernels have a
  single block named "cells"). A between-subject contrast named anything but
  the term name → "unknown" → no warning. `kernel_info$factor_scope` + `cells`
  already carry what's needed to classify structurally.
- `R/dkge-aggregate.R:603-609`: `dkge_aggregate_align` never checks row
  identity; `:490` skips label check when K has rownames but no colnames;
  `:706-716` permuting a subset of `group_vars` can change the joint row set
  mid-run (B=200 died with "Kernel dimensions must match"); `:302` NA in
  grouping vars gives an unrelated error.
- `R/dkge-plot.R:167-232, 332-360`: `dkge_design_basis(normalize="unit_l2")`
  normalizes Euclidean but scores are `CᵀKU`; grand-mean coordinate inflated
  ~2.45× vs interaction terms purely by the metric. `saliences(scale="unit")`
  divides an already K-unit vector by its Euclidean norm. Add `unit_K`.
- `R/dkge-contrast-validated.R:108`: `merge` on `contrast` cartesian-blows-up
  with duplicate names (verified 4 rows for 2 contrasts); `cbind` suffices.
- `R/design-kernel.R:185, 199-200` / `R/dkge-effect-grid.R:27-42, 104-106`:
  typo'd `terms` name → "missing value where TRUE/FALSE needed" (previously
  silently ignored); `list(L=3, levels=c("a","b"))` → dimnames error;
  continuous list spec silently coerced to a 2-level nominal factor.
- `R/dkge-effect-grid.R:65-71` duplicates `.dkge_effect_grid_cell_labels`
  line-for-line.
- `R/dkge-fit.R:51` vs `R/dkge-procrustes.R:5`: `.dkge_kernel_roots` defined
  twice with different third elements; collation order decides. Pre-existing
  but now load-bearing.
- Backward compat: `dkge_fit()`, `dkge()`, `dkge_subject_model()` gained args
  mid-signature; missingness defaults in `dkge_loso_contrast`/`dkge_cv_*` now
  inherit from the fit (correct, but unannounced — NEWS entry).

## 4. Efficiency

- **`R/dkge-aggregate.R:263-321`** — `.dkge_aggregate_from_spec` builds n·q
  one-row data frames + string-mask scans per row. Measured 0.93 s/draw at
  n=100, q=20 → ~16 min of pure aggregation for B=999. A `rowsum()` version on
  a once-stacked matrix is bit-identical and ~400× faster. Better still,
  precompute the Gram `G = VVᵀ` once so permutation never touches P.
- **`R/dkge-bootstrap.R:183-194`** — `.dkge_repool_fit()` returns non-NULL for
  every current fit, so the O(1) matvec branch is dead and every replicate
  re-pools S q×q moments: 129× slower (0.645 s vs 0.005 s / 500 reps, S=40,
  q=20) even when pooling is provably linear. The right predicate
  (`nonlinear_pool`) already exists in `dkge-bootstrap.R:303` / `dkge-analytic.R:78`.
- `R/dkge-moments.R:145`: four q×q products per call; precompute `A = R Khalf`,
  use `crossprod(A, M %*% A)`. `.dkge_fold_weight_context` recomputes all
  subject moments from raw betas per fold (O(S²q²P)) and discards `contribs`.
- `.dkge_pool_effect_moments:186,191` recomputes `tcrossprod(mask)` per subject
  per bootstrap draw; sample-independent — cache.
- `R/dkge-fit-core.R:239-241, 650-661`: with `debias="none"`, `effect_moments`
  and `effect_moments_raw` are byte-identical and `noise_moments` all zero
  (~1 MB waste at q=40,S=30); `fit$subjects` now retains every design matrix.
- `R/dkge-between-infer.R:228`: `coef_rows` backsolve done per permutation even
  under `scope="global"`.
- Both aggregate resampling loops and the between-permutation loop are serial;
  the package uses `.dkge_lapply`/future.apply elsewhere.

## 5. Tests — why the suite is green

Genuinely strong: `test-data.R` "unobserved row values cannot leak" (1e12
poison vs hand value); `test-design-kernel.R` block-factor structure and grid
ordering; `test-effect-reliability.R:29-84` pairwise/precision/legacy oracles;
`test-trial-subject.R` QR, noise-subtraction, trace, cross-moment oracles.

Self-fulfilling / masking:

- `test-aggregate.R:190` hand-swaps `singular_values` alongside `U` — doing
  exactly the bookkeeping the code omits (B2). Remove that line; assert
  `singular_values == sqrt(colSums(scores_feature^2))`.
- `test-aggregate.R:149-151` restate the implementation; no rank>1 resampling
  test (hides B2); no singleton stratum (hides B1); "null calibration" at
  `:280` would pass with a no-op permutation.
- `test-fit.R:120-155` recomputes shrink with the buggy formula (B6);
  `:52-83` re-derives `Khalf[obs,obs]` (locks in the divergent convention);
  `:85-118` tests `.dkge_accumulate_chat`, which the package no longer calls.
- `test-plot-helpers.R`, `test-plot-suite.R:32-42`: `expect_s3_class` and row
  counts only; principal angles tested only with `K=I`, identical bases.
- `test-contrast.R:322-352` passes only because contrast names equal kernel
  term names.
- `test-effect-reliability.R:137-140` compares two implementations under
  review to each other (they diverge by construction under adaptive weights).

Missing golden tests that currently pass and should be pinned: new FL fast
path vs old refit; `.dkge_repool_fit` vs linear contrib path;
`.dkge_fold_weight_context(fit, 1:S)` == `fit$Chat`; `dkge_aggregate_fit` vs
direct SVD; `Chat == XXᵀ` under partial coverage.

## 6. Dead code, docs, hygiene

- Dead: `.dkge_accumulate_chat`, `.dkge_pair_counts_from_provenance`
  (`dkge-fit.R:216-311`); `.dkge_between_weight_response`,
  `.dkge_between_rrr_fast_eval` (`between-infer.R:400-445`);
  `Chat_sym` (`moments.R:325`); `folds.R:333`; `contrast.R:169 method`;
  `plot.R:667` recompute of `fold$basis_aligned`; `contrast_scale` back-door in
  aggregate `...`.
- R CMD check NOTEs: add `salience, score, contrast, projection, group` to
  `globals.R`; add `^AGENTS\.md$` to `.Rbuildignore`.
- Docs: `dkge_aggregate_bootstrap` inherits "B = number of permutation draws";
  stratified-by-default bootstrap undocumented; `@param K` on
  `dkge_aggregate_fit` says design_kernel() but q here = aggregate rows (a
  different space from the CLAUDE.md convention — say so); none of the 8 new
  plot exports nor `dkge_effect_grid` have `@examples`; new topics/vignette
  missing from `_pkgdown.yml`; fit-object field inventory at
  `dkge-fit-core.R:4-53` omits ~12 new fields; `Roxygen: markdown` not set so
  `[fn()]` links render literal.
- Repo: `testdata/` (126 MB) and `artifacts/` (79 MB) are untracked and **not
  gitignored** — one `git add .` away from a repo disaster; `inst/examples/`
  ships 2 MB of PDFs in the tarball; `.mote/`, `.omc/`, `.claude/` untracked.

## 7. Suggested order of work

1. B1–B3 + `rowsum` aggregation (one file, one PR); add rank>1 and singleton-
   stratum tests that fail first.
2. Decide §2 convention; then B4, B5, B6, empty-mask weight, delete dead
   accumulator, fix subject-weights, reword vignette; add `Chat == XXᵀ` and
   duplicate-ID tests.
3. B7–B9 + `unit_K` normalization + globals; add value-asserting plot tests.
4. Estimability structural classifier; `cbind` in validated summary.
5. Bootstrap linear-path predicate; K/effects dimname validation +
   `dkge_data(effects=)`.
6. `.gitignore` testdata/artifacts; `.Rbuildignore` AGENTS.md; pkgdown index;
   NEWS.

---

# Resolution — 2026-08-16

Every item above was addressed in a fix pass (three parallel slices by file
ownership, then an adversarial re-verification pass and a fresh-context code
review of the fixes, then a second fix round for what those found). All
CONFIRMED audit items were re-verified numerically after the fixes; the
default full-coverage path is unchanged vs `953d2a6` to machine precision
(`fit$K/Khalf/Btil/weights` exactly 0; `Chat` 4e-16 rel; LOSO contrast 4e-15).

## Blockers B1–B9 — all fixed and re-verified by an independent agent
B1 `sample.int` on the stratum indices + <2-subject warning · B2 rotated
`singular_values` from `scores_feature` · B3 signed statistic, `abs` only in
the two-sided p (new `alternative`), bootstrap `excludes_zero` · B4 duplicate
IDs rejected · B5 full-`Khalf` blocks (`Chat == XXᵀ` to 2e-16 rel) · B6
`diag(rescaled)` · B7 K-orthonormalize with rank detection before principal
angles · B8 name-preserving group matching · B9 signed projection
`⟨K u_j, rowMeans(B̃_s)⟩` (documented; flips under negation).

## Convention decision (§2) — implemented as recommended
Effect-space masking then full `Khalf` everywhere; subject weights switched to
it; empty-mask subject → weight 0 (+warning); `Khalf[obs,obs]` and the dead
accumulator deleted; vignette reworded (coupling K spreads observed energy by
design; `block_factors` for strict separation).

## Also fixed
Everything in §3–§6: kernel/effect dimname alignment + `dkge_data(effects=)`,
structural estimability (authoritative over contrast names), `cbind` in the
validated summary (+ positional indexing), `unit_K` normalization (default,
incl. identity fallback), design_kernel/effect_grid validation (terms, rho,
levels, unknown fields, `[[` accessors, `f$l`↔`f$levels` partial-match bug),
`x[-integer(0),]` trap, rescale under effect weights, dead code (accumulator,
pair-counts helper, between helpers, `Chat_sym`, …), duplicate
`.dkge_kernel_roots`, `rowsum` aggregation (bit-identical, ~100–400×), Gram
precompute `A = R Khalf`, fold-context moment reuse, pool cache, bootstrap
linear path via shared `.dkge_fit_pool_is_nonlinear`, `coef_rows` skip,
`parallel=` on aggregate/between resampling, args moved to end of signatures,
NEWS.md, pkgdown index, `.gitignore`/`.Rbuildignore`, `tests/testthat.R`
runner, Roxygen markdown mode.

## Found and fixed along the way (not in the original audit)
* `dkge_procrustes_K()` returned the TRANSPOSE of the optimal rotation
  (`U Vᵀ` instead of `V Uᵀ`); invisible at rank ≤ 2, permutes components at
  rank ≥ 3. Affects every K-Procrustes call site. Pre-existing in committed code.
* Default `w_method = "mfa_sigma1"` used a randomly seeded power iteration →
  fits (and component signs) depended on ambient RNG state. Now exact.
* `design_kernel()` read the length-scale as `f$l`, partial-matching `f$levels`.
* `pair_counts` became a sample-weighted mass but policies assumed integer
  counts (multiplier bootstrap); rescale/shrink under `effect_weights` were
  double-normalized; reduced-model refit in `dkge_between_permute()` ignored
  `tol`; `dkge_update_weights()` dropped the new fit settings on refit;
  `subjects[[s]]$id` never stamped; principal angles fabricated 0° on
  rank-deficient input; adaptive voxel weights silently recycled across
  heterogeneous `P_s` (now a clear error); several vacuous/self-referential
  tests replaced with oracle-based ones.

## Left as flagged (pre-existing, out of scope)
* `design_kernel()` with a single factor emits a doubled effect space with
  duplicated labels (default `terms = list("A","A")`); several tests depend on
  it. Kernel/effect name reconciliation warns and skips in that case.
* Freedman–Lane between-subject p-values are mildly anti-conservative in a
  300-sim null (type-I ≈ 0.07 at 0.05) — documented in `@details`; not a
  regression.

Final: `devtools::test()` 2261 pass / 0 fail; `R CMD check` 0 errors, 0 notes
(local clang/R-header install warning only).
