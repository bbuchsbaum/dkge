# Review of new paths (code + docs) — 2026-08-17

Scope: everything in the working tree since `953d2a6` after the 2026-08-16 fix
pass — partial effect spaces / effect reliability / trialwise, aggregate
decomposition, between-subject rotation, component diagnostics/plots,
contrast estimability, and the seven touched vignettes. Five independent
fresh-context reviewers took disjoint slices against a scratch install of the
working tree; every CONFIRMED item was reproduced numerically, and the
top-severity items were re-run by hand (scripts under the session scratchpad
`rev-*/`).

Baseline: focused test files green; `R CMD check` per author: 0 E / 1 W
(toolchain) / 2,178 tests.

## 0. Verdict

The 2026-08-16 fix pass held: **every** prior audit blocker (B1–B9), the
convention decision, and the "also fixed" list were re-verified OK, including
the `dkge_procrustes_K` transpose fix (all call sites consistent), the
deterministic `mfa_sigma1` weight, and the FL fast path vs refit oracle. The
rotation permutation is correctly implemented (proper Haar, reduced-model
residual space, `crossprod(Y)` preserved, exact vs the F oracle at P=1, at
nominal size in four independent null arms incl. non-Gaussian). The frozen
8,100-replicate audit reproduces bit-exactly and all quoted numbers match.

New problems are again at the seams: **three silent-wrong-number defects
reachable from defaults or documented entry points**, one invalid null, and a
handful of vignette statements that the code contradicts.

## 1. Blockers (wrong numbers, reachable via documented paths)

| # | Where | What | Fix |
|---|-------|------|-----|
| N1 | `R/dkge-aggregate.R:1043` | `dkge_aggregate_permute()` K-Procrustes-aligns **every null draw to the observed fit** before taking the statistic; observed is evaluated unrotated. Aligned `sv₁ = √Σ d_k² R²_k1 ≤ d₁`, so nulls only shrink. Type-I at α=.05, rank=3, null data: `singular_value` **0.275**, `component_score_contrast` **0.135**, `salience_contrast` **0.000**; unaligned 0.06/0.10/0.07; rank=1 exact (0.050). Repro: p_aligned 0.030 vs p_unaligned 0.194 on pure-null data. | Do not align inside the permutation loop; two-sided `abs` handles sign. Keep alignment in the bootstrap (observed fit is a legitimate reference). Keep `alignment_summary` as diagnostic only. |
| N2 | `R/dkge-fit-core.R:147-152` vs `:201-203, :784` | Kernel is permuted to match data effect order but **`kernel_info` (cells, cell_labels, blocks) is not**; `dkge_design_basis()` (`R/dkge-plot.R:292-300`) and `.dkge_contrast_structural_scope()` (`R/dkge-contrast.R:243-253`) index `cells` positionally. Verified: `task` column of design basis wrong for 3/6 rows; `dkge_component_contrast_scores` reports −0.148 where oracle is 1.963; between contrast classified `"within"` (LOSO warning does not fire) and within classified `"between"` (spurious warning). Reachable whenever subject beta rownames differ from kernel cell order — the normal case under partial coverage / `dkge_data()` union order. | Permute `kernel_info$cells/cell_labels` and remap `$blocks` in `.dkge_align_kernel_effects()`, **and/or** match by name in both consumers (`idx <- match(fit$effects, info$cell_labels)`) with fallback to identity basis / NULL scope. Add a permuted-rowname fixture test. |
| N3 | `R/dkge-trial-subject.R:255-259, 284`; `R/dkge-moments.R:152` | `dkge_trial_subject_chunks()` stores an **unweighted** `noise_trace`; `.dkge_noise_trace()` short-circuits on it and drops `Omega`/voxel weights. `noise_trace_scope="unweighted"` is written, never read. `debias="analytic"`: dense vs chunked `Chat` differ (rel 3.3e-3 with `omega` spread 0.2–3; grows with spread). Contradicts `?dkge_trial_subject_chunks` and vignette :316. | Drop `noise_trace` from the chunked constructor (recomputable from `residual_variance`), or honour `noise_trace_scope`. Add a chunk test with non-unit `omega`. |
| N4 | `R/dkge-moments.R:360` → `R/dkge-fit.R:246,253` | Unweighted branch of `missingness="rescale"/"shrink"` divides `Σ_s w_s M_s` by **unweighted** `pair_counts` (effect-weighted branch correctly uses `pair_weight`, `:337`). Under default `w_method="mfa_sigma1"` per-pair scale depends on which subjects observed the pair — the thing rescale exists to remove (entry ratios 0.57–1.34 in a 3-subject repro). Newly reachable from `dkge_fit()`. | Divide by `pair_weight` (already computed at `R/dkge-moments.R:323`) in both branches. |
| N5 | `R/dkge-data.R:626`; `R/dkge-align-data.R:30,38-42` | Duplicated effect labels (e.g. `design_kernel()` single-factor default `terms=list("A","A")`, or repeated colnames) are not rejected: `match()` duplicates metadata rows (`effect_n`, `effect_noise_cov`, precision) while betas stay put; in the union path one real beta row is **lost** and another zero-filled with `obs_mask` still TRUE. | Reject `anyDuplicated(effects)` in `.align_effects()`/`dkge_data()`; error in `.dkge_align_subjects_to_union()`. Closes the open single-factor hazard from memory. |
| N6 | `R/dkge-aggregate.R:116,135` / `:128` / `:176,183` | `dkge_aggregate_target()`: (a) `cell_data=NULL` builds then discards the inferred cell factor → all cells collapse to group means (doc says otherwise); (b) subject rows matched positionally when `cell_data` given without `cell_id_col`, ignoring rownames; (c) feature columns of subject 1 are stamped onto every other subject — a column-permuted subject is silently mislabeled. | (a) default `cell_vars` to inferred cell; (b) match rownames whenever present; (c) require/reorder by colnames, error otherwise. |

## 2. Significant (contract or reproducibility violations, narrower reach)

- **Serial ≠ parallel for ≥2 terms** — `R/dkge-between-infer.R:336-341`: `future_lapply(future.seed=TRUE)` advances `.Random.seed`; descriptors drawn inside per-term `lapply`, so term 2+ differ (both methods; p 0.115 vs 0.154 in repro). Falsifies roxygen :27-31 and NEWS. Test covers one term only. Fix: draw all descriptors before the term loop, or `.dkge_seed_enter/_exit` around `.dkge_apply()`.
- **Power verdict is raw-size** — rotation 0.430 vs FL 0.580 compares an exact test against one with size ≈0.086 at nominal .05. At size-matched threshold FL ≈ 0.487; predeclared *relative* gate would pass (0.430 ≥ 0.367). Keep the predeclared gate outcome; add the size-adjusted comparison in report/NEWS/roxygen :70-74/vignette and drop "expect a possible power cost".
- `missingness="none"/"mask"` + `effect_weights` restores cohort mass (`R/dkge-moments.R:338`) → pairs observed by k of S subjects inflated by S/k (mean imputation of absent subjects; verified ×6). `?dkge_fit` :355 says all policies coincide with unweighted under unit precision — false for these two. Reword and/or gate restoration on full coverage.
- Aggregate bootstrap percentile interval for `singular_value` routinely excludes the observed value (upward resampling bias; 99.3% of draws above observed) and `excludes_zero` is vacuous for non-negative statistics — document or offer BC/BCa. `degeneracy_tol=1e-6` only catches exact ties.
- `dkge_design_basis()` dies on non-syntactic factor names (`~ my factor`) and on 1-level factors (`R/dkge-plot.R:316-326`).
- `dkge_preprocess_blocks()` (`R/dkge-project.R:66-67`) stamps *training* parcel names onto new data. `dkge_transform_block()` ignores `fit$voxel_weights_subject`, contradicting new comment at `R/dkge-fit-core.R:706-715` (rel 0.14 mismatch with a voxel prior).
- `miss_args` names never validated (`gama=3` silently ignored) — same class as the fixed `f$l` bug.
- Rotation single-block guard uses `nlevels()` → trips on unused factor level (`droplevels`). `dkge_subject_model(na.action=na.omit)` errors misleadingly; `nuisance=` stored but never consumed; energy guard scaled by QR `tol` (`R/dkge-between-infer.R:471`).
- `contrast_estimability` junk row names; `dkge_plot_subject_contrib()` panels order subjects differently; block separator lines assume index-ordered blocks; user `term` attr dropped in `.dkge_validate_design_basis()`; `basis$basis` partial-matching.

## 3. Vignette / doc factual errors (ranked)

Blocking (reader would act on a false statement):
1. `vignettes/dkge-workflow.Rmd:150` — `boot_q$summary[[1]]$medoid$sd[1]` is `NULL` (no voxel operator → spliced top level). Use `$sd[1]` (0.2156).
2. `vignettes/dkge-concepts.Rmd:149-152, 254-255` — `crossprod(fit_I$U, fit_K$U)` can never be near-identity (K-orthonormal columns; norms 1.5–2.0). Use K-metric cross-product or `dkge_principal_angles_K()`.
3. `vignettes/dkge-architecture.Rmd:104-107` — cohort-mass exception list omits `"shrink"` (`R/dkge-moments.R:357`).
4. `vignettes/dkge-weighting.Rmd:163-193` — both fits use `w_method="none"` so printed weights are `1 1 1`; prose claims reduced leverage. Use `mfa_sigma1` or print `sdev`.
5. `vignettes/dkge.Rmd:46` — no `print.dkge`; multivarious printer shows neither eigenspectrum nor subject count.
6. `vignettes/dkge-concepts.Rmd:183` (and runtime estimability warning, partial-effect :370) — "subject-label permutation / bootstrap of U" name no exported callable for a `dkge_fit`.

Non-blocking: workflow :116-117 `comp$transport` shape wrong **and** `R/dkge-components.R:80` indexes medoid clusters by component index (drops clusters 3–4 at rank 2 — pre-existing bug, also affects weighting :255-273); partial-effect :157-158 subject weights are computed after `R`/`Khalf` mixing, not on zero rows; :134 `mask` thresholds `pair_ess` under effect weights; :316 chunked noise trace claim (N3); between-subjects :94-97 drops rank-clipping caveat; :150-151 NA blocks not rotation-specific; :195-217 two frozen studies juxtaposed without stating shared DGP, power arm B=199 not 399; weighting :63 "manual overrides" name nothing; concepts :84 "cell-mean PCA" is uncentered by default; `w_method="none"` used throughout partial-effect vignette without saying it departs from the default; `dkge.Rmd:92` signature `dkge(subjects, K, rank)` misorders positionals.

Roxygen: `R/dkge-plot.R:136-143` cites `dkge_k_orthonormalize()` but code uses `svd(Khalf U)`; :818-826 projection/energy are not a signed/unsigned pair (voxel weights differ); `R/dkge-contrast.R:143-152` name-based downstream matching does not exist (implementing it is the fix for N2); `.dkge_apply_missingness()` docs still say "compressed covariance"; `.dkge_subject_weights()` overclaims "same quantity the eigensolve sees"; between-infer :44-53 "requested rank" and clipping caveat over-applied to rotation; `dkge_subject_model` `@param subject_ids`/`nuisance`/`na.action`; no `@examples` on any of the six `dkge_aggregate_*` exports nor `dkge_between_permute/_rrr`, `dkge_subject_model`.

## 4. Test-quality notes

- Green because untested: `test-aggregate.R:807` calibration runs `rank=1` (only exact case); `:651` asserts `min|boot stat| > 0.25·observed`, which *rewards* the N1 alignment attraction; every `kernel_info` fixture uses kernel-order rownames (N2 invisible); `test-trial-subject-chunks.R:22` pins the unweighted trace (N3); `test-effect-reliability.R:228-248` full coverage only (`none/mask` inflation invisible); `:137-140` context vs repool is self-calling; parallel-identity test single-term.
- Genuine oracles worth keeping: procrustes 3-cycle + trace test; FL vs dense refit; rotation compressed vs dense; frozen replicate rows reproduce bit-exactly; `Chat == XXᵀ` under partial coverage; hand `⟨Ku, rowMeans⟩` projection.
- Add: permutation type-I at rank ≥ 2; permuted-rowname design-basis/estimability; chunk with `omega`; rescale with non-unit subject weights under partial coverage; duplicated effect labels; ≥2-term serial/parallel identity; rownames/colnames mismatch in aggregate targets.

## 5. Hygiene

- `git status`: `man/dkge_between_permute.Rd`, `dkge_between_rrr.Rd`, `dkge_make_target.Rd`, `dkge_subject_model.Rd`, `dkge_term_map.Rd` are **staged as deleted** while regenerated copies sit untracked — a plain commit drops them. `git add man/` before committing.
- `_pkgdown.yml`: `dkge_trial_subject_chunks` is the only export not indexed (build_reference errors); also missing from NEWS.
- `dev/` not in `.Rbuildignore` (CRAN NOTE); `inst/extdata` ships 1.5 MB of raw p-value CSVs; `.claude/` untracked & un-ignored.
- Evidence reports record `HEAD 953d2a6` but ran uncommitted code; record a working-tree digest.

## 6. Suggested order

1. N1 (one function, one test that fails first) · N2 (align `kernel_info` + name matching, permuted fixture) · N3 · N4 · N5 · N6.
2. Parallel seed fix + 2-term test; `miss_args` validation; design-basis robustness; `dkge_preprocess_blocks` names.
3. Vignette blockers 1–6 and the roxygen contradictions; size-adjusted power paragraph.
4. Hygiene: `git add man/`, pkgdown/NEWS entry, `.Rbuildignore dev`.
