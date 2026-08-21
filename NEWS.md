# dkge 0.0.0.9000

## New features

* `dkge_signflip_maxT()` now explicitly exposes both max-T FWER-adjusted `p`
  and per-location unadjusted `p_unadj`. The latter is the raw value reported
  by `dkge_infer()` in `p_values`, while `p_adjusted` retains max-T control;
  result axes carry stable subject, feature, and permutation labels.
* The partial-effect-space article is restored and now executes coverage,
  pair-ESS, weighted chunked-debiasing, and estimability-warning contracts.
  The frozen 8,100-result between-subject rotation calibration was also rerun
  from exact clean snapshot `5b185e9f`; replicate and summary artifacts were
  byte-identical to the historical run, whose missing dirty-tree digest is now
  explicitly labeled non-certifying on its own.
* **Partial effect spaces.** Subjects may observe only a subset of the design
  effects (`observed_rows` on `dkge_subject()`; `dkge_effect_grid()` for a
  canonical cell grid). Coverage provenance (`obs_mask`, `pair_counts`) flows
  into `dkge_fit()`, which gains `missingness = c("none", "mask", "rescale",
  "shrink")` and `miss_args`. Masking acts in raw effect space *before* the
  `R`/`K^{1/2}` congruence; a coupling kernel therefore spreads observed energy
  into unobserved coordinates by design — use `block_factors` in
  `design_kernel()` for strict separation. See `?dkge_subject` and `?dkge_fit`.
* **Effect-reliability weighting and debiasing.** `dkge_effect_weights()`
  (`"none"`, `"count"`, `"precision"`) and `dkge_fit(debias =
  c("none", "analytic", "split_half"))` on top of a new central moment-pooling
  engine. `dkge_trial_subject()` and `dkge_trial_subject_chunks()` build
  subjects from trial-level designs with the sufficient statistics debiasing
  needs.
* **Aggregate (cell-mean) decomposition.** `dkge_aggregate_target()`,
  `dkge_aggregate_fit()`, `dkge_aggregate_align()`, `dkge_aggregate_stat()`,
  `dkge_aggregate_permute()`, `dkge_aggregate_bootstrap()`: a PLS-style
  decomposition of group-by-cell aggregate maps in the kernel metric with
  permutation and bootstrap inference (`alternative`, `parallel`).
* **Component diagnostics and plots.** `dkge_component_saliences()`,
  `dkge_component_contrast_scores()`, `dkge_design_basis()`,
  `dkge_subject_component_projections()`, `dkge_principal_angles_K()`, and
  `dkge_plot_*` counterparts.
* `dkge_data()`/`dkge()` gain `effects=` to pin the effect union order (e.g. to
  `dkge_effect_grid()$cell_labels`). `dkge_fit()` reorders a kernel whose
  dimnames are a permutation of the data's effects and errors on a set
  mismatch.
* `dkge_fit(effect_scaling = "none")` keeps effect rows on the input scale.
* `dkge_between_permute()` gains `method = "rotation"` for Haar rotation in
  the reduced model's residual space. It is finite-sample exact under the
  documented matrix-normal row-sphericity model, preserves `crossprod(Y)`, and
  fails closed for weights, multiple blocks, or fewer than two residual
  dimensions. It remains opt-in: a frozen 8,100-result audit fixed the earlier
  Freedman–Lane null inflation under a global Gaussian null, but missed one
  nuisance-null promotion gate and both predeclared power gates (Freedman–Lane's
  raw power advantage is mostly size inflation). The function
  also gains `parallel=`; both methods use compressed evaluation and give
  seeded serial/parallel-identical results.
* `dkge_subject_model()` gains `subject_id_col`.
* Contrast estimability (`between`/`within`/`mixed`) is now classified
  structurally from the kernel's factor scopes, so `dkge_contrast(method =
  "loso")` warns on between-subject contrasts regardless of how they are named.

## Bug fixes

* Sinkhorn transport now separates joint couplings from value-application
  operators. Intensive fields preserve constants, extensive values preserve
  total mass, and null support is represented by zero columns or rows in the
  corresponding application operator. Fitted reliability is not applied twice,
  and each solve reports convergence, iteration, marginal-error, and cache-hit
  diagnostics. Warm-start keys digest the complete numerical problem and
  non-converged states are not cached; the legacy `sinkhorn_cpp` method name is
  a deprecated alias because the main path already uses C++.
* K-Procrustes now reports the achieved proper-rotation objective when
  reflections are forbidden, validates PSD kernels and K-orthonormal inputs,
  accepts arbitrary eigenvector sign reflections in fold/analytic alignment,
  and validates consensus controls. `dkge_sim_toy()` also samples named term
  blocks by position, fixing the scalar-`sample()` ambiguity that could plant
  duplicate, metric-singular components.
* Design-kernel metadata now declares cell and effect coordinate spaces and
  names both axes of the cell-to-effect map. Kernel permutations align only the
  declared K axis; design-basis, target, plotting, classification, and
  estimability consumers rematch names and fail closed on ambiguous mappings.
  In particular, equal cell/effect dimensions are no longer treated as proof
  that cell metadata indexes effect rows.
* Analytic LOSO now reports structural fallback causes before numerical ones.
  Pair-normalized effect/missingness pooling is again labelled
  `pair_normalized_pooling`, covariance-aware moments are labelled
  `covariance_aware_moment`, and a later large perturbation cannot mask either
  primary cause.
* `dkge_procrustes_K()` returned the transpose of the optimal rotation; the
  error was invisible at rank ≤ 2 (transpositions are involutions) but permuted
  components at rank ≥ 3. All K-Procrustes call sites (bootstrap, analytic,
  folds, aggregate alignment, neuralign adapter) are affected.
* `design_kernel()`: the RBF length-scale was read as `f$l`, which partially
  matched `f$levels`, breaking ordinal/circular/continuous factors carrying
  level labels.
* Duplicate subject IDs are now rejected in `dkge_data()`.
* Duplicated effect labels are rejected on `dkge_subject()`, `dkge_data()`,
  and the union-alignment path instead of dropping a beta row.
* `dkge_aggregate_permute()` evaluates the statistic on the unaligned null
  refit (alignment is diagnostic only). `dkge_aggregate_bootstrap()` gains
  `interval = c("percentile", "basic")` and returns `excludes_zero = NA` for
  the non-negative `"singular_value"` statistic.
* `dkge_aggregate_target()` keeps an inferred cell factor when
  `cell_data = NULL`, matches rows/columns by name, and warns about unused
  `values` names.
* Chunked trial subjects no longer let an unweighted `noise_trace` short-circuit
  analytic debiasing when Omega or voxel weights are present.
* Rescale/shrink missingness divide by subject-weight mass (`pair_weight`)
  rather than unweighted pair counts; `none`/`mask` restore observed pair mass
  under partial coverage (no `S/k` inflation).
* `dkge_between_permute(terms = NULL)` skips terms listed in
  `design$nuisance`. Rotation accepts a factor with unused block levels.
* `dkge_subject_model(na.action = na.omit)` keeps subject IDs aligned with
  retained rows.
* Subject weights under the default `w_method = "mfa_sigma1"` retain the
  legacy power-iteration numerical contract inside a frozen private RNG scope.
  Fits are reproducible and leave the caller's RNG unchanged, while matching
  the canonical pre-extension default weights and downstream results.
* Contrast estimability now prefers structural evidence over a contrast's
  name, so a between-subject contrast named after a within-subject term is
  still flagged.
* `dkge_update_weights()` refits now preserve `effect_scaling`,
  `effect_weights`, `debias`, `missingness`, and `miss_args`.

## Breaking / behaviour changes

* A one-factor `design_kernel(terms = NULL)` now contains its main-effect term
  once. Previously the identical main effect and full interaction were both
  added, which doubled an unnormalised cell kernel and duplicated the
  effect-basis block and labels. Multi-factor defaults are unchanged; callers
  that intentionally need the old cell-kernel scale can set the sole term's
  `rho` to 2 explicitly.
* Cross-fitting helpers (`dkge_loso_contrast()`, `dkge_cv_*`, k-fold builders)
  now inherit `missingness` from the fit instead of defaulting to `"none"`.
  Numerically identical for full-coverage fits.
* `dkge_design_basis()` / `dkge_component_contrast_scores()` default to
  `normalize = "unit_K"` (unit K-norm) instead of unit Euclidean norm;
  `dkge_component_saliences(scale = "unit")` now means unit K-norm.
* `fit$subjects` no longer retains per-subject `beta`/`design`/`omega`
  (only debiasing sufficient statistics); `fit$effect_moments_raw` aliases
  `fit$effect_moments` and `fit$noise_moments` is `NULL` when
  `debias = "none"`.
* New arguments to `dkge_fit()`, `dkge()`, and `dkge_subject_model()` are
  appended after all pre-existing ones; positional calls from earlier
  releases are unaffected.
* `dkge_make_target(type = "transported_maps", centroids = NULL)` now keeps
  every supplied contrast (previously contrasts after the first were silently
  dropped) and labels features `<contrast>:<index>`, so `ncol(Y)` grows with
  the number of contrasts.
* `dkge_fit()`: a design kernel whose labels cannot be reconciled with the
  data's effect names (duplicated labels) now warns; a kernel with
  `rownames != colnames` is an error.
* Adaptive/prior voxel weighting (`dkge_weights()`) now errors when subjects
  have different numbers of voxels/clusters instead of silently recycling one
  subject's weights onto another.
* Roxygen markdown mode is enabled package-wide (documentation now renders
  cross-reference links and code spans correctly).
