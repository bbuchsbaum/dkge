# Changelog

## dkge 0.0.0.9000

### New features

- **Model-level spatial regularization.**
  [`dkge_spatial_regularizer()`](https://bbuchsbaum.github.io/dkge/reference/dkge_spatial_regularizer.md)
  builds graph Laplacians with `adjoin` (or accepts precomputed
  Laplacians), binds them to shared or subject-specific beta-column
  domains, and applies sparse `(I + lambda * L)^{-1}` solves inside the
  pooled moment and all reconstructed component, contrast, bootstrap,
  transport, and prediction fields. The regularizer can therefore change
  the learned q-space solution rather than merely blur its display.
  `lambda = 0` is exactly the unsmoothed fit;
  [`dkge_cv_spatial_grid()`](https://bbuchsbaum.github.io/dkge/reference/dkge_cv_spatial_grid.md)
  scores candidate penalties against raw held-out fields in a fixed
  validation geometry and uses the smoothest one-SE choice. Analytic
  noise debiasing fails closed because its diagonal residual-variance
  contract does not identify smoothing-induced spatial covariance;
  split-half debiasing remains supported.
- [`dkge_contrast_diagnostics()`](https://bbuchsbaum.github.io/dkge/reference/dkge_contrast_diagnostics.md)
  preflights planned contrasts against the fitted kernel support,
  reporting retained/null fractions and distinct contrasts that produce
  nearly proportional kernel queries. Numerical support tolerance is now
  separate from practical query collinearity; the returned summary names
  the maximally correlated query pair and reports the query-norm range.
- [`dkge_component_contrasts()`](https://bbuchsbaum.github.io/dkge/reference/dkge_component_contrasts.md)
  constructs the component-isolating contrast matrix `R %*% U`. This
  makes the contrast path distinct from the dual salience/read-out basis
  `K %*% U`, which must not be fed back as a contrast.
- [`dkge_signflip_maxT()`](https://bbuchsbaum.github.io/dkge/reference/dkge_signflip_maxT.md)
  now explicitly exposes both max-T FWER-adjusted `p` and per-location
  unadjusted `p_unadj`. The latter is the raw value reported by
  [`dkge_infer()`](https://bbuchsbaum.github.io/dkge/reference/dkge_infer.md)
  in `p_values`, while `p_adjusted` retains max-T control; result axes
  carry stable subject, feature, and permutation labels.
- The partial-effect-space article is restored and now executes
  coverage, pair-ESS, weighted chunked-debiasing, and
  estimability-warning contracts.
- **Partial effect spaces.** Subjects may observe only a subset of the
  design effects (`observed_rows` on
  [`dkge_subject()`](https://bbuchsbaum.github.io/dkge/reference/dkge_subject.md);
  [`dkge_effect_grid()`](https://bbuchsbaum.github.io/dkge/reference/dkge_effect_grid.md)
  for a canonical cell grid). Coverage provenance (`obs_mask`,
  `pair_counts`) flows into
  [`dkge_fit()`](https://bbuchsbaum.github.io/dkge/reference/dkge_fit.md),
  which gains `missingness = c("none", "mask", "rescale", "shrink")` and
  `miss_args`. Masking acts in raw effect space *before* the
  `R`/`K^{1/2}` congruence; a coupling kernel therefore spreads observed
  energy into unobserved coordinates by design — use `block_factors` in
  [`design_kernel()`](https://bbuchsbaum.github.io/dkge/reference/design_kernel.md)
  for strict separation. See
  [`?dkge_subject`](https://bbuchsbaum.github.io/dkge/reference/dkge_subject.md)
  and
  [`?dkge_fit`](https://bbuchsbaum.github.io/dkge/reference/dkge_fit.md).
- **Effect-reliability weighting and debiasing.**
  [`dkge_effect_weights()`](https://bbuchsbaum.github.io/dkge/reference/dkge_effect_weights.md)
  (`"none"`, `"count"`, `"precision"`) and
  `dkge_fit(debias = c("none", "analytic", "split_half"))` on top of a
  new central moment-pooling engine.
  [`dkge_trial_subject()`](https://bbuchsbaum.github.io/dkge/reference/dkge_trial_subject.md)
  and
  [`dkge_trial_subject_chunks()`](https://bbuchsbaum.github.io/dkge/reference/dkge_trial_subject_chunks.md)
  build subjects from trial-level designs with the sufficient statistics
  debiasing needs.
- **Aggregate (cell-mean) decomposition.**
  [`dkge_aggregate_target()`](https://bbuchsbaum.github.io/dkge/reference/dkge_aggregate_target.md),
  [`dkge_aggregate_fit()`](https://bbuchsbaum.github.io/dkge/reference/dkge_aggregate_fit.md),
  [`dkge_aggregate_align()`](https://bbuchsbaum.github.io/dkge/reference/dkge_aggregate_align.md),
  [`dkge_aggregate_stat()`](https://bbuchsbaum.github.io/dkge/reference/dkge_aggregate_stat.md),
  [`dkge_aggregate_permute()`](https://bbuchsbaum.github.io/dkge/reference/dkge_aggregate_permute.md),
  [`dkge_aggregate_bootstrap()`](https://bbuchsbaum.github.io/dkge/reference/dkge_aggregate_bootstrap.md):
  a PLS-style decomposition of group-by-cell aggregate maps in the
  kernel metric with permutation and bootstrap inference (`alternative`,
  `parallel`).
- **Component diagnostics and plots.**
  [`dkge_component_saliences()`](https://bbuchsbaum.github.io/dkge/reference/dkge_component_saliences.md),
  [`dkge_component_contrast_scores()`](https://bbuchsbaum.github.io/dkge/reference/dkge_component_contrast_scores.md),
  [`dkge_design_basis()`](https://bbuchsbaum.github.io/dkge/reference/dkge_design_basis.md),
  [`dkge_subject_component_projections()`](https://bbuchsbaum.github.io/dkge/reference/dkge_subject_component_projections.md),
  [`dkge_principal_angles_K()`](https://bbuchsbaum.github.io/dkge/reference/dkge_principal_angles_K.md),
  and `dkge_plot_*` counterparts.
- [`dkge_data()`](https://bbuchsbaum.github.io/dkge/reference/dkge_data.md)/[`dkge()`](https://bbuchsbaum.github.io/dkge/reference/dkge.md)
  gain `effects=` to pin the effect union order (e.g. to
  `dkge_effect_grid()$cell_labels`).
  [`dkge_fit()`](https://bbuchsbaum.github.io/dkge/reference/dkge_fit.md)
  reorders a kernel whose dimnames are a permutation of the data’s
  effects and errors on a set mismatch.
- `dkge_fit(effect_scaling = "none")` keeps effect rows on the input
  scale.
- [`dkge_between_permute()`](https://bbuchsbaum.github.io/dkge/reference/dkge_between_permute.md)
  gains `method = "rotation"` for Haar rotation in the reduced model’s
  residual space. It is finite-sample exact under the documented
  matrix-normal row-sphericity model, preserves `crossprod(Y)`, and
  fails closed for weights, multiple blocks, or fewer than two residual
  dimensions. It remains opt-in: a frozen 8,100-result audit fixed the
  earlier Freedman–Lane null inflation under a global Gaussian null, but
  missed one nuisance-null promotion gate and both predeclared power
  gates (Freedman–Lane’s raw power advantage is mostly size inflation).
  The function also gains `parallel=`; both methods use compressed
  evaluation and give seeded serial/parallel-identical results.
- [`dkge_subject_model()`](https://bbuchsbaum.github.io/dkge/reference/dkge_subject_model.md)
  gains `subject_id_col`.
- Contrast estimability (`between`/`within`/`mixed`) is now classified
  structurally from the kernel’s factor scopes, so
  `dkge_contrast(method = "loso")` warns on between-subject contrasts
  regardless of how they are named.

### Bug fixes

- Numerical-rank decisions now use one scale-equivariant
  transformed-moment tolerance across full fits, folds, LOSO, CV, and
  aggregate fits. Zero-signal aggregates report rank zero, and
  rank-deficient resamples contribute exact zeros for unavailable
  components instead of aborting; ridge remains inside the kernel range;
  fold weight dimensions and custom assessment partitions fail closed
  instead of recycling or silently duplicating subjects.
- Scalar permutation and transport controls are validated before
  coercion. Sinkhorn solves only on positive mass support, rejects
  malformed native warm starts, and re-expands structural zero rows and
  columns. Classification permutations require an independently fixed
  penalty and full representation recomputation for cell-based modes.
- Kernel roots now preserve exact positive-semidefinite support instead
  of jittering `null(K)` into artificial inverse directions. Fits cap
  latent rank at `rank(K)`, keep ridge inside `image(K)`, expose
  rank/nullity/condition diagnostics, and reject contrasts that the
  chosen kernel cannot represent.
- Kernel cross-validation now requires full-rank candidates by default
  and scores every candidate in one fixed validation geometry
  (identity/effect space by default). Intentional quotient models
  require `kernel_rank_policy = "allow_singular"`; excluded candidates
  remain in audit rows, and score saturation is flagged, preventing a
  low-rank kernel from grading itself only on the directions it kept.
- Full-rank but spectrally concentrated kernels now expose
  participation-ratio effective rank, effective-rank fraction, and
  leading-eigenvalue share in fit, contrast, and CV diagnostics. CV
  warns—without changing its predictive selection—when the selected
  kernel is effectively too concentrated for the selected latent rank.
- Sinkhorn transport now separates joint couplings from
  value-application operators. Intensive fields preserve constants,
  extensive values preserve total mass, fitted reliability is not
  applied twice, and each solve reports convergence, iteration,
  marginal-error, and cache-hit diagnostics. Warm-start keys digest the
  complete numerical problem and non-converged states are not cached;
  the legacy `sinkhorn_cpp` method name is a deprecated alias because
  the main path already uses C++.
- K-Procrustes now reports the achieved proper-rotation objective when
  reflections are forbidden, validates PSD kernels and K-orthonormal
  inputs, accepts arbitrary eigenvector sign reflections in
  fold/analytic alignment, and validates consensus controls.
  [`dkge_sim_toy()`](https://bbuchsbaum.github.io/dkge/reference/dkge_sim_toy.md)
  also samples named term blocks by position, fixing the
  scalar-[`sample()`](https://rdrr.io/r/base/sample.html) ambiguity that
  could plant duplicate, metric-singular components.
- Design-kernel metadata now declares cell and effect coordinate spaces
  and names both axes of the cell-to-effect map. Kernel permutations
  align only the declared K axis; design-basis, target, plotting,
  classification, and estimability consumers rematch names and fail
  closed on ambiguous mappings. In particular, equal cell/effect
  dimensions are no longer treated as proof that cell metadata indexes
  effect rows.
- Analytic LOSO now reports structural fallback causes before numerical
  ones. Pair-normalized effect/missingness pooling is again labeled
  `pair_normalized_pooling`, covariance-aware moments are labeled
  `covariance_aware_moment`, and a later large perturbation cannot mask
  either primary cause.
- [`dkge_procrustes_K()`](https://bbuchsbaum.github.io/dkge/reference/dkge_procrustes_K.md)
  returned the transpose of the optimal rotation; the error was
  invisible at rank ≤ 2 (transpositions are involutions) but permuted
  components at rank ≥ 3. All K-Procrustes call sites (bootstrap,
  analytic, folds, aggregate alignment, neuralign adapter) are affected.
- [`design_kernel()`](https://bbuchsbaum.github.io/dkge/reference/design_kernel.md):
  the RBF length-scale was read as `f$l`, which partially matched
  `f$levels`, breaking ordinal/circular/continuous factors carrying
  level labels.
- [`dkge_data()`](https://bbuchsbaum.github.io/dkge/reference/dkge_data.md)
  now rejects duplicate subject IDs.
- [`dkge_subject()`](https://bbuchsbaum.github.io/dkge/reference/dkge_subject.md),
  [`dkge_data()`](https://bbuchsbaum.github.io/dkge/reference/dkge_data.md),
  and the union-alignment path now reject duplicated effect labels
  instead of dropping a beta row.
- [`dkge_aggregate_permute()`](https://bbuchsbaum.github.io/dkge/reference/dkge_aggregate_permute.md)
  evaluates the statistic on the unaligned null refit (alignment is
  diagnostic only).
  [`dkge_aggregate_bootstrap()`](https://bbuchsbaum.github.io/dkge/reference/dkge_aggregate_bootstrap.md)
  gains `interval = c("percentile", "basic")` and returns
  `excludes_zero = NA` for the non-negative `"singular_value"`
  statistic.
- [`dkge_aggregate_target()`](https://bbuchsbaum.github.io/dkge/reference/dkge_aggregate_target.md)
  keeps an inferred cell factor when `cell_data = NULL`, matches
  rows/columns by name, and warns about unused `values` names.
- Chunked trial subjects no longer let an unweighted `noise_trace`
  short-circuit analytic debiasing when Omega or voxel weights are
  present.
- Rescale/shrink missingness divide by subject-weight mass
  (`pair_weight`) rather than unweighted pair counts; `none`/`mask`
  restore observed pair mass under partial coverage (no `S/k`
  inflation).
- `dkge_between_permute(terms = NULL)` skips terms listed in
  `design$nuisance`. Rotation accepts a factor with unused block levels.
- `dkge_subject_model(na.action = na.omit)` keeps subject IDs aligned
  with retained rows.
- Subject weights under the default `w_method = "mfa_sigma1"` retain the
  legacy power-iteration numerical contract inside a frozen private RNG
  scope. Fits are reproducible and leave the caller’s RNG unchanged,
  while matching the canonical pre-extension default weights and
  downstream results.
- Contrast estimability now prefers structural evidence over a
  contrast’s name, so a between-subject contrast named after a
  within-subject term is still flagged.
- [`dkge_update_weights()`](https://bbuchsbaum.github.io/dkge/reference/dkge_update_weights.md)
  refits now preserve `effect_scaling`, `effect_weights`, `debias`,
  `missingness`, and `miss_args`.

### Breaking / behavior changes

- [`kernel_roots()`](https://bbuchsbaum.github.io/dkge/reference/kernel_roots.md)
  now defaults to `jitter = 0` and returns a Moore–Penrose inverse
  square root for singular kernels. A positive `jitter` explicitly
  regularizes `K`; it no longer happens invisibly on every call.
  `dkge_cv_*()` scores are therefore comparable across kernels but are
  not numerically interchangeable with the former candidate-specific
  K-metric scores.
- `design_kernel(terms = NULL)` on a one-factor design now contains its
  main-effect term once. Previously the identical main effect and full
  interaction were both added, which doubled an unnormalized cell kernel
  and duplicated the effect-basis block and labels. Multi-factor
  defaults are unchanged; callers that intentionally need the old
  cell-kernel scale can set the sole term’s `rho` to 2 explicitly.
- Cross-fitting helpers
  ([`dkge_loso_contrast()`](https://bbuchsbaum.github.io/dkge/reference/dkge_loso_contrast.md),
  `dkge_cv_*`, k-fold builders) now inherit `missingness` from the fit
  instead of defaulting to `"none"`. Numerically identical for
  full-coverage fits.
- [`dkge_design_basis()`](https://bbuchsbaum.github.io/dkge/reference/dkge_design_basis.md)
  /
  [`dkge_component_contrast_scores()`](https://bbuchsbaum.github.io/dkge/reference/dkge_component_contrast_scores.md)
  default to `normalize = "unit_K"` (unit K-norm) instead of unit
  Euclidean norm; `dkge_component_saliences(scale = "unit")` now means
  unit K-norm.
- `fit$subjects` no longer retains per-subject `beta`/`design`/`omega`
  (only debiasing sufficient statistics); `fit$effect_moments_raw`
  aliases `fit$effect_moments` and `fit$noise_moments` is `NULL` when
  `debias = "none"`.
- New arguments to
  [`dkge_fit()`](https://bbuchsbaum.github.io/dkge/reference/dkge_fit.md),
  [`dkge()`](https://bbuchsbaum.github.io/dkge/reference/dkge.md), and
  [`dkge_subject_model()`](https://bbuchsbaum.github.io/dkge/reference/dkge_subject_model.md)
  are appended after all pre-existing ones; positional calls from
  earlier releases are unaffected.
- `dkge_make_target(type = "transported_maps", centroids = NULL)` now
  keeps every supplied contrast (previously contrasts after the first
  were silently dropped) and labels features `<contrast>:<index>`, so
  `ncol(Y)` grows with the number of contrasts.
- [`dkge_fit()`](https://bbuchsbaum.github.io/dkge/reference/dkge_fit.md):
  a design kernel whose labels cannot be reconciled with the data’s
  effect names (duplicated labels) now warns; a kernel with
  `rownames != colnames` is an error.
- [`dkge_weights()`](https://bbuchsbaum.github.io/dkge/reference/dkge_weights.md)
  now errors when subjects have different numbers of voxels/clusters
  instead of silently recycling one subject’s weights onto another.
- Roxygen markdown mode is enabled package-wide (documentation now
  renders cross-reference links and code spans correctly).
