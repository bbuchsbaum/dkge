# DKGE Concepts: What Is Estimated and What Can Be Claimed?

After fitting DKGE, you may see a stable-looking component, a strong
planned contrast, or a localized cluster map. Those are different
results. Before interpreting any of them, you need to know which
quantity DKGE decomposed, how the design kernel changed its geometry,
and which level of evidence your next procedure can support.

This page supplies that mental model. For runnable setup and spatial
transport, start with
[`vignette("dkge")`](https://bbuchsbaum.github.io/dkge/articles/dkge.md)
and
[`vignette("dkge-workflow")`](https://bbuchsbaum.github.io/dkge/articles/dkge-workflow.md).

## What enters the decomposition?

For subject $`s`$, let $`B_s`$ be the `q` by `P_s` beta matrix. In the
simplest unweighted case, its raw effect-space moment is

``` math
M_s = B_s B_s^\top.
```

Spatial weights replace this with $`B_s \Omega_s B_s^\top`$. DKGE first
builds these small `q` by `q` moments, then pools them across subjects.
This is why the core fit scales with the number of design effects rather
than with a feature-by-feature covariance matrix.

The pooled raw moment $`M`$ is transformed as

``` math
\widehat C = K^{1/2} R^\top M R K^{1/2},
```

where $`R`$ is the pooled design ruler and $`K`$ is the design kernel.
DKGE eigendecomposes $`\widehat C`$, then maps its retained eigenvectors
back through $`K^{-1/2}`$ to obtain the group basis $`U`$. The columns
of $`U`$ are K-orthonormal: $`U^\top K U = I`$.

You can inspect each stage on a fitted object:

``` r

toy <- dkge_sim_toy(
  factors = list(condition = list(L = 2), load = list(L = 3)),
  active_terms = c("condition", "load"),
  S = 5, P = 20, snr = 5, seed = 2024
)
fit <- dkge(toy$B_list, toy$X_list, K = toy$K, rank = 2)

c(
  effects = nrow(fit$effect_moment),
  transformed_rows = nrow(fit$Chat),
  retained_components = ncol(fit$U)
)
#>             effects    transformed_rows retained_components 
#>                   5                   5                   2
```

**Input:** subject effect moments.

**Output:** a low-rank basis in the kernel metric.

**Next operation:** inspect $`K U`$ with
[`dkge_component_saliences()`](https://bbuchsbaum.github.io/dkge/reference/dkge_component_saliences.md),
or project a prespecified contrast with
[`dkge_contrast()`](https://bbuchsbaum.github.io/dkge/reference/dkge_contrast.md).

## What does the kernel change?

The kernel is a metric, not a fixed set of component directions. It
changes which effect-space directions count as large or smooth during
the eigensolve; the empirical moment still determines the fitted
components.

An identity kernel is the essential baseline. It applies no
kernel-imposed coupling among effects, although the full subject-level
fit still uses GLM effect rows and, by default, the pooled design ruler
$`R`$. A structured kernel regularizes the fit toward directions favored
by its factor, ordinal, or block structure.

Fit both versions on the same data:

``` r

fit_identity <- dkge(
  toy$B_list, toy$X_list, K = diag(nrow(toy$K)), rank = 2,
  w_method = "none"
)
fit_structured <- dkge(
  toy$B_list, toy$X_list, K = toy$K, rank = 2,
  w_method = "none"
)

data.frame(
  component = 1:2,
  identity = dkge_variance_explained(fit_identity)$prop_var,
  structured = dkge_variance_explained(fit_structured)$prop_var
)
#>   component identity structured
#> 1         1    0.615      0.568
#> 2         2    0.385      0.432
```

Compare retained subspaces in the structured metric:

``` r

angles <- dkge_principal_angles_K(
  fit_identity$U, fit_structured$U, fit_structured$K
)
round(angles * 180 / pi, 1)
#> [1] 0 0
```

Small angles mean the two fits retained similar subspaces under that
metric; large angles mean the structured kernel materially changed them.
A raw `crossprod(fit_identity$U, fit_structured$U)` is not a valid
substitute because the bases are not Euclidean-orthonormal in the same
geometry.

Both fits disable subject block weighting, isolating the kernel change
in this comparison. Neither outcome proves the structured kernel is
correct. Report the kernel as a modeling choice and the identity
comparison as a sensitivity analysis.

## What happens before the kernel when effects are incomplete?

Coverage, effect precision, and finite-trial debiasing alter the raw
pooled moment $`M`$ before the $`R`$ and $`K`$ transforms. They are not
cosmetic weights on already fitted components.

This ordering matters:

1.  mark which effect rows each subject actually observed;
2.  construct subject moments, optionally subtracting analytic noise or
    using a split-half cross-moment;
3.  pool observed pairs with the selected effect-precision and
    missingness policy; and
4.  apply $`R`$, $`K`$, and the eigensolve.

An absent cell is therefore not an observed zero.
[`dkge_data()`](https://bbuchsbaum.github.io/dkge/reference/dkge_data.md)
zero-fills internal aligned matrices only after recording observation
masks, and the fit uses those masks when pooling. Policies such as
`missingness = "rescale"` or `"shrink"` change the estimand; they must
be chosen and reported, not used as silent repairs.

See
[`vignette("dkge-partial-effect-spaces")`](https://bbuchsbaum.github.io/dkge/articles/dkge-partial-effect-spaces.md)
for the runnable workflow and
[`vignette("dkge-weighting")`](https://bbuchsbaum.github.io/dkge/articles/dkge-weighting.md)
for the distinction among effect, subject, spatial, and transport
weights.

## Which claim does each analysis support?

DKGE results live at several levels. Match the tool and inference unit
to the claim rather than treating “significant DKGE” as a single
outcome.

| Level | Question | Typical tool | What the result does not establish |
|----|----|----|----|
| Component diagnostics | How much fitted variation is retained, and how sensitive is the subspace to refitting choices? | eigenspectrum, identity comparison, [`dkge_plot_subspace_stability()`](https://bbuchsbaum.github.io/dkge/reference/dkge_plot_subspace_stability.md) | inferential reliability or that a particular contrast is non-zero |
| Aggregate component inference | Is a prespecified aggregate component statistic unusual under its resampling null? | [`dkge_aggregate_permute()`](https://bbuchsbaum.github.io/dkge/reference/dkge_aggregate_permute.md), [`dkge_aggregate_bootstrap()`](https://bbuchsbaum.github.io/dkge/reference/dkge_aggregate_bootstrap.md) | validity for the subject-wise q-space fit or an untested contrast |
| Contrast | Does a prespecified effect-space contrast produce a reproducible field? | `dkge_contrast(method = "loso")`, bootstrap or analytic contrast inference | spatial localization without transport and multiplicity control |
| Feature | Where is a component or contrast expressed after alignment? | [`dkge_component_stats()`](https://bbuchsbaum.github.io/dkge/reference/dkge_component_stats.md), [`dkge_transport_contrasts_to_medoid()`](https://bbuchsbaum.github.io/dkge/reference/dkge_transport_contrasts_to_medoid.md) | that the latent direction was selected without bias |
| Between-subject term | Is a subject-level covariate associated with a named multivariate target? | [`dkge_between_rrr()`](https://bbuchsbaum.github.io/dkge/reference/dkge_between_rrr.md), [`dkge_between_permute()`](https://bbuchsbaum.github.io/dkge/reference/dkge_between_permute.md) | a test of the DKGE component itself, or a causal effect without identification assumptions |

For population claims, the subject is the resampling unit. Clusters,
anchors, voxels, and factorial cells within a subject are not
independent subjects. Cross-fitting protects held-out scoring from basis
reuse, but it is not a replacement for uncertainty estimation or
multiple-testing correction.

A prespecified contrast does not need to pass through a
component-significance gate first. Component diagnostics, q-space
contrast inference, transported feature inference, aggregate
decomposition, and between-subject term tests answer different
questions; use the branch that matches the estimand.

Likewise, a contrast can describe a controlled model comparison in
observational data without identifying a causal effect. Causal language
requires design and identification assumptions outside DKGE itself.

## What should you decide before fitting?

Use this sequence:

1.  **Name the estimand.** Decide whether you want subject-level group
    structure from
    [`dkge()`](https://bbuchsbaum.github.io/dkge/reference/dkge.md) /
    [`dkge_fit()`](https://bbuchsbaum.github.io/dkge/reference/dkge_fit.md)
    or an aggregate cell-mean decomposition from
    [`dkge_aggregate_fit()`](https://bbuchsbaum.github.io/dkge/reference/dkge_aggregate_fit.md).
2.  **Audit the input contract.** Verify effect names, observed-cell
    masks, cluster weights, and the subject-level unit of analysis.
3.  **Fit the identity baseline.** Save its eigenspectrum and leading
    saliences.
4.  **Add only justified kernel structure.** List every requested term
    in
    [`design_kernel()`](https://bbuchsbaum.github.io/dkge/reference/design_kernel.md)
    and document its scaling.
5.  **Compare fits in the K metric.** Use
    [`dkge_principal_angles_K()`](https://bbuchsbaum.github.io/dkge/reference/dkge_principal_angles_K.md)
    and compare interpreted saliences, not raw Euclidean cross-products.
6.  **Prespecify the claim level.** Component, contrast, feature, or
    between-subject; then use inference designed for that level.
7.  **Transport only with defensible features.** Record coordinates,
    masses, mapper diagnostics, and reference choice.

## See also

- [`vignette("dkge-workflow")`](https://bbuchsbaum.github.io/dkge/articles/dkge-workflow.md)
  — end-to-end pipeline
- [`vignette("dkge-design-kernels")`](https://bbuchsbaum.github.io/dkge/articles/dkge-design-kernels.md)
  — kernel construction and CV-based selection
- [`vignette("dkge-contrasts-inference")`](https://bbuchsbaum.github.io/dkge/articles/dkge-contrasts-inference.md)
  — LOSO, analytic, and bootstrap inference details
- [`vignette("dkge-components")`](https://bbuchsbaum.github.io/dkge/articles/dkge-components.md)
  — component interpretation and rotation
- [`vignette("dkge-vs-pls")`](https://bbuchsbaum.github.io/dkge/articles/dkge-vs-pls.md)
  — detailed comparison with classical PLSC
- [`vignette("dkge-between-subjects")`](https://bbuchsbaum.github.io/dkge/articles/dkge-between-subjects.md)
  — subject-level RRR and permutation tests
- [`vignette("dkge-partial-effect-spaces")`](https://bbuchsbaum.github.io/dkge/articles/dkge-partial-effect-spaces.md)
  — partial cell coverage, effect-precision weighting, and trialwise
  debiasing
