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

The decomposition is easiest to understand as a sequence:

``` text
subject beta blocks -> small effect-by-effect moments -> pooled raw moment -> design/kernel geometry -> components
```

The first arrow summarizes spatial covariation within each subject. The
middle step combines subjects. Only then does the design kernel say
which directions in effect space should count as nearby or important.
This ordering is why missingness and precision must be handled before
the kernel transform.

## What enters the decomposition?

Terms used throughout this page and the rest of the suite are defined in
`?dkge-glossary`: effect space, design kernel, salience, cross-fitting,
reference support, functional correspondence, template, aligned subject
rows, inference, rendering, medoid, transport, and the estimand each
analysis targets.

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

In plain language, (M) records which effects tend to have large spatial
patterns together. (R) puts subjects’ design estimates on a common
ruler. (K) defines the scientific geometry among effect directions. The
eigensolve then finds a small set of directions that summarize the
transformed moment.

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

The three matrices below make the transformation concrete. Each panel
has its own color scale because the scientific question is the pattern
within a matrix, not equality of raw numerical ranges across stages.

![Three heatmaps showing the pooled raw effect moment, the design
kernel, and the transformed moment used for the
eigendecomposition.](dkge-concepts_files/figure-html/moment-panels-1.png)

The left panel comes from the data; the middle panel encodes the chosen
design geometry; the right panel is the matrix actually decomposed. A
strong pattern in the right panel is therefore a joint consequence of
empirical covariance and modeling choices, not a feature discovered
independently of the kernel.

**Input:** subject effect moments.

**Output:** a low-rank basis in the kernel metric.

**Next operation:** inspect $`K U`$ with
[`dkge_component_saliences()`](https://bbuchsbaum.github.io/dkge/reference/dkge_component_saliences.md),
or project a prespecified contrast with
[`dkge_contrast()`](https://bbuchsbaum.github.io/dkge/reference/dkge_contrast.md).

## What does the kernel change?

The kernel sets the metric. Component directions are estimated under it:
the kernel changes which effect-space directions count as large or
smooth during the eigensolve; the empirical moment still determines the
fitted components.

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
for the coverage contract,
[`vignette("dkge-unbalanced-trialwise")`](https://bbuchsbaum.github.io/dkge/articles/dkge-unbalanced-trialwise.md)
for the runnable workflow, and
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
| Feature | Where is a prespecified contrast expressed after identified alignment? | [`dkge_transport_contrasts_to_reference()`](https://bbuchsbaum.github.io/dkge/reference/dkge_transport_contrasts_to_reference.md), [`dkge_infer_aligned()`](https://bbuchsbaum.github.io/dkge/reference/dkge_infer_aligned.md) | that the latent direction was selected without bias or that an approximate alignment is exact |
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

## How the pieces fit together

A DKGE analysis moves through four stages, and each one has a public
entry point:

1.  **Harmonize** subject data with
    [`dkge_subject()`](https://bbuchsbaum.github.io/dkge/reference/dkge_subject.md)
    and
    [`dkge_data()`](https://bbuchsbaum.github.io/dkge/reference/dkge_data.md),
    so every subject’s beta rows carry the same effect labels.
2.  **Fit** a group embedding with
    [`dkge()`](https://bbuchsbaum.github.io/dkge/reference/dkge.md),
    which compresses to effect space and solves for the shared basis
    under the design kernel.
3.  **Cross-fit** a prespecified contrast with
    [`dkge_contrast()`](https://bbuchsbaum.github.io/dkge/reference/dkge_contrast.md),
    so no subject’s own data shapes the basis used to score it.
4.  **Align** the resulting fields onto an identified support with typed
    functional features,
    [`dkge_prepare_alignment()`](https://bbuchsbaum.github.io/dkge/reference/dkge_prepare_alignment.md),
    and
    [`dkge_transport_contrasts_to_reference()`](https://bbuchsbaum.github.io/dkge/reference/dkge_transport_contrasts_to_reference.md).

The effect-space basis, reference support, fitted correspondence,
functional template, aligned subject rows, group inference, and renderer
are distinct objects. In particular, a coordinate grid has no functional
mapping until features have been learned on it. See
[`vignette("dkge-functional-alignment")`](https://bbuchsbaum.github.io/dkge/articles/dkge-functional-alignment.md).

[`dkge_pipeline()`](https://bbuchsbaum.github.io/dkge/reference/dkge_pipeline.md)
runs all four in one call. `CONTRIBUTING.md` in the package repository
carries the implementation map for contributors: which source file owns
each stage, and the rules for changing them.

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

The next practical pages are
[`vignette("dkge-design-kernels")`](https://bbuchsbaum.github.io/dkge/articles/dkge-design-kernels.md)
for kernel construction,
[`vignette("dkge-contrasts-inference")`](https://bbuchsbaum.github.io/dkge/articles/dkge-contrasts-inference.md)
for contrast uncertainty, and
[`vignette("dkge-between-subjects")`](https://bbuchsbaum.github.io/dkge/articles/dkge-between-subjects.md)
for subject-level associations.
