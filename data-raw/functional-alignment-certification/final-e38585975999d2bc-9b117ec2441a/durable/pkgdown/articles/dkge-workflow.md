# DKGE Workflow: From Subject Betas to a Shared Map

Suppose each subject has a GLM beta matrix in a common effect space but
a different spatial parcellation. Your practical goal is to learn a
shared effect-space representation, evaluate a planned contrast on
held-out subjects, and place those subject-specific contrast fields on
one reference parcellation.

This page follows that one path. It uses synthetic data so the code is
self-contained; replace the beta matrices, design matrices, and
centroids with your own. If `q by P_s` and “design kernel” are not yet
familiar, begin with
[`vignette("dkge")`](https://bbuchsbaum.github.io/dkge/articles/dkge.md).

The workflow has two coordinate systems. The fit and contrast are
defined in the shared **effect space**; the final map is defined in a
shared **spatial space**. Keeping those stages separate is the central
organizing idea:

``` text
validate inputs -> optionally regularize spatial fields inside the fit -> fit in effect space -> name a contrast -> cross-fit subject fields -> transport to one map
```

Each section below ends by naming the object you now have and the
operation that becomes valid next.

## What are the inputs?

We use a 2 by 3 factorial design with a planted condition and load
signal. Six subjects have different numbers of clusters, which is the
reason transport is needed later.

``` r
toy <- dkge_sim_toy(
  factors = list(condition = list(L = 2), load = list(L = 3)),
  active_terms = c("condition", "load"),
  S = 6, P = c(16, 18, 15, 17, 16, 19), snr = 4, seed = 2024
)
vapply(toy$B_list, dim, integer(2))
#>      [,1] [,2] [,3] [,4] [,5] [,6]
#> [1,]    5    5    5    5    5    5
#> [2,]   16   18   15   17   16   19
```

The output confirms a common number of effect rows and subject-specific
cluster counts. In real data, beta-row names must match design-column
names.
[`dkge_data()`](https://bbuchsbaum.github.io/dkge/reference/dkge_data.md)
validates that contract and preserves subject IDs.

![Bar chart showing that all six subjects share five effects but have
different numbers of
clusters.](dkge-workflow_files/figure-html/input-shapes-1.png)

The bars differ because each subject has a different spatial partition.
The dashed line marks the five-row effect dimension that is common to
every beta matrix. DKGE can pool information along that common dimension
without pretending that cluster 7 in one subject is cluster 7 in
another.

``` r
bundle <- dkge_data(
  toy$B_list, designs = toy$X_list, subject_ids = toy$subject_ids
)
c(subjects = length(bundle$subject_ids), effects = bundle$q)
#> subjects  effects 
#>        6        5
```

**Next operation:** fit a low-rank group basis. If your subjects lack
effects or have very unequal cell precision, stop here and read
[`vignette("dkge-partial-effect-spaces")`](https://bbuchsbaum.github.io/dkge/articles/dkge-partial-effect-spaces.md)
first.

If unsmoothed beta columns need a spatial prior, define it before
fitting with
[`dkge_spatial_regularizer()`](https://bbuchsbaum.github.io/dkge/reference/dkge_spatial_regularizer.md)
and pass it through `spatial =`. Although DKGE’s eigensolve is in effect
space, the regularized beta fields change the pooled effect-space moment
and can therefore change the learned basis. This is not the same as
smoothing a rendered map after fitting; see
[`vignette("dkge-spatial-regularization")`](https://bbuchsbaum.github.io/dkge/articles/dkge-spatial-regularization.md)
for graph construction, held-out selection of `lambda`, and inference
boundaries.

## What is the baseline fit?

Start with an identity kernel and no spatial regularizer. It retains the
GLM/design scaling but adds neither kernel-imposed similarity between
effects nor graph-imposed similarity between spatial units.

``` r
K_identity <- diag(bundle$q)
fit_identity <- dkge(bundle, K = K_identity, rank = 2, w_method = "none")
dkge_plot_scree(fit_identity)
```

![Scree plot for the two-component identity-kernel baseline
fit.](dkge-workflow_files/figure-html/fit-identity-1.png)

This fit returns a two-column group basis `fit_identity$U`. The scree
plot describes variation represented by those components; it does not
establish their reliability or scientific importance.

**Output:** an unconstrained effect-space baseline.

**Next operation:** fit the scientifically motivated kernel and compare
it with this baseline.

## What changes when the design kernel is used?

The simulator supplies the structured kernel used to plant the signal. A
real analysis would construct this object explicitly with
[`design_kernel()`](https://bbuchsbaum.github.io/dkge/reference/design_kernel.md)
and prespecify the included terms. Both fits disable subject block
weighting so the comparison isolates the kernel change.

``` r
fit_structured <- dkge(bundle, K = toy$K, rank = 2, w_method = "none")

comparison <- data.frame(
  component = 1:2,
  identity = dkge_variance_explained(fit_identity)$prop_var,
  structured = dkge_variance_explained(fit_structured)$prop_var
)
round(comparison, 3)
#>   component identity structured
#> 1         1    0.599       0.58
#> 2         2    0.401       0.42
```

![Grouped bar chart comparing variance proportions for identity and
structured design-kernel
fits.](dkge-workflow_files/figure-html/compare-kernel-variance-1.png)

Compare the retained subspaces as well as the scree values. Because the
bases use a kernel metric, a raw
[`crossprod()`](https://rdrr.io/pkg/Matrix/man/matmult-methods.html) is
not a valid rotation comparison.

``` r
angles <- dkge_principal_angles_K(
  fit_identity$U, fit_structured$U, fit_structured$K
)
round(angles * 180 / pi, 1)
#> [1] 0 0
```

Here component 1 accounts for 59.9% of the retained variation under the
identity kernel and 58% under the structured kernel. Both principal
angles round to 0 and 0 degrees: the kernel changes how variation is
allocated within this seeded two-dimensional subspace, but not the
subspace itself at displayed precision. This is a sensitivity
comparison, not a test of which kernel is true.

**Output:** a structured fit plus an explicit record of how it differs
from the identity baseline.

**Next operation:** inspect effect saliences, then specify the contrast
you actually want to estimate.

## Which effect directions define the components?

``` r
dkge_plot_effect_loadings(fit_structured, comps = 1:2)
```

![Heatmap of effect-space saliences for the first two components under
the structured design
kernel.](dkge-workflow_files/figure-html/show-effect-loadings-1.png)

The heatmap displays $`K U`$: effect-space saliences, not voxel
coefficients. Look for coherent relative patterns within a component. Do
not compare a component’s arbitrary sign across independent fits without
alignment.

## How do you obtain held-out contrast fields?

Build the contrast in fitted effect order. Here the simulator identifies
the coordinate for the planted condition term.

``` r
condition_contrast <- numeric(bundle$q)
condition_contrast[toy$active_cols$condition] <- 1

condition_loso <- dkge_contrast(
  fit_structured,
  list(condition = condition_contrast),
  method = "loso"
)
condition_loso
#> DKGE Contrasts
#> --------------
#> Method: loso
#> Contrasts: 1
#> Subjects: 6
#> Contrast names: condition
```

The output is a named list of six subject-specific cluster vectors. LOSO
cross-fitting keeps each subject out of the basis used to score that
subject. It does not make a randomized causal contrast out of an
observational design, nor does it supply group uncertainty by itself.

**Output:** one cross-fitted contrast field per subject, still indexed
by that subject’s clusters.

**Next operation:** transport those fields before stacking or performing
feature-level group inference.

## How do subject fields reach one parcellation?

Transport needs real spatial features, normally cluster centroids in a
common coordinate system. The synthetic coordinates below exist only to
demonstrate the interface; their resulting map has no anatomical
interpretation.

``` r
centroids <- lapply(seq_along(toy$B_list), function(s) {
  p <- ncol(toy$B_list[[s]])
  cbind(x = seq_len(p), y = rep(s / 100, p), z = rep(0, p))
})
```

These coordinates are enough for geometry-only resampling, but this
synthetic example contains no independent functional acquisition from
which to identify parcel correspondence. It therefore stops before
feature-level group alignment rather than silently treating subject 1 or
a bare grid as functional truth.

For the complete path—independent kernel-image features, held-out
reference selection, a group functional template, aligned subject rows,
inference, and rendering—continue with
[`vignette("dkge-functional-alignment")`](https://bbuchsbaum.github.io/dkge/articles/dkge-functional-alignment.md).
The key output there is an `S` by `Q` matrix for each contrast: subjects
are the inference units, and columns are locations on one identified
support.

## What should you save and report?

For a reproducible analysis, retain:

- effect names and the design matrices used to produce each beta matrix;
- the identity and structured kernels, their construction parameters,
  and the chosen rank;
- the identity-versus-structured comparison;
- cross-fitting folds or LOSO specification;
- cluster coordinates, masses, typed feature provenance, mapper
  settings, reference-selection criterion, template trajectory, and
  transport diagnostics; and
- the inference unit and multiplicity procedure for every reported
  claim.

Continue with
[`vignette("dkge-concepts")`](https://bbuchsbaum.github.io/dkge/articles/dkge-concepts.md)
for the fitted moment and inference levels,
[`vignette("dkge-contrasts-inference")`](https://bbuchsbaum.github.io/dkge/articles/dkge-contrasts-inference.md)
for uncertainty, and
[`vignette("dkge-functional-alignment")`](https://bbuchsbaum.github.io/dkge/articles/dkge-functional-alignment.md)
for cross-subject correspondence.
[`vignette("dkge-weighting")`](https://bbuchsbaum.github.io/dkge/articles/dkge-weighting.md)
distinguishes spatial, effect, subject, and transport weights.
