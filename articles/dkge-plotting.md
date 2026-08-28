# Plotting DKGE Fits

You have a fitted DKGE model and want to see what it found. Four panels
answer that, each a descriptive question about the fit: how much
variation was retained, which effects define a component, whether one
subject dominates, and how stable the subspace is under refitting. None
supplies a p-value or validates the scientific model on its own, and a
polished figure is easy to mistake for evidence that it does.

You need a fitted `dkge` object. This page builds one from a small
simulation so that every panel is reproducible.

## Prerequisites

``` r

library(dkge)
library(ggplot2)
library(patchwork)
```

`ggplot2` and `patchwork` are required. `ggrepel` is optional: the
information-map functions label top anchors with it when it is
installed, and with base text otherwise.

## Simulate a toy dataset

Each subject gets a diagonal design matrix, one column per effect, plus
modest signal and noise so the latent components are recoverable without
being trivial.
[`dkge_sim_toy()`](https://bbuchsbaum.github.io/dkge/reference/dkge_sim_toy.md)
is not used here because the plots need deterministic ground truth for a
fixed four-effect design.

``` r

q <- 4L   # number of effects
P <- 30L  # clusters (voxels) per subject
S <- 8L   # subjects

make_subject <- function(id) {
  design <- diag(q)
  colnames(design) <- paste0('eff', seq_len(q))
  signal <- matrix(rnorm(q * P, sd = 0.4), nrow = q)
  noise <- matrix(rnorm(q * P, sd = 0.2), nrow = q)
  dkge_subject(signal + noise, design = design, id = paste0('sub', id))
}

subjects <- lapply(seq_len(S), make_subject)
K <- diag(q)
fit <- dkge(subjects, K = K, rank = 3, w_method = 'mfa_sigma1')
```

## Auxiliary objects

To explain the stability plot’s input contract, we create perturbed
versions of the fitted rank-three basis and re-orthonormalize them in
the $`K`$-metric. These are controlled sensitivity perturbations, not
cross-validation folds. In an analysis, supply bases actually refitted
inside the relevant folds.

``` r

fit_U <- fit[['U']]
bases <- replicate(4, {
  perturb <- matrix(rnorm(length(fit_U), sd = 0.02), nrow = nrow(fit_U))
  dkge_k_orthonormalize(fit_U + perturb, fit[['K']])
}, simplify = FALSE)
base_labels <- paste0('base', seq_along(bases))
```

The optional information-map panels require real outputs from
[`dkge_info_map_haufe()`](https://bbuchsbaum.github.io/dkge/reference/dkge_info_map_haufe.md)
or
[`dkge_info_map_loco()`](https://bbuchsbaum.github.io/dkge/reference/dkge_info_map_loco.md).
We omit them here because random vectors would demonstrate only the
plotting signature while looking like scientific evidence. See
[`vignette("dkge-classification")`](https://bbuchsbaum.github.io/dkge/articles/dkge-classification.md)
for the classifier workflow that precedes those attribution methods.

## Individual plots

### Scree

The scree plot shows each retained component’s share of the fitted
variation. The marked rank below is chosen for illustration; in an
analysis you would obtain it from
[`dkge_cv_rank_loso()`](https://bbuchsbaum.github.io/dkge/reference/dkge_cv_rank_loso.md)
or
[`dkge_cv_kernel_rank()`](https://bbuchsbaum.github.io/dkge/reference/dkge_cv_kernel_rank.md)
rather than by eye.

``` r

one_se_pick <- 3
dkge_plot_scree(fit, one_se_pick = one_se_pick)
```

![Scree plot of fitted variation with the illustrative rank-three
selection marked.](dkge-plotting_files/figure-html/scree-1.png)

### Effect-space loadings

``` r

dkge_plot_effect_loadings(fit, comps = 1:3, zscore = TRUE)
```

![Heatmap of standardized effect-space saliences for three DKGE
components.](dkge-plotting_files/figure-html/loadings-1.png)

Rows are effects, columns are components. Read down a column: effects
with similar color and magnitude move together in that component, and
effects with opposite colors define a contrast-like pattern. Because
`zscore = TRUE` standardizes loadings within each effect, compare cells
along a row rather than across rows; the display shows where each effect
is expressed, not how large its raw loading is.

### Subject contributions

[`dkge_plot_subject_contrib()`](https://bbuchsbaum.github.io/dkge/reference/dkge_plot_subject_contrib.md)
returns two linked panels. The left panel shows the subject-level
weights used while fitting the model (`fit$weights`). Those weights
depend on the `w_method` argument passed to
[`dkge()`](https://bbuchsbaum.github.io/dkge/reference/dkge.md). This
fit uses `w_method = "mfa_sigma1"`, the MFA-style default, so the bars
vary around one: a subject whose block carries a large leading singular
value is downweighted. `w_method = "none"` would give every subject unit
weight and a flat row of bars. The heatmap on the right shows how much
norm (“energy”) each subject contributes to the selected components
after those weights have been applied.

``` r

contrib <- dkge_plot_subject_contrib(fit, comps = 1:3)
contrib$weights + contrib$energy + patchwork::plot_layout(widths = c(1, 2))
```

![Linked bar and heatmap panels showing subject fitting weights and
component contribution
energy.](dkge-plotting_files/figure-html/contrib-1.png)

### Subspace stability

This diagnostic compares each supplied basis (controlled perturbations
here) against a consensus basis by computing principal angles in the
$`K`$ metric. Smaller angles mean that a base reproduces the consensus
component more closely, so parallel lines near zero indicate a stable
subspace, whereas large excursions highlight folds or components that
deviate materially.

``` r

dkge_plot_subspace_stability(bases, K = fit[['K']], labels = base_labels)
```

![Principal-angle stability curves comparing four perturbed bases with
the consensus basis.](dkge-plotting_files/figure-html/stability-1.png)

## Combine the fit diagnostics

``` r

dkge_plot_suite(fit,
                one_se_pick = one_se_pick,
                comps = 1:3,
                bases = bases,
                consensus = fit[['U']],
                base_labels = base_labels,
                top = 5)
```

![Dashboard combining scree, effect saliences, subject contributions,
and subspace stability
diagnostics.](dkge-plotting_files/figure-html/suite-1.png)

[`dkge_plot_suite()`](https://bbuchsbaum.github.io/dkge/reference/dkge_plot_suite.md)
arranges the available panels and leaves the attribution row empty when
Haufe/LOCO results are absent. Before reporting the dashboard, replace
the sensitivity perturbations with actual refitted bases and interpret
each panel next to the numerical diagnostic that produced it.

## Saving the dashboard

``` r

dkge_plot_suite(fit,
                bases = bases,
                consensus = fit[['U']],
                base_labels = base_labels,
                save_path = "dkge_dashboard.png",
                width = 10,
                height = 10)
```

## What these panels do and do not establish

Use the scree plot to describe retained variation, the salience heatmap
to name effect directions, the contribution panels to detect subject
dominance, and principal angles to summarize refit sensitivity. Save the
figure only after the objects behind those panels come from the analysis
being reported. Visual coherence improves communication; it does not
upgrade descriptive diagnostics into inferential evidence.

## Where to go next

- [`vignette("dkge-components")`](https://bbuchsbaum.github.io/dkge/articles/dkge-components.md)
  reads the same saliences as scientific claims rather than as panels.
- [`vignette("dkge-contrasts-inference")`](https://bbuchsbaum.github.io/dkge/articles/dkge-contrasts-inference.md)
  is the page that supplies the evidence these panels deliberately do
  not.
- [`vignette("dkge-dense-rendering")`](https://bbuchsbaum.github.io/dkge/articles/dkge-dense-rendering.md)
  maps a component into a shared brain space once you have one worth
  mapping.
