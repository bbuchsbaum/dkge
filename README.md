# dkge

<!-- badges: start -->
[![R-CMD-check](https://github.com/bbuchsbaum/dkge/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bbuchsbaum/dkge/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/bbuchsbaum/dkge/actions/workflows/pkgdown.yaml/badge.svg)](https://bbuchsbaum.github.io/dkge/)
<!-- badges: end -->

Design-Kernel Group Embedding (DKGE) turns subject-level GLM outputs into a shared, design-aware latent space. It preserves the structure of experimental designs, supports cross-validated contrasts, and provides transport utilities for mapping parcellated fields onto common anchor or voxel representations.

## What it does
- **Design kernels** encode factorial structure, effect-space smoothness, and interactions, which control how effects align across subjects.
- **Model-level spatial regularization** uses sparse graph-Laplacian solves to smooth subject fields inside the pooled moment, so the spatial prior can change the learned basis as well as its reconstructed maps.
- **Contrasts and inference** use leave-one-subject-out (LOSO) or K-fold cross-fitting. Rank-truncated cohort-trained inference is labelled approximate and requires explicit opt-in; aligned-map bootstraps require typed correspondence provenance.
- **Functional alignment and rendering** use typed independent response signatures, auditable reference selection, or an iterative group template to map subject fields onto one identified support. A subject is a *medoid* only when selected by a stated criterion; a bare MNI grid supplies coordinates, not functional correspondence. Rendering is downstream of alignment.
- **Classifier localization** cross-fits latent classifiers and returns decoder, Haufe, and LOCO maps.
- **Component interpretation** projects new data, rotates components, and summarizes variance explained.

## Installation
```r
# install.packages("remotes")
remotes::install_github("bbuchsbaum/dkge")
```
The package depends on `RcppArmadillo`, `future`, `multivarious`, and other CRAN libraries; these install automatically.

## Getting started
```r
library(dkge)

# simulate three subjects with four effects and five clusters
set.seed(1)
betas <- replicate(3, matrix(rnorm(4 * 5), 4, 5), simplify = FALSE)
designs <- replicate(3, qr.Q(qr(matrix(rnorm(60 * 4), 60, 4))), simplify = FALSE)

# fit DKGE with an identity kernel and rank 2
fit <- dkge(betas, designs, kernel = diag(4), rank = 2)

# project subjects into component space
scores <- dkge_project_btil(fit, fit$Btil)
str(scores, max.level = 1)
```
Start with `vignette("dkge")`, then `vignette("dkge-workflow")`. The full set:

**Start here** — `vignette("dkge")`, `vignette("dkge-workflow")`, `vignette("dkge-concepts")`

**Core analysis** — `vignette("dkge-design-kernels")`, `vignette("dkge-contrasts-inference")`, `vignette("dkge-components")`, `vignette("dkge-classification")`

**Study designs** — `vignette("dkge-partial-effect-spaces")`, `vignette("dkge-unbalanced-trialwise")`, `vignette("dkge-between-subjects")`

**Weighting** — `vignette("dkge-weighting")`, `vignette("dkge-adaptive-weighting")`

**Spatial mapping** — `vignette("dkge-functional-alignment")`, `vignette("dkge-spatial-regularization")`, `vignette("dkge-dense-rendering")`, `vignette("dkge-anchors")`, `vignette("dkge-performance")`

**Extras** — `vignette("dkge-plotting")`, `vignette("dkge-cpca")`, `vignette("dkge-vs-pls")`

## Helper constructors

DKGE now provides small helper constructors that validate common orchestration inputs.
They shorten calls to `dkge_pipeline()` and prediction helpers while keeping backward
compatibility with raw lists.

```r
kernel <- diag(nrow(betas[[1]]))
contrasts <- c(1, -1, 0, 0)
inference <- dkge_inference_spec(
  B = 1000,
  tail = "two.sided",
  allow_approximate_alignment = TRUE
)
cls_spec <- dkge_classification_spec(targets = ~ condition, method = "lda")

results <- dkge_pipeline(
  betas = betas,
  designs = designs,
  kernel = kernel,
  contrasts = contrasts,
  inference = inference,
  classification = cls_spec
)
```

The opt-in is explicit because the rank-truncated latent span is estimated from
the same cohort. The returned inference object labels that estimator
`"approximate"`; omit the opt-in to fail closed. All requested contrasts and
support locations share one max-T family.

For cross-subject functional correspondence, use the typed
`dkge_prepare_alignment()` / `dkge_transport_contrasts_to_reference()` workflow
in `vignette("dkge-functional-alignment")`; the pipeline does not infer a
functional mapping from coordinates alone.

To score new subjects without manually assembling `B_list`, use
`dkge_predict_subjects()`:

```r
pred <- dkge_predict_subjects(fit, betas = new_subjects, contrasts = my_contrasts)
```


## Documentation & support
Project home and documentation: <https://github.com/bbuchsbaum/dkge>. Issues and feature requests are welcome on the [GitHub tracker](https://github.com/bbuchsbaum/dkge/issues).

## Development

Pull requests are encouraged. See [CONTRIBUTING.md](CONTRIBUTING.md) for
packaging conventions, the architecture map, and how to propose changes to it.

## License
MIT License. See `LICENSE` for details.

<!-- albersdown:theme-note:start -->
## Albers theme
This package uses the albersdown theme. Existing vignette theme hooks are replaced so `albers.css` and local `albers.js` render consistently on CRAN and GitHub Pages. The defaults are configured via `params$family` and `params$preset` (family = 'red', preset = 'interaction'). The pkgdown site uses `template: { package: albersdown }` together with generated `pkgdown/extra.css` and `pkgdown/extra.js` so the theme is linked and activated on site pages.
<!-- albersdown:theme-note:end -->
