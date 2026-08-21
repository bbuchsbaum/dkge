# dkge

<!-- badges: start -->
[![R-CMD-check](https://github.com/bbuchsbaum/dkge/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bbuchsbaum/dkge/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/bbuchsbaum/dkge/actions/workflows/pkgdown.yaml/badge.svg)](https://bbuchsbaum.github.io/dkge/)
<!-- badges: end -->

Design-Kernel Group Embedding (DKGE) turns subject-level GLM outputs into a shared, design-aware latent space. It preserves the structure of experimental designs, supports cross-validated contrasts, and provides transport utilities for mapping parcellated fields onto common anchor or voxel representations.

## What it does
- **Design kernels** encode factorial structure, smoothness, and interactions, which control how effects align across subjects.
- **Contrasts and inference** use leave-one-subject-out (LOSO) or K-fold cross-fitting, analytic approximations, and bootstrap utilities for medoid or voxel maps.
- **Transport and rendering** map parcellated fields to a common space with barycentric kNN or C++-accelerated Sinkhorn mappers, anchor graph smoothing, and voxel decoders. The *medoid* is the reference subject whose parcellation the others are mapped onto; it is an index you supply, defaulting to subject 1.
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

**Spatial mapping** — `vignette("dkge-dense-rendering")`, `vignette("dkge-anchors")`, `vignette("dkge-performance")`

**Extras** — `vignette("dkge-plotting")`, `vignette("dkge-cpca")`, `vignette("dkge-vs-pls")`

## Helper constructors

DKGE now provides small helper constructors that validate common orchestration inputs.
They shorten calls to `dkge_pipeline()` and prediction helpers while keeping backward
compatibility with raw lists.

```r
kernel <- diag(nrow(betas[[1]]))
transport <- dkge_transport_spec(
  centroids = centroids,
  sizes = sizes,
  medoid = 2
)
inference <- dkge_inference_spec(B = 1000, tail = "two.sided")
cls_spec <- dkge_classification_spec(targets = ~ condition, method = "lda")

results <- dkge_pipeline(
  betas = betas,
  designs = designs,
  kernel = kernel,
  contrasts = contrasts,
  transport = transport,
  inference = inference,
  classification = cls_spec
)
```

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
