# dkge

Design-Kernel Group Embedding (DKGE) turns subject-level GLM outputs into a shared, design-aware latent space. It preserves the structure of experimental designs, supports cross-validated contrasts, and provides transport utilities for mapping parcellated fields onto common anchor or voxel representations.

## Key capabilities
- **Flexible design kernels** – encode factorial structure, smoothness, and interactions to control how effects align across subjects.
- **Robust contrasts and inference** – LOSO/K-fold cross-fitting, analytic approximations, and bootstrap utilities for medoid or voxel maps.
- **Transport & rendering** – barycentric kNN and C++-accelerated Sinkhorn mappers with warm starts, anchor graph smoothing, and voxel decoders.
- **Component interpretability** – convenience helpers for projecting new data, rotating components, and summarising variance explained.

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
fit <- dkge(betas, designs, K = diag(4), rank = 2)

# project subjects into component space
scores <- dkge_project_btil(fit, fit$Btil)
str(scores, max.level = 1)
```
See the vignettes for full workflows:

- `vignette("dkge-workflow")`
- `vignette("dkge-design-kernels")`
- `vignette("dkge-contrasts-inference")`
- `vignette("dkge-dense-rendering")`
- `vignette("dkge-components")`
- `vignette("dkge-performance")`
- `vignette("dkge-weighting")`

## Documentation & support
Rendered articles and function reference are available at the pkgdown site: <https://bbuchsbaum.github.io/dkge/>. Issues and feature requests are welcome on the [GitHub tracker](https://github.com/bbuchsbaum/dkge/issues).

## Development
- Pull requests are encouraged; please accompany user-facing changes with tests and documentation.
- For large feature work, open an issue to discuss design choices before implementation.

## License
MIT License. See `LICENSE` for details.
