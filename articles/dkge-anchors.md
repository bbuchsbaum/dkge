# Feature-Anchored DKGE

## Overview

This vignette describes how to run DKGE when each subject is observed on
a different set of stimuli but every item carries a feature vector in a
shared space (e.g., a 100-dimensional embedding). The feature-anchored
workflow replaces discrete-cell completion with a common anchor basis
derived directly from the feature space.

We will:

1.  build an anchor descriptor from subject-specific features and item
    kernels,
2.  fit DKGE through the shared `dkge_input` interface, and
3.  evaluate cross-fitted contrasts that operate in the anchor basis.

``` r
library(dkge)
library(Matrix)
set.seed(1)
```

## Simulated feature-aligned data

In this example we simulate three subjects. Each subject has their own
set of item features sampled around four latent prototypes. We generate
subject-specific beta maps by projecting the item responses through SVD
loadings and adding gaussian noise.

``` r
# Number of latent anchors and feature dimension
d <- 20L
anchors_true <- matrix(rnorm(4 * d), 4, d)

make_subject <- function(n_items, n_vox, seed) {
  set.seed(seed)
  # Each subject observes item features around the latent anchors
  centers <- anchors_true[sample.int(nrow(anchors_true), n_items, replace = TRUE), , drop = FALSE]
  features <- centers + matrix(rnorm(n_items * d, sd = 0.4), n_items, d)

  # Build an item similarity kernel (e.g., RSA over betas)
  latent <- matrix(rnorm(n_items * 3), n_items, 3)
  beta_loadings <- matrix(rnorm(3 * n_vox), 3, n_vox)
  betas <- latent %*% beta_loadings + matrix(rnorm(n_items * n_vox, sd = 0.2), n_items, n_vox)
  item_kernel <- betas %*% t(betas)

  list(features = features,
       item_kernel = item_kernel)
}

subjects <- list(
  s1 = make_subject(30, 120, seed = 11),
  s2 = make_subject(45, 120, seed = 12),
  s3 = make_subject(35, 120, seed = 13)
)

features_list <- lapply(subjects, `[[`, "features")
K_item_list  <- lapply(subjects, `[[`, "item_kernel")
```

## Build an anchor descriptor

We choose 16 anchors via the default d-kpp selector. The descriptor
records both the anchor configuration and the DKGE options to be used
after congruence.

``` r
anchor_input <- dkge_input_anchor(
  features_list = features_list,
  K_item_list = K_item_list,
  anchors = list(L = 16, method = "dkpp", seed = 99L),
  dkge_args = list(w_method = "none")
)
```

## Fit DKGE through the shared interface

[`dkge_fit_from_input()`](https://bbuchsbaum.github.io/dkge/reference/dkge_fit_from_input.md)
converts the descriptor into aligned anchor kernels and calls the
standard fitter. The resulting object is a regular `dkge` fit containing
the anchor provenance.

``` r
fit_anchor <- dkge_fit_from_input(anchor_input)
fit_anchor
#> Multiblock Bi-Projector object:
#>   Projection matrix dimensions:  48 x 16 
#>   Block indices:
#>     Block 1: 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16
#>     Block 2: 17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32
#>     Block 3: 33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48
fit_anchor$provenance$anchors$coverage
#>   subject      p50      p90      p95
#> 1      s1 2.228517 2.568386 2.630797
#> 2      s2 2.262534 2.635445 2.727658
#> 3      s3 2.382655 2.674812 2.714045
```

## Cross-fitted contrasts in the anchor basis

Contrasts are vectors over the anchor index. For illustration we
consider the first latent axis and compute LOSO contrasts.

``` r
contrast_vec <- rep(0, fit_anchor$provenance$anchors$L)
contrast_vec[1] <- 1

res_contrast <- dkge_contrast(fit_anchor,
                               contrasts = list(anchor1 = contrast_vec),
                               method = "loso")
res_contrast$values$anchor1
#> $s1
#>  [1]  0.8578792172  0.1872228773 -0.0575353725  0.0210228694  0.0040815713
#>  [6]  0.0026493434  0.0084611867  0.0060824492 -0.0036419263 -0.0071511478
#> [11]  0.0020458780  0.0051067140  0.0017552751  0.0019931698 -0.0007561988
#> [16]  0.0008524834
#> 
#> $s2
#>  [1]  0.6501518261  0.0408124609  0.0263281662 -0.0020778356  0.0214882891
#>  [6] -0.0063855690  0.0001019468  0.0048395348 -0.0024434787 -0.0009357923
#> [11] -0.0032645685  0.0008471532  0.0003527032 -0.0011270490  0.0016957389
#> [16]  0.0031749291
#> 
#> $s3
#>  [1]  0.7341606369  0.0011397863  0.1100017265 -0.0159284837  0.0019253123
#>  [6]  0.0105681143 -0.0110783417  0.0204310852 -0.0192341142 -0.0088613586
#> [11]  0.0107199583  0.0040217500 -0.0039543365  0.0025623464  0.0002294024
#> [16] -0.0010917674
```

## Using the pipeline helper

The same workflow integrates with
[`dkge_pipeline()`](https://bbuchsbaum.github.io/dkge/reference/dkge_pipeline.md)
by supplying the descriptor via the new `input` argument. All downstream
services (contrasts, classification, inference, transport) operate
exactly as with design-level inputs.

``` r
pipeline_res <- dkge_pipeline(input = anchor_input,
                               contrasts = list(anchor1 = contrast_vec),
                               method = "analytic",
                               inference = NULL)
summary(pipeline_res$contrasts)
#>           Length Class  Mode     
#> values     1     -none- list     
#> method     1     -none- character
#> contrasts  1     -none- list     
#> metadata  13     -none- list
```

## Classification targets

Anchor-based fits do not store the design-factor mapping that
[`dkge_targets()`](https://bbuchsbaum.github.io/dkge/reference/dkge_targets.md)
relies on. When you want to classify anchor effects you must provide
explicit weight matrices (rows = classes, columns = anchors) or
pre-built `dkge_target` objects. The helpers
[`dkge_anchor_targets_from_prototypes()`](https://bbuchsbaum.github.io/dkge/reference/dkge_anchor_targets_from_prototypes.md)
and
[`dkge_anchor_targets_from_directions()`](https://bbuchsbaum.github.io/dkge/reference/dkge_anchor_targets_from_directions.md)
turn feature-space prototypes or directions into the required matrices.

``` r
anchors_mat <- fit_anchor$provenance$anchors$anchors
proto_list <- list(
  classA = anchors_mat[c(1, 2), , drop = FALSE],
  classB = anchors_mat[c(3, 4), , drop = FALSE]
)
target_matrix <- dkge_anchor_targets_from_prototypes(anchors_mat, proto_list)
target_matrix
#>                [,1]         [,2]         [,3]         [,4]         [,5]
#> classA 7.069816e-01 7.069816e-01 3.284593e-12 5.711037e-08 7.286561e-04
#> classB 3.885725e-10 5.673129e-08 7.070586e-01 7.070586e-01 7.422089e-09
#>                [,6]         [,7]         [,8]         [,9]        [,10]
#> classA 5.171745e-03 6.703027e-11 5.793457e-03 6.889859e-08 3.608127e-03
#> classB 3.606838e-11 3.967735e-03 1.593540e-09 1.825886e-03 5.129250e-07
#>               [,11]        [,12]        [,13]        [,14]        [,15]
#> classA 1.798528e-09 1.091749e-08 4.974672e-03 1.048068e-02 1.206072e-02
#> classB 2.632166e-03 9.580194e-03 7.723220e-07 1.189718e-10 1.058633e-11
#>               [,16]
#> classA 2.170235e-07
#> classB 4.312809e-03

# Ready for classification
cls <- dkge_classify(fit_anchor,
                     targets = target_matrix,
                     method = "lda",
                     folds = 2)
cls$summary
#> NULL
```

## Diagnostics and provenance

Anchor coverage, leverage, and bandwidth settings are stored under
`fit$provenance$anchors`. These diagnostics are useful for checking
whether the median heuristic and chosen number of anchors provide
adequate coverage across subjects.

``` r
dkge_anchor_diagnostics(fit_anchor)
#> $summary
#> $summary$method
#> [1] "dkpp"
#> 
#> $summary$sigma
#> [1] 6.084196
#> 
#> $summary$L
#> [1] 16
#> 
#> $summary$mean_item_count
#> [1] 36.66667
#> 
#> 
#> $coverage
#>   subject      p50      p90      p95
#> 1      s1 2.228517 2.568386 2.630797
#> 2      s2 2.262534 2.635445 2.727658
#> 3      s3 2.382655 2.674812 2.714045
#> 
#> $leverage
#>       anchor  leverage
#> 1   anchor_1 1.7532216
#> 2   anchor_2 0.8433377
#> 3   anchor_3 0.9832618
#> 4   anchor_4 0.8512535
#> 5   anchor_5 1.4294242
#> 6   anchor_6 1.2357775
#> 7   anchor_7 0.6482376
#> 8   anchor_8 1.4155902
#> 9   anchor_9 0.5546306
#> 10 anchor_10 0.7626024
#> 11 anchor_11 0.5788092
#> 12 anchor_12 0.3621954
#> 13 anchor_13 0.6351782
#> 14 anchor_14 1.5522472
#> 15 anchor_15 1.2370644
#> 16 anchor_16 1.1571685
```

## Summary

The feature-anchored path allows DKGE to align subjects with disjoint
item sets without imputing missing cells. By sampling a fold-safe set of
anchors, whitening the anchor basis, and delegating back to the core
DKGE routines, the approach keeps the computational core untouched while
offering a modular front-end suitable for modern embedding-based
experiments.

### Special cases: shared and mixed item sets

- **All subjects share the same items.** The anchor pipeline reduces to
  a re-basing of the common item kernel because every subject projects
  onto identical feature rows. You can keep the anchor path for
  consistency (the whitening step simply orthonormalises the shared
  basis) or, if preferred, fall back to
  [`dkge_fit_from_kernels()`](https://bbuchsbaum.github.io/dkge/reference/dkge_fit_from_kernels.md)
  with the shared item kernel—both produce comparable PSD inputs for the
  core fitter.
- **Subgroups with identical items.** When subsets of participants see
  the same stimulus sequence, they automatically receive identical
  anchor projections because
  [`dkge_build_anchor_kernels()`](https://bbuchsbaum.github.io/dkge/reference/dkge_build_anchor_kernels.md)
  selects anchors once per fold from the pooled training subjects.
  Coverage and leverage diagnostics in `fit$provenance$anchors` reveal
  these overlaps; large leverage spikes indicate anchors dominated by a
  subgroup and may motivate a smaller `L` or tighter bandwidth.
