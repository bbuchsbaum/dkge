# Feature-Anchored DKGE

The ordinary DKGE workflow assumes that subjects’ effect rows can be
named on one discrete grid. That assumption fails when each subject saw
a different set of stimuli. If every item instead carries a comparable
feature vector—such as an embedding—you can align subjects through
representative locations in that feature space.

## Why use feature anchors?

Feature anchors solve the *item alignment* problem. They do not solve
spatial alignment of voxels or prove that the feature representation is
scientifically adequate.

``` text
different item sets + shared item features -> pooled feature anchors -> aligned item kernels -> DKGE basis
```

This page will:

1.  build an anchor descriptor from subject-specific features and item
    kernels,
2.  fit DKGE through the shared `dkge_input` interface, and
3.  evaluate LOSO contrasts conditional on that fixed anchor
    representation.

``` r
library(dkge)
#> Registered S3 method overwritten by 'genpca':
#>   method                   from        
#>   transfer.cross_projector multivarious
library(Matrix)
set.seed(1)
```

## Simulated feature-aligned data

The example simulates three subjects. Each has its own item features,
sampled around four latent prototypes, and its own clusters, so no two
subjects share a parcellation.

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

The default d-kpp selector picks 16 anchors. The descriptor records the
anchor configuration and the per-subject kernels built against it.

``` r
anchor_input <- dkge_input_anchor(
  features_list = features_list,
  K_item_list = K_item_list,
  anchors = list(L = 16, method = "dkpp", seed = 99L),
  dkge_args = list(w_method = "none")
)
```

This call chooses anchors from the pooled subjects; no fold-specific
anchor construction is requested. The anchor representation is therefore
fixed before the LOSO analysis below. That is a useful conditional
analysis, but it is not end-to-end cross-validation of feature
preprocessing.

## Fit DKGE through the shared interface

[`dkge_fit_from_input()`](https://bbuchsbaum.github.io/dkge/reference/dkge_fit_from_input.md)
converts the descriptor into aligned anchor kernels and calls the
standard fitter, so the anchor path reuses the ordinary fit rather than
forking it.

``` r
fit_anchor <- dkge_fit_from_input(anchor_input)
fit_anchor
#> <dkge>
#>   Subjects: 3 
#>   Effects: 16 
#>   Rank: 16 
#>   Subject weighting: none (tau = 0.3) 
#>   Weight range: 1 to 1 (median = 1, CV = 0) 
#>   Effective subject mass: 3 of 3 usable
fit_anchor$provenance$anchors$coverage
#>   subject      p50      p90      p95
#> 1      s1 2.343931 2.524823 2.591260
#> 2      s2 2.292074 2.637105 2.725748
#> 3      s3 2.284481 2.691827 2.746810
```

## LOSO basis refitting on a fixed anchor representation

Contrasts are vectors over the anchor index. For illustration we select
the first anchor coordinate and refit the DKGE basis without each
held-out subject. The anchors themselves remain the pooled anchors
created above.

``` r
contrast_vec <- rep(0, fit_anchor$provenance$anchors$L)
contrast_vec[1] <- 1

res_contrast <- dkge_contrast(fit_anchor,
                               contrasts = list(anchor1 = contrast_vec),
                               method = "loso")
res_contrast$values$anchor1
#> $s1
#>  [1]  0.6661474804  0.0380264465  0.0672318690 -0.0164033904  0.0068233108
#>  [6]  0.0044237672  0.0026184187 -0.0040477714 -0.0043873945 -0.0011417468
#> [11] -0.0006777050 -0.0070584071  0.0005356340 -0.0002535858  0.0001932203
#> [16]  0.0021004627
#> 
#> $s2
#>  [1]  0.795526890 -0.092946075 -0.217101867  0.019278187  0.019240180
#>  [6] -0.006793277  0.002019734 -0.005389905  0.002999691 -0.001541060
#> [11]  0.002125122 -0.000996077 -0.001134298 -0.003550054  0.001032554
#> [16]  0.003181416
#> 
#> $s3
#>  [1] -0.0657289971  0.4226379838  0.1312078198  0.0457933474 -0.0535131315
#>  [6] -0.0032367720  0.0158027103  0.0251838360  0.0025116447 -0.0013197438
#> [11]  0.0014247376 -0.0087628490 -0.0045179459 -0.0002163051 -0.0017240536
#> [16]  0.0016514285
```

## Using the pipeline helper

[`dkge_pipeline()`](https://bbuchsbaum.github.io/dkge/reference/dkge_pipeline.md)
takes the same descriptor through its `input` argument, so the anchor
workflow runs end to end without assembling the stages by hand.

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
#> metadata  18     -none- list
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
#> classA 7.070029e-01 7.070029e-01 1.428454e-12 3.353575e-07 2.324685e-03
#> classB 3.538618e-13 3.352883e-07 7.068548e-01 7.068548e-01 3.464425e-13
#>                [,6]         [,7]         [,8]         [,9]        [,10]
#> classA 5.322435e-03 1.340088e-11 9.871332e-08 2.192997e-03 3.887870e-03
#> classB 2.231057e-08 3.966592e-03 2.379476e-02 1.560241e-11 2.463435e-06
#>               [,11]        [,12]        [,13]        [,14]        [,15]
#> classA 4.922622e-03 2.805816e-09 6.376889e-03 9.929082e-07 3.975926e-03
#> classB 5.235459e-11 2.631415e-03 3.209766e-11 1.111844e-02 5.092563e-12
#>               [,16]
#> classA 1.262806e-02
#> classB 1.584097e-07

# Ready for classification
cls <- dkge_classify(fit_anchor,
                     targets = target_matrix,
                     method = "lda",
                     folds = 2)
cls$summary
#> NULL
```

Here `folds = 2` cross-validates the classifier conditional on the
fitted anchor representation and DKGE fit. It does not make anchor
selection or basis fitting fold-specific, so treat this chunk as an API
example rather than an estimate of end-to-end generalization
performance.

## Diagnostics and provenance

Anchor coverage, leverage, and bandwidth settings are stored under
`fit_anchor$provenance$anchors`. Read them to see whether a few subjects
dominate the anchor basis.

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
#> 1      s1 2.343931 2.524823 2.591260
#> 2      s2 2.292074 2.637105 2.725748
#> 3      s3 2.284481 2.691827 2.746810
#> 
#> $leverage
#>       anchor  leverage
#> 1   anchor_1 1.3457778
#> 2   anchor_2 1.5800768
#> 3   anchor_3 1.1389105
#> 4   anchor_4 1.0412132
#> 5   anchor_5 0.6859295
#> 6   anchor_6 1.4390011
#> 7   anchor_7 0.6417655
#> 8   anchor_8 1.3194554
#> 9   anchor_9 1.1064288
#> 10 anchor_10 0.9256246
#> 11 anchor_11 0.7204953
#> 12 anchor_12 0.5705827
#> 13 anchor_13 0.7818146
#> 14 anchor_14 0.6982803
#> 15 anchor_15 0.9106550
#> 16 anchor_16 1.0939890
```

## Summary

The feature-anchored path aligns subjects with disjoint item sets
without imputing missing cells. It selects a shared anchor
representation, orthonormalizes that representation, and delegates the
aligned kernels to the core DKGE fitter. For confirmatory prediction,
define anchors independently or wrap anchor construction and model
fitting in an outer resampling loop; the convenience workflow shown here
fixes anchors using all subjects.

### Special cases: shared and mixed item sets

- **All subjects share the same items.** The anchor pipeline
  re-expresses the common item kernel because every subject projects
  identical feature rows into the same anchor coordinates. You can keep
  the anchor path for consistency or fit the shared item kernel with
  [`dkge_fit_from_kernels()`](https://bbuchsbaum.github.io/dkge/reference/dkge_fit_from_kernels.md);
  these are different representations and should be compared empirically
  if the choice matters.
- **Subgroups with identical items.** Subjects with identical feature
  rows receive projections into the same pooled anchor basis. Coverage
  and leverage diagnostics in `fit_anchor$provenance$anchors` help
  reveal uneven representation; large leverage spikes indicate anchors
  dominated by a subset and may motivate a smaller `L` or tighter
  bandwidth. These diagnostics do not substitute for nested
  preprocessing when the target is out-of-sample performance.

## Where to go next

- [`vignette("dkge-dense-rendering")`](https://bbuchsbaum.github.io/dkge/articles/dkge-dense-rendering.md)
  — the other half of the spatial story: once subjects are aligned
  through anchors, this is how cluster values become a dense field on a
  shared support.
- [`vignette("dkge-performance")`](https://bbuchsbaum.github.io/dkge/articles/dkge-performance.md)
  — mapper choice and warm starts, which is where anchor count and
  bandwidth start to cost you runtime.
- [`vignette("dkge-workflow")`](https://bbuchsbaum.github.io/dkge/articles/dkge-workflow.md)
  — the ordinary discrete-grid path, for contrast with the
  feature-anchored one.
