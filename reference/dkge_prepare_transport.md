# Prepare subject-to-medoid transport operators for reuse

Computes and caches subject-to-medoid transport matrices so downstream
routines (e.g. bootstraps) can reuse a fixed consensus mapping without
re-solving the transport problem on every call.

## Usage

``` r
dkge_prepare_transport(
  fit,
  centroids = NULL,
  loadings = NULL,
  betas = NULL,
  sizes = NULL,
  mapper = "sinkhorn",
  medoid = 1L,
  reference_selection = NULL,
  alignment_features = NULL,
  preprocessing = NULL,
  ...
)
```

## Arguments

- fit:

  A `dkge` object.

- centroids:

  List of subject centroid matrices. Defaults to the centroids stored on
  `fit` or `fit$input`.

- loadings:

  Optional list of subject loadings (`P_s x r`). When omitted, they are
  recomputed from `fit$Btil` or the supplied `betas`.

- betas:

  Optional list of subject betas used to recompute loadings when
  `loadings` is `NULL`.

- sizes:

  Optional list of cluster masses.

- mapper:

  Mapper specification or shorthand passed to
  [`dkge_mapper_spec()`](https://bbuchsbaum.github.io/dkge/reference/dkge_mapper_spec.md).

- medoid:

  Index (1-based) of the reference subject.

- reference_selection:

  Optional typed selection from
  [`dkge_select_reference_subject()`](https://bbuchsbaum.github.io/dkge/reference/dkge_select_reference_subject.md).
  When supplied, its selected subject is authoritative and `medoid` is
  only a compatibility alias.

- alignment_features:

  Optional typed feature object used both for reference selection and
  mapper fitting.

- preprocessing:

  Optional immutable provenance for feature construction.
  Caller-authored provenance is accepted only when it exactly matches a
  typed `alignment_features` object. Loose loadings/betas are always
  recorded as descriptive and inferentially ineligible.

- ...:

  Additional mapper arguments such as `epsilon` or `lambda_spa`.

## Value

A list containing cached application `operators`, joint transport
`plans`, solver `diagnostics`, `mapper_spec`, `feature_list`,
`size_list`, `feature_ref`, `size_ref`, `centroids`, and `medoid`.
