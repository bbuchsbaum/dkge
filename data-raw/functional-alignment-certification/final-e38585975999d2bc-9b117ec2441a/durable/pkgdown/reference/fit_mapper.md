# Fit a mapper on subject/reference features

Fit a mapper on subject/reference features

## Usage

``` r
fit_mapper(spec, ...)
```

## Arguments

- spec:

  Mapper specification created with
  [`dkge_mapper_spec()`](https://bbuchsbaum.github.io/dkge/reference/dkge_mapper_spec.md)
  or
  [`dkge_mapper()`](https://bbuchsbaum.github.io/dkge/reference/dkge_mapper.md).

- ...:

  Strategy-specific arguments (e.g. feature matrices, spatial
  coordinates, weights).

## Value

A fitted mapping object (see strategy-specific classes).
