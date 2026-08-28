# Construct a renderer for an identified reference support

A renderer never fits subject correspondence. It only attaches display
mechanics to an already identified support and consumes aligned maps or
an already computed statistic.

## Usage

``` r
dkge_renderer(support, provenance = NULL)
```

## Arguments

- support:

  A
  [`dkge_reference_support()`](https://bbuchsbaum.github.io/dkge/reference/dkge_reference_support.md)
  object.

- provenance:

  Optional rendering provenance.

## Value

An immutable `dkge_renderer` object.
