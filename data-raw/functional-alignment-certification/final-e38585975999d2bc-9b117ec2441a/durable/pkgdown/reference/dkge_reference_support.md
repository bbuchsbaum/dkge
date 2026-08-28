# Define the identified support on which aligned maps are represented

A reference support owns only anatomical/display information:
coordinates, labels or topology, and an optional decoder. Functional
features belong to
[`dkge_functional_template()`](https://bbuchsbaum.github.io/dkge/reference/dkge_functional_template.md)
and are deliberately excluded here. A bare MNI lattice is therefore a
valid support but is not functional correspondence.

## Usage

``` r
dkge_reference_support(
  coordinates,
  labels = NULL,
  topology = NULL,
  decoder = NULL,
  support_id = NULL,
  provenance = NULL
)
```

## Arguments

- coordinates:

  Finite numeric matrix with one row per target location.

- labels:

  Optional target labels in the same row order.

- topology:

  Optional topology or adjacency metadata.

- decoder:

  Optional pre-fitted support-to-display decoder.

- support_id:

  Optional stable identifier. A content-derived identifier is used by
  default.

- provenance:

  Optional list describing how the support was obtained.

## Value

An immutable `dkge_reference_support` object.
