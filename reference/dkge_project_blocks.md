# Project new blocks into DKGE score space

Project new blocks into DKGE score space

## Usage

``` r
dkge_project_blocks(fit, B_list, Omega_list = NULL, w = NULL)
```

## Arguments

- fit:

  A \`dkge\` object.

- B_list:

  List of subject beta matrices.

- Omega_list:

  Optional list of spatial weights aligned with \`B_list\`.

- w:

  Optional vector of subject weights for the new data.

## Value

Matrix of projected scores (qxrank).
