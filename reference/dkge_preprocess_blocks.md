# Preprocess multiple blocks into the DKGE training space

Preprocess multiple blocks into the DKGE training space

## Usage

``` r
dkge_preprocess_blocks(fit, B_list, Omega_list = NULL, w = NULL)
```

## Arguments

- fit:

  A \`dkge\` object returned by \[dkge_fit()\].

- B_list:

  List of subject beta matrices.

- Omega_list:

  Optional list of spatial weights aligned with \`B_list\`.

- w:

  Optional vector of subject weights for the new data.

## Value

qx(sum P_s) matrix matching the training block layout.
