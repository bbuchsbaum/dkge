# Apply a fitted mapper to new source values

Apply a fitted mapper to new source values

## Usage

``` r
predict_mapper(mapping, new_source_vals, ...)
```

## Arguments

- mapping:

  Mapping object returned by \[fit_mapper()\].

- new_source_vals:

  Matrix or vector of new source values (P_s x K).

- ...:

  Optional arguments used by specific strategies.

## Value

Matrix (Q x K) or vector of mapped values.
