# Declare the provenance of an inferential transport operator

Transport learned from sign-sensitive subject loadings can change under
a sign flip. This constructor records why an operator may be held fixed,
or whether it is rebuilt for every randomization. The declaration is
checked by \[dkge_infer()\]; it does not turn an unsafe operator into a
valid one.

## Usage

``` r
dkge_transport_provenance(class, source_sha = NULL, details = list())
```

## Arguments

- class:

  One of \`"independent_training"\`, \`"prespecified_frozen"\`,
  \`"geometry_only"\`, \`"sign_invariant"\`, \`"fully_recomputed"\`, or
  \`"conditional_approximate"\`.

- source_sha:

  Optional immutable source or training-artifact SHA.

- details:

  Optional named list with supporting provenance details.

## Value

A versioned \`dkge_transport_provenance\` record.
