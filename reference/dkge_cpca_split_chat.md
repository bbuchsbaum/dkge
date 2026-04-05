# Split compressed covariance into design/residual parts (CPCA inside-span)

Split compressed covariance into design/residual parts (CPCA
inside-span)

## Usage

``` r
dkge_cpca_split_chat(Chat, T, K)
```

## Arguments

- Chat:

  qxq compressed covariance expressed in the \\K^{1/2}\\ metric.

- T:

  qxq0 matrix describing the design-aligned subspace.

- K:

  qxq design kernel.

## Value

List with \`Chat_design\`, \`Chat_resid\`, and the projector \`P_hat\`.
