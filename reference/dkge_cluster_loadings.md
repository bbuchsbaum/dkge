# Cluster-to-latent loadings for DKGE subjects

Computes \\A_s = B_s^\top K U\\, the linear map that pulls latent-space
vectors back to subject cluster space. Multiplying \\A_s\\ by a latent
coefficient vector produces a Haufe-style decoder at the subject level.

## Usage

``` r
dkge_cluster_loadings(fit)
```

## Arguments

- fit:

  Fitted object of class \`dkge\` containing \`Btil\`, \`K\`, and \`U\`.

## Value

A list of \`P_s x r\` matrices of cluster loadings.
