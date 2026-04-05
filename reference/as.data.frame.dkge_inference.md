# Convert DKGE inference results to a tidy data frame

Convert DKGE inference results to a tidy data frame

## Usage

``` r
# S3 method for class 'dkge_inference'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- x:

  A \`dkge_inference\` object

- row.names:

  NULL or a character vector giving the row names

- optional:

  Logical; if TRUE, setting row names is optional

- ...:

  Additional arguments passed to \[base::data.frame()\], including
  \`stringsAsFactors\`

## Value

Data frame with columns \`contrast\`, \`cluster\`, \`statistic\`,
\`p_value\`, \`p_adjusted\`, and \`significant\`
