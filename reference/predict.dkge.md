# Predict contrasts for new subjects using a DKGE fit

S3 front-end that forwards to \[dkge_predict()\] while accepting
\`newdata\` lists with \`betas\`/\`B_list\` and \`contrasts\` entries.

## Usage

``` r
# S3 method for class 'dkge'
predict(object, newdata = NULL, ...)

# S3 method for class 'dkge_model'
predict(object, newdata = NULL, ...)
```

## Arguments

- object:

  A \`dkge\` fit.

- newdata:

  Optional list with elements \`betas\`/\`B_list\` and \`contrasts\`.

- ...:

  Additional arguments passed to \[dkge_predict()\].

## Value

Output of \[dkge_predict()\].
